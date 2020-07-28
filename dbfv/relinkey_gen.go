package dbfv

import (
	"errors"

	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// RKGProtocol is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type RKGProtocol struct {
	context                  *dbfvContext
	tmpPoly1                 *ring.Poly
	tmpPoly2                 *ring.Poly
	polypool                 *ring.Poly
	gaussianSampler          *ring.GaussianSampler
	ternarySamplerMontgomery *ring.TernarySampler
}

type RKGShare [][2]*ring.Poly

// MarshalBinary encodes the target element on a slice of bytes.
func (share *RKGShare) MarshalBinary() ([]byte, error) {
	//we have modulus * bitLog * Len of 1 ring rings
	rLength := ((*share)[0])[0].GetDataLen(true)
	data := make([]byte, 1+2*rLength*uint64(len(*share)))
	if len(*share) > 0xFF {
		return []byte{}, errors.New("RKGShare : uint8 overflow on length")
	}
	data[0] = uint8(len(*share))

	//write all of our rings in the data.
	//write all the polys
	ptr := uint64(1)
	for _, elem := range *share {
		_, err := elem[0].WriteTo(data[ptr : ptr+rLength])
		if err != nil {
			return []byte{}, err
		}
		ptr += rLength
		_, err = elem[1].WriteTo(data[ptr : ptr+rLength])
		if err != nil {
			return []byte{}, err
		}
		ptr += rLength
	}

	return data, nil

}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RKGShare) UnmarshalBinary(data []byte) error {
	lenShare := data[0]
	rLength := (len(data) - 1) / (2 * int(lenShare))

	if *share == nil {
		*share = make([][2]*ring.Poly, lenShare)
	}
	ptr := (1)
	for i := (0); i < int(lenShare); i++ {
		if (*share)[i][0] == nil || (*share)[i][1] == nil {
			(*share)[i][0] = new(ring.Poly)
			(*share)[i][1] = new(ring.Poly)
		}

		err := (*share)[i][0].UnmarshalBinary(data[ptr : ptr+rLength])
		if err != nil {
			return err
		}
		ptr += rLength
		err = (*share)[i][1].UnmarshalBinary(data[ptr : ptr+rLength])
		if err != nil {
			return err
		}
		ptr += rLength

	}

	return nil
}

// AllocateShares allocates the shares of the EKG protocol.
func (ekg *RKGProtocol) AllocateShares() (r1 RKGShare, r2 RKGShare) {
	r1 = make([][2]*ring.Poly, ekg.context.params.Beta())
	r2 = make([][2]*ring.Poly, ekg.context.params.Beta())
	for i := uint64(0); i < ekg.context.params.Beta(); i++ {
		r1[i][0] = ekg.context.contextQP.NewPoly()
		r1[i][1] = ekg.context.contextQP.NewPoly()
		r2[i][0] = ekg.context.contextQP.NewPoly()
		r2[i][1] = ekg.context.contextQP.NewPoly()
	}
	return
}

// NewEkgProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewEkgProtocol(params *bfv.Parameters) *RKGProtocol {

	context := newDbfvContext(params)

	ekg := new(RKGProtocol)
	ekg.context = context

	ekg.tmpPoly1 = ekg.context.contextQP.NewPoly()
	ekg.tmpPoly2 = ekg.context.contextQP.NewPoly()
	ekg.polypool = ekg.context.contextQP.NewPoly()
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ekg.ternarySamplerMontgomery = ring.NewTernarySampler(prng, ekg.context.contextQP, 0.5, true)
	ekg.gaussianSampler = ring.NewGaussianSampler(prng, ekg.context.contextQP, params.Sigma(), uint64(6*params.Sigma()))

	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first rounds).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *RKGProtocol) NewEphemeralKey() (ephemeralKey *ring.Poly) {
	ephemeralKey = ekg.ternarySamplerMontgomery.ReadNew()
	ekg.context.contextQP.NTT(ephemeralKey, ephemeralKey)
	return ephemeralKey
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(u, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShare) {

	var index uint64
	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	ringContext := ekg.context.contextQP

	ekg.polypool.Copy(sk)

	ringContext.MulScalarBigint(ekg.polypool, ekg.context.contextP.ModulusBigint, ekg.polypool)

	ringContext.InvMForm(ekg.polypool, ekg.polypool)

	for i := uint64(0); i < ekg.context.params.Beta(); i++ {

		// h = e
		ekg.gaussianSampler.Read(shareOut[i][0])
		ringContext.NTT(shareOut[i][0], shareOut[i][0])

		// h = sk*CrtBaseDecompQi + e
		for j := uint64(0); j < ekg.context.params.Alpha(); j++ {

			index = i*ekg.context.params.Alpha() + j
			qi := ringContext.Modulus[index]
			tmp0 := ekg.polypool.Coeffs[index]
			tmp1 := shareOut[i][0].Coeffs[index]

			for w := uint64(0); w < ekg.context.contextQP.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index == uint64(len(ekg.context.params.Qi())-1) {
				break
			}
		}
		// h = sk*CrtBaseDecompQi + -u*a + e
		ekg.context.contextQP.MulCoeffsMontgomeryAndSub(u, crp[i], shareOut[i][0])

		// Second Element
		// e_2i
		ekg.gaussianSampler.Read(shareOut[i][1])
		ringContext.NTT(shareOut[i][1], shareOut[i][1])
		// s*a + e_2i
		ringContext.MulCoeffsMontgomeryAndAdd(sk, crp[i], shareOut[i][1])
	}

	ekg.polypool.Zero() // TODO: check if we can remove this one

	return
}

// AggregateShareRoundOne adds share1 and share2 on shareOut.
func (ekg *RKGProtocol) AggregateShareRoundOne(share1, share2, shareOut RKGShare) {

	for i := uint64(0); i < ekg.context.params.Beta(); i++ {
		ekg.context.contextQP.Add(share1[i][0], share2[i][0], shareOut[i][0])
		ekg.context.contextQP.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}

}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Upon receiving the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(round1 RKGShare, u, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShare) {

	ringContext := ekg.context.contextQP

	// (u_i - s_i)
	ringContext.Sub(u, sk, ekg.tmpPoly1)

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := uint64(0); i < ekg.context.params.Beta(); i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// (AggregateShareRoundTwo samples) * sk
		ringContext.MulCoeffsMontgomery(round1[i][0], sk, shareOut[i][0])

		// (AggregateShareRoundTwo samples) * sk + e_1i
		ekg.gaussianSampler.Read(ekg.tmpPoly2)
		ringContext.NTT(ekg.tmpPoly2, ekg.tmpPoly2)
		ringContext.Add(shareOut[i][0], ekg.tmpPoly2, shareOut[i][0])

		// second part
		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		ekg.gaussianSampler.Read(shareOut[i][1])
		ringContext.NTT(shareOut[i][1], shareOut[i][1])
		ringContext.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1, round1[i][1], shareOut[i][1])
	}

}

// AggregateShareRoundTwo is the first part of the third and last round of the RKGProtocol protocol. Upon receiving the j-1 elements, each party
// computes :
//
// [sum(s_j * (-u*a + s*w + e) + e_j1), sum(s_j*a + e_j2)]
//
// = [s * (-u*a + s*w + e) + e_1, s*a + e_2].
func (ekg *RKGProtocol) AggregateShareRoundTwo(share1, share2, shareOut RKGShare) {

	for i := uint64(0); i < ekg.context.params.Beta(); i++ {
		ekg.context.contextQP.Add(share1[i][0], share2[i][0], shareOut[i][0])
		ekg.context.contextQP.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}
}

// GenRelinearizationKey finalizes the protocol and returns the common EvaluationKey.
func (ekg *RKGProtocol) GenRelinearizationKey(round1 RKGShare, round2 RKGShare, evalKeyOut *bfv.EvaluationKey) {

	ringContext := ekg.context.contextQP

	key := evalKeyOut.Get()[0].Get()
	for i := uint64(0); i < ekg.context.params.Beta(); i++ {

		ringContext.Add(round2[i][0], round2[i][1], key[i][0])
		key[i][1].Copy(round1[i][1])

		ringContext.MForm(key[i][0], key[i][0])
		ringContext.MForm(key[i][1], key[i][1])

	}
}
