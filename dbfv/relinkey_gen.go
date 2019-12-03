package dbfv

import (
	"errors"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// RKGProtocol is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type RKGProtocol struct {
	context  *dbfvContext
	tmpPoly1 *ring.Poly
	tmpPoly2 *ring.Poly
	polypool *ring.Poly
}

// RKGShareRoundOne is a struct storing the round one RKG shares.
type RKGShareRoundOne []*ring.Poly

// RKGShareRoundTwo is a struct storing the round two RKG shares.
type RKGShareRoundTwo [][2]*ring.Poly

// RKGShareRoundThree is a struct storing the round three RKG shares.
type RKGShareRoundThree []*ring.Poly

// MarshalBinary encodes the target element on a slice of bytes.
func (share *RKGShareRoundOne) MarshalBinary() ([]byte, error) {
	rLength := (*share)[0].GetDataLen(true)
	data := make([]byte, 1+rLength*uint64(len(*share)))
	data[0] = uint8(len(*share))

	pointer := uint64(1)
	for _, s := range *share {
		tmp, err := s.WriteTo(data[pointer : pointer+rLength])
		if err != nil {
			return []byte{}, err
		}
		pointer += tmp
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RKGShareRoundOne) UnmarshalBinary(data []byte) error {
	//share.modulus = data[0]
	lenShare := data[0]
	rLength := len(data[1:]) / int(lenShare)
	if *share == nil {
		*share = make([]*ring.Poly, lenShare)
	}
	ptr := 1
	for i := uint8(0); i < lenShare; i++ {
		if (*share)[i] == nil {
			(*share)[i] = new(ring.Poly)
		}
		err := (*share)[i].UnmarshalBinary(data[ptr : ptr+rLength])
		if err != nil {
			return err
		}
		ptr += rLength
	}

	return nil
}

// MarshalBinary encodes the target element on a slice of bytes.
func (share *RKGShareRoundTwo) MarshalBinary() ([]byte, error) {
	//we have modulus * bitLog * Len of 1 ring rings
	rLength := ((*share)[0])[0].GetDataLen(true)
	data := make([]byte, 1+2*rLength*uint64(len(*share)))
	if len(*share) > 0xFF {
		return []byte{}, errors.New("RKGShareRoundTwo : uint8 overflow on length")
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
func (share *RKGShareRoundTwo) UnmarshalBinary(data []byte) error {
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

// MarshalBinary encodes the target element on a slice of bytes.
func (share *RKGShareRoundThree) MarshalBinary() ([]byte, error) {
	rLength := (*share)[0].GetDataLen(true)
	data := make([]byte, 1+rLength*uint64(len(*share)))
	data[0] = uint8(len(*share))

	pointer := uint64(1)
	for _, s := range *share {
		tmp, err := s.WriteTo(data[pointer : pointer+rLength])
		if err != nil {
			return []byte{}, err
		}
		pointer += tmp
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RKGShareRoundThree) UnmarshalBinary(data []byte) error {
	//share.modulus = data[0]
	lenShare := data[0]
	rLength := len(data[1:]) / int(lenShare)
	if *share == nil {
		*share = make([]*ring.Poly, lenShare)
	}
	ptr := 1
	for i := uint8(0); i < lenShare; i++ {
		if (*share)[i] == nil {
			(*share)[i] = new(ring.Poly)
		}
		err := (*share)[i].UnmarshalBinary(data[ptr : ptr+rLength])
		if err != nil {
			return err
		}
		ptr += rLength
	}

	return nil
}

// AllocateShares allocates the shares of the EKG protocol.
func (ekg *RKGProtocol) AllocateShares() (r1 RKGShareRoundOne, r2 RKGShareRoundTwo, r3 RKGShareRoundThree) {
	r1 = make([]*ring.Poly, ekg.context.params.Beta())
	r2 = make([][2]*ring.Poly, ekg.context.params.Beta())
	r3 = make([]*ring.Poly, ekg.context.params.Beta())
	for i := uint64(0); i < ekg.context.params.Beta(); i++ {
		r1[i] = ekg.context.contextQP.NewPoly()
		r2[i][0] = ekg.context.contextQP.NewPoly()
		r2[i][1] = ekg.context.contextQP.NewPoly()
		r3[i] = ekg.context.contextQP.NewPoly()
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

	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first rounds).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *RKGProtocol) NewEphemeralKey(p float64) (ephemeralKey *ring.Poly) {
	return ekg.context.contextQP.SampleTernaryMontgomeryNTTNew(p)
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(u, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShareRoundOne) {

	var index uint64

	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	ekg.polypool.Copy(sk)

	ekg.context.contextQP.MulScalarBigint(ekg.polypool, ekg.context.contextP.ModulusBigint, ekg.polypool)

	ekg.context.contextQP.InvMForm(ekg.polypool, ekg.polypool)

	for i := uint64(0); i < ekg.context.params.Beta(); i++ {

		// h = e
		ekg.context.gaussianSampler.SampleNTT(shareOut[i])

		// h = sk*CrtBaseDecompQi + e
		for j := uint64(0); j < ekg.context.params.Alpha(); j++ {

			index = i*ekg.context.params.Alpha() + j
			qi := ekg.context.contextQP.Modulus[index]
			tmp0 := ekg.polypool.Coeffs[index]
			tmp1 := shareOut[i].Coeffs[index]

			for w := uint64(0); w < ekg.context.contextQP.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index == uint64(len(ekg.context.params.LogQi)-1) {
				break
			}
		}

		// h = sk*CrtBaseDecompQi + -u*a + e
		ekg.context.contextQP.MulCoeffsMontgomeryAndSub(u, crp[i], shareOut[i])
	}

	ekg.polypool.Zero()

	return
}

// AggregateShareRoundOne adds share1 and share2 on shareOut.
func (ekg *RKGProtocol) AggregateShareRoundOne(share1, share2, shareOut RKGShareRoundOne) {

	for i := uint64(0); i < ekg.context.params.Beta(); i++ {
		ekg.context.contextQP.Add(share1[i], share2[i], shareOut[i])
	}

}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Upon receiving the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(round1 RKGShareRoundOne, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShareRoundTwo) {

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := uint64(0); i < ekg.context.params.Beta(); i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// (AggregateShareRoundTwo samples) * sk
		ekg.context.contextQP.MulCoeffsMontgomery(round1[i], sk, shareOut[i][0])

		// (AggregateShareRoundTwo samples) * sk + e_1i
		ekg.context.gaussianSampler.SampleNTT(ekg.tmpPoly1)
		ekg.context.contextQP.Add(shareOut[i][0], ekg.tmpPoly1, shareOut[i][0])

		// Second Element
		// e_2i
		ekg.context.gaussianSampler.SampleNTT(shareOut[i][1])
		// s*a + e_2i
		ekg.context.contextQP.MulCoeffsMontgomeryAndAdd(sk, crp[i], shareOut[i][1])
	}

}

// AggregateShareRoundTwo is the first part of the third and last round of the RKGProtocol protocol. Upon receiving the j-1 elements, each party
// computes :
//
// [sum(s_j * (-u*a + s*w + e) + e_j1), sum(s_j*a + e_j2)]
//
// = [s * (-u*a + s*w + e) + e_1, s*a + e_2].
func (ekg *RKGProtocol) AggregateShareRoundTwo(share1, share2, shareOut RKGShareRoundTwo) {

	for i := uint64(0); i < ekg.context.params.Beta(); i++ {
		ekg.context.contextQP.Add(share1[i][0], share2[i][0], shareOut[i][0])
		ekg.context.contextQP.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}

}

// GenShareRoundThree is the second pard of the third and last round of the RKGProtocol protocol. Each party operates a key-switch on [s*a + e_2],
// by computing :
//
// [(u_i - s_i)*(s*a + e_2)]
//
// and broadcasts the result to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundThree(round2 RKGShareRoundTwo, u, sk *ring.Poly, shareOut RKGShareRoundThree) {

	// (u_i - s_i)
	ekg.context.contextQP.Sub(u, sk, ekg.tmpPoly1)

	for i := uint64(0); i < ekg.context.params.Beta(); i++ {

		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		ekg.context.gaussianSampler.SampleNTT(shareOut[i])
		ekg.context.contextQP.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1, round2[i][1], shareOut[i])
	}
}

// AggregateShareRoundThree adds share1 and share2 on shareOut.
func (ekg *RKGProtocol) AggregateShareRoundThree(share1, share2, shareOut RKGShareRoundThree) {
	for i := uint64(0); i < ekg.context.params.Beta(); i++ {
		ekg.context.contextQP.Add(share1[i], share2[i], shareOut[i])
	}
}

// GenRelinearizationKey finalizes the protocol and returns the common EvaluationKey.
func (ekg *RKGProtocol) GenRelinearizationKey(round2 RKGShareRoundTwo, round3 RKGShareRoundThree, evalKeyOut *bfv.EvaluationKey) {

	key := evalKeyOut.Get()[0].Get()
	for i := uint64(0); i < ekg.context.params.Beta(); i++ {

		ekg.context.contextQP.Add(round2[i][0], round3[i], key[i][0])
		key[i][1].Copy(round2[i][1])

		ekg.context.contextQP.MForm(key[i][0], key[i][0])
		ekg.context.contextQP.MForm(key[i][1], key[i][1])

	}
}
