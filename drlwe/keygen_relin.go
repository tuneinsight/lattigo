package drlwe

import (
	"errors"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// RelinearizationKeyGenerator is an interface describing the local steps of a generic RLWE RKG protocol
type RelinearizationKeyGenerator interface {
	AllocateShares() (ephKey *rlwe.SecretKey, r1 *RKGShare, r2 *RKGShare, rkgCRP RKGCRP)
	GenShareRoundOne(sk *rlwe.SecretKey, crp RKGCRP, ephKeyOut *rlwe.SecretKey, shareOut *RKGShare)
	GenShareRoundTwo(ephSk, sk *rlwe.SecretKey, round1 *RKGShare, crp RKGCRP, shareOut *RKGShare)
	AggregateShares(share1, share2, shareOut *RKGShare)
	GenRelinearizationKey(round1 *RKGShare, round2 *RKGShare, relinKeyOut *rlwe.RelinearizationKey)
}

// RKGProtocol is the structure storing the parameters and and precomputations for the collective relinearization key generation protocol.
type RKGProtocol struct {
	params           rlwe.Parameters
	pBigInt          *big.Int
	gaussianSamplerQ *ring.GaussianSampler
	ternarySamplerQ  *ring.TernarySampler // sampling in Montgomerry form

	tmpPoly1 ring.PolyQP
	tmpPoly2 ring.PolyQP
}

// RKGShare is a share in the RKG protocol
type RKGShare struct {
	Value [][2]ring.PolyQP
}

// NewRKGProtocol creates a new RKG protocol struct
func NewRKGProtocol(params rlwe.Parameters, ephSkPr float64) *RKGProtocol {
	rkg := new(RKGProtocol)
	rkg.params = params

	var err error
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	rkg.pBigInt = params.PBigInt()
	rkg.gaussianSamplerQ = ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))
	rkg.ternarySamplerQ = ring.NewTernarySampler(prng, params.RingQ(), ephSkPr, false)
	rkg.tmpPoly1 = *params.RingQP().NewPoly()
	rkg.tmpPoly2 = *params.RingQP().NewPoly()
	return rkg
}

// AllocateShares allocates the shares of the EKG protocol.
func (ekg *RKGProtocol) AllocateShares() (ephSk *rlwe.SecretKey, r1 *RKGShare, r2 *RKGShare, rkgCRP RKGCRP) {
	ephSk = rlwe.NewSecretKey(ekg.params)
	r1, r2 = new(RKGShare), new(RKGShare)
	r1.Value = make([][2]ring.PolyQP, ekg.params.Beta())
	r2.Value = make([][2]ring.PolyQP, ekg.params.Beta())
	rkgCRP = make([]ring.PolyQP, ekg.params.Beta())
	for i := 0; i < ekg.params.Beta(); i++ {
		r1.Value[i][0] = *ekg.params.RingQP().NewPoly()
		r1.Value[i][1] = *ekg.params.RingQP().NewPoly()
		r2.Value[i][0] = *ekg.params.RingQP().NewPoly()
		r2.Value[i][1] = *ekg.params.RingQP().NewPoly()
		rkgCRP[i] = *ekg.params.RingQP().NewPoly()
	}
	return
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(sk *rlwe.SecretKey, crp RKGCRP, ephSkOut *rlwe.SecretKey, shareOut *RKGShare) {
	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	ringQ := ekg.params.RingQ()
	ringQP := ekg.params.RingQP()
	levelP := ekg.params.PCount() - 1

	ringQ.MulScalarBigint(sk.Value.Q, ekg.pBigInt, ekg.tmpPoly1.Q)
	ringQ.InvMForm(ekg.tmpPoly1.Q, ekg.tmpPoly1.Q)

	ekg.ternarySamplerQ.Read(ephSkOut.Value.Q)
	ringQP.ExtendBasisSmallNormAndCenter(ephSkOut.Value.Q, levelP, &ephSkOut.Value)
	ringQP.NTT(&ephSkOut.Value, &ephSkOut.Value)
	ringQP.MForm(&ephSkOut.Value, &ephSkOut.Value)

	for i := 0; i < ekg.params.Beta(); i++ {
		// h = e
		ekg.gaussianSamplerQ.Read(shareOut.Value[i][0].Q)
		ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][0].Q, levelP, &shareOut.Value[i][0])
		ringQP.NTT(&shareOut.Value[i][0], &shareOut.Value[i][0])

		// h = sk*CrtBaseDecompQi + e
		for j := 0; j < ekg.params.PCount(); j++ {
			index := i*ekg.params.PCount() + j

			// Handles the case where nb pj does not divides nb qi
			if index >= ekg.params.QCount() {
				break
			}

			qi := ringQ.Modulus[index]
			skP := ekg.tmpPoly1.Q.Coeffs[index]
			h := shareOut.Value[i][0].Q.Coeffs[index]

			for w := 0; w < ringQ.N; w++ {
				h[w] = ring.CRed(h[w]+skP[w], qi)
			}
		}

		// h = sk*CrtBaseDecompQi + -u*a + e
		ringQP.MulCoeffsMontgomeryAndSub(&ephSkOut.Value, &crp[i], &shareOut.Value[i][0])

		// Second Element
		// e_2i
		ekg.gaussianSamplerQ.Read(shareOut.Value[i][1].Q)
		ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][1].Q, levelP, &shareOut.Value[i][1])
		ringQP.NTT(&shareOut.Value[i][1], &shareOut.Value[i][1])
		// s*a + e_2i
		ringQP.MulCoeffsMontgomeryAndAdd(&sk.Value, &crp[i], &shareOut.Value[i][1])
	}
}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Upon receiving the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(ephSk, sk *rlwe.SecretKey, round1 *RKGShare, crp RKGCRP, shareOut *RKGShare) {

	ringQP := ekg.params.RingQP()
	levelP := ekg.params.PCount() - 1

	// (u_i - s_i)
	ringQP.Sub(&ephSk.Value, &sk.Value, &ekg.tmpPoly1)

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := 0; i < ekg.params.Beta(); i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// (AggregateShareRoundTwo samples) * sk
		ringQP.MulCoeffsMontgomeryConstant(&round1.Value[i][0], &sk.Value, &shareOut.Value[i][0])

		// (AggregateShareRoundTwo samples) * sk + e_1i
		ekg.gaussianSamplerQ.Read(ekg.tmpPoly2.Q)
		ringQP.ExtendBasisSmallNormAndCenter(ekg.tmpPoly2.Q, levelP, &ekg.tmpPoly2)
		ringQP.NTT(&ekg.tmpPoly2, &ekg.tmpPoly2)
		ringQP.Add(&shareOut.Value[i][0], &ekg.tmpPoly2, &shareOut.Value[i][0])

		// second part
		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		ekg.gaussianSamplerQ.Read(shareOut.Value[i][1].Q)
		ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][1].Q, levelP, &shareOut.Value[i][1])
		ringQP.NTT(&shareOut.Value[i][1], &shareOut.Value[i][1])
		ringQP.MulCoeffsMontgomeryAndAdd(&ekg.tmpPoly1, &round1.Value[i][1], &shareOut.Value[i][1])
	}

}

// AggregateShares combines two RKG shares into a single one
func (ekg *RKGProtocol) AggregateShares(share1, share2, shareOut *RKGShare) {

	for i := 0; i < ekg.params.Beta(); i++ {
		ekg.params.RingQP().Add(&share1.Value[i][0], &share2.Value[i][0], &shareOut.Value[i][0])
		ekg.params.RingQP().Add(&share1.Value[i][1], &share2.Value[i][1], &shareOut.Value[i][1])
	}
}

// GenRelinearizationKey computes the generated RLK from the public shares and write the result in evalKeyOut
func (ekg *RKGProtocol) GenRelinearizationKey(round1 *RKGShare, round2 *RKGShare, evalKeyOut *rlwe.RelinearizationKey) {
	for i := 0; i < ekg.params.Beta(); i++ {
		ekg.params.RingQP().Add(&round2.Value[i][0], &round2.Value[i][1], &evalKeyOut.Keys[0].Value[i][0])
		evalKeyOut.Keys[0].Value[i][1].Copy(&round1.Value[i][1])

		ekg.params.RingQP().MForm(&evalKeyOut.Keys[0].Value[i][0], &evalKeyOut.Keys[0].Value[i][0])
		ekg.params.RingQP().MForm(&evalKeyOut.Keys[0].Value[i][1], &evalKeyOut.Keys[0].Value[i][1])
	}
}

// MarshalBinary encodes the target element on a slice of bytes.
func (share *RKGShare) MarshalBinary() ([]byte, error) {
	//we have modulus * bitLog * Len of 1 ring rings
	data := make([]byte, 1+2*share.Value[0][0].GetDataLen(true)*len(share.Value))
	if len(share.Value) > 0xFF {
		return []byte{}, errors.New("RKGShare : uint8 overflow on length")
	}
	data[0] = uint8(len(share.Value))

	//write all of our rings in the data.
	//write all the polys
	ptr := 1
	var inc int
	var err error
	for _, elem := range share.Value {

		if inc, err = elem[0].WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc

		if inc, err = elem[1].WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc
	}

	return data, nil

}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RKGShare) UnmarshalBinary(data []byte) (err error) {
	share.Value = make([][2]ring.PolyQP, data[0])
	ptr := 1
	var inc int
	for i := range share.Value {
		if inc, err = share.Value[i][0].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}
		ptr += inc

		if inc, err = share.Value[i][1].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}
		ptr += inc
	}

	return nil
}
