package drlwe

import (
	"errors"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math/big"
)

// RelinearizationKeyGenerator is an interface describing the local steps of a generic RLWE RKG protocol
type RelinearizationKeyGenerator interface {
	AllocateShares() (ephKey *rlwe.SecretKey, r1 *RKGShare, r2 *RKGShare)
	GenShareRoundOne(sk *rlwe.SecretKey, crp []*ring.Poly, ephKeyOut *rlwe.SecretKey, shareOut *RKGShare)
	GenShareRoundTwo(ephSk, sk *rlwe.SecretKey, round1 *RKGShare, crp []*ring.Poly, shareOut *RKGShare)
	AggregateShares(share1, share2, shareOut *RKGShare)
	GenRelinearizationKey(round1 *RKGShare, round2 *RKGShare, relinKeyOut *rlwe.RelinearizationKey) // TODO type for generic eval key
}

// RKGProtocol is the structure storing the parameters and and precomputations for the collective relinearization key generation protocol.
type RKGProtocol struct {
	params           rlwe.Parameters
	pBigInt          *big.Int
	ringQ            *ring.Ring
	ringP            *ring.Ring
	gaussianSamplerQ *ring.GaussianSampler
	ternarySamplerQ  *ring.TernarySampler // sampling in Montgomerry form

	tmpPoly1 [2]*ring.Poly
	tmpPoly2 [2]*ring.Poly
}

// RKGShare is a share in the RKG protocol
type RKGShare struct {
	Value [][2][2]*ring.Poly
}

// NewRKGProtocol creates a new RKG protocol struct
func NewRKGProtocol(params rlwe.Parameters, ephSkPr float64) *RKGProtocol {
	rkg := new(RKGProtocol)
	rkg.params = params
	rkg.ringQ = params.RingQ()
	rkg.ringP = params.RingP()

	var err error
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err) // TODO error
	}

	rkg.pBigInt = params.PBigInt()
	rkg.gaussianSamplerQ = ring.NewGaussianSampler(prng, rkg.ringQ, params.Sigma(), int(6*params.Sigma()))
	rkg.ternarySamplerQ = ring.NewTernarySampler(prng, rkg.ringQ, ephSkPr, false)
	rkg.tmpPoly1 = [2]*ring.Poly{rkg.ringQ.NewPoly(), rkg.ringP.NewPoly()}
	rkg.tmpPoly2 = [2]*ring.Poly{rkg.ringQ.NewPoly(), rkg.ringP.NewPoly()}
	return rkg
}

// AllocateShares allocates the shares of the EKG protocol.
func (ekg *RKGProtocol) AllocateShares() (ephSk *rlwe.SecretKey, r1 *RKGShare, r2 *RKGShare) {
	ephSk = rlwe.NewSecretKey(ekg.params)
	r1, r2 = new(RKGShare), new(RKGShare)
	r1.Value = make([][2][2]*ring.Poly, ekg.params.Beta())
	r2.Value = make([][2][2]*ring.Poly, ekg.params.Beta())
	for i := 0; i < ekg.params.Beta(); i++ {
		r1.Value[i][0][0] = ekg.ringQ.NewPoly()
		r1.Value[i][0][1] = ekg.ringP.NewPoly()
		r1.Value[i][1][0] = ekg.ringQ.NewPoly()
		r1.Value[i][1][1] = ekg.ringP.NewPoly()

		r2.Value[i][0][0] = ekg.ringQ.NewPoly()
		r2.Value[i][0][1] = ekg.ringP.NewPoly()
		r2.Value[i][1][0] = ekg.ringQ.NewPoly()
		r2.Value[i][1][1] = ekg.ringP.NewPoly()
	}
	return
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(sk *rlwe.SecretKey, crp [][2]*ring.Poly, ephSkOut *rlwe.SecretKey, shareOut *RKGShare) {
	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	ringQ := ekg.ringQ
	ringP := ekg.ringP

	ringQ.MulScalarBigint(sk.Value[0], ekg.pBigInt, ekg.tmpPoly1[0])
	ringQ.InvMForm(ekg.tmpPoly1[0], ekg.tmpPoly1[0])

	ekg.ternarySamplerQ.Read(ephSkOut.Value[0])
	extendBasisSmallNormAndCenter(ringQ, ringP, ephSkOut.Value[0], ephSkOut.Value[1])
	ringQ.NTT(ephSkOut.Value[0], ephSkOut.Value[0])
	ringP.NTT(ephSkOut.Value[1], ephSkOut.Value[1])
	ringQ.MForm(ephSkOut.Value[0], ephSkOut.Value[0])
	ringP.MForm(ephSkOut.Value[1], ephSkOut.Value[1])

	for i := 0; i < ekg.params.Beta(); i++ {
		// h = e
		ekg.gaussianSamplerQ.Read(shareOut.Value[i][0][0])
		extendBasisSmallNormAndCenter(ringQ, ringP, shareOut.Value[i][0][0], shareOut.Value[i][0][1])
		ringQ.NTT(shareOut.Value[i][0][0], shareOut.Value[i][0][0])
		ringP.NTT(shareOut.Value[i][0][1], shareOut.Value[i][0][1])

		// h = sk*CrtBaseDecompQi + e
		for j := 0; j < ekg.params.PCount(); j++ {
			index := i*ekg.params.PCount() + j

			// Handles the case where nb pj does not divides nb qi
			if index >= ekg.params.QCount() {
				break
			}

			qi := ringQ.Modulus[index]
			skP := ekg.tmpPoly1[0].Coeffs[index]
			h := shareOut.Value[i][0][0].Coeffs[index]

			for w := 0; w < ringQ.N; w++ {
				h[w] = ring.CRed(h[w]+skP[w], qi)
			}
		}

		// h = sk*CrtBaseDecompQi + -u*a + e
		ringQ.MulCoeffsMontgomeryAndSub(ephSkOut.Value[0], crp[i][0], shareOut.Value[i][0][0])
		ringP.MulCoeffsMontgomeryAndSub(ephSkOut.Value[1], crp[i][1], shareOut.Value[i][0][1])

		// Second Element
		// e_2i
		ekg.gaussianSamplerQ.Read(shareOut.Value[i][1][0])
		extendBasisSmallNormAndCenter(ringQ, ringP, shareOut.Value[i][1][0], shareOut.Value[i][1][1])
		ringQ.NTT(shareOut.Value[i][1][0], shareOut.Value[i][1][0])
		ringP.NTT(shareOut.Value[i][1][1], shareOut.Value[i][1][1])
		// s*a + e_2i
		ringQ.MulCoeffsMontgomeryAndAdd(sk.Value[0], crp[i][0], shareOut.Value[i][1][0])
		ringP.MulCoeffsMontgomeryAndAdd(sk.Value[1], crp[i][1], shareOut.Value[i][1][1])
	}
}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Upon receiving the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(ephSk, sk *rlwe.SecretKey, round1 *RKGShare, crp [][2]*ring.Poly, shareOut *RKGShare) {

	ringQ := ekg.ringQ
	ringP := ekg.ringP

	// (u_i - s_i)
	ringQ.Sub(ephSk.Value[0], sk.Value[0], ekg.tmpPoly1[0])
	ringP.Sub(ephSk.Value[1], sk.Value[1], ekg.tmpPoly1[1])

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := 0; i < ekg.params.Beta(); i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// (AggregateShareRoundTwo samples) * sk
		ekg.ringQ.MulCoeffsMontgomeryConstant(round1.Value[i][0][0], sk.Value[0], shareOut.Value[i][0][0])
		ekg.ringP.MulCoeffsMontgomeryConstant(round1.Value[i][0][1], sk.Value[1], shareOut.Value[i][0][1])

		// (AggregateShareRoundTwo samples) * sk + e_1i
		ekg.gaussianSamplerQ.Read(ekg.tmpPoly2[0])
		extendBasisSmallNormAndCenter(ringQ, ringP, ekg.tmpPoly2[0], ekg.tmpPoly2[1])
		ringQ.NTT(ekg.tmpPoly2[0], ekg.tmpPoly2[0])
		ringP.NTT(ekg.tmpPoly2[1], ekg.tmpPoly2[1])
		ringQ.Add(shareOut.Value[i][0][0], ekg.tmpPoly2[0], shareOut.Value[i][0][0])
		ringP.Add(shareOut.Value[i][0][1], ekg.tmpPoly2[1], shareOut.Value[i][0][1])

		// second part
		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		ekg.gaussianSamplerQ.Read(shareOut.Value[i][1][0])
		extendBasisSmallNormAndCenter(ringQ, ringP, shareOut.Value[i][1][0], shareOut.Value[i][1][1])
		ringQ.NTT(shareOut.Value[i][1][0], shareOut.Value[i][1][0])
		ringP.NTT(shareOut.Value[i][1][1], shareOut.Value[i][1][1])
		ringQ.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1[0], round1.Value[i][1][0], shareOut.Value[i][1][0])
		ringP.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1[1], round1.Value[i][1][1], shareOut.Value[i][1][1])
	}

}

// AggregateShares combines two RKG shares into a single one
func (ekg *RKGProtocol) AggregateShares(share1, share2, shareOut *RKGShare) {

	for i := 0; i < ekg.params.Beta(); i++ {
		ekg.ringQ.Add(share1.Value[i][0][0], share2.Value[i][0][0], shareOut.Value[i][0][0])
		ekg.ringP.Add(share1.Value[i][0][1], share2.Value[i][0][1], shareOut.Value[i][0][1])
		ekg.ringQ.Add(share1.Value[i][1][0], share2.Value[i][1][0], shareOut.Value[i][1][0])
		ekg.ringP.Add(share1.Value[i][1][1], share2.Value[i][1][1], shareOut.Value[i][1][1])
	}
}

// GenRelinearizationKey computes the generated RLK from the public shares and write the result in evalKeyOut
func (ekg *RKGProtocol) GenRelinearizationKey(round1 *RKGShare, round2 *RKGShare, evalKeyOut *rlwe.RelinearizationKey) {
	for i := 0; i < ekg.params.Beta(); i++ {
		ekg.ringQ.Add(round2.Value[i][0][0], round2.Value[i][1][0], evalKeyOut.Keys[0].Value[i][0][0])
		ekg.ringP.Add(round2.Value[i][0][1], round2.Value[i][1][1], evalKeyOut.Keys[0].Value[i][0][1])
		evalKeyOut.Keys[0].Value[i][1][0].Copy(round1.Value[i][1][0])
		evalKeyOut.Keys[0].Value[i][1][1].Copy(round1.Value[i][1][1])

		ekg.ringQ.MForm(evalKeyOut.Keys[0].Value[i][0][0], evalKeyOut.Keys[0].Value[i][0][0])
		ekg.ringP.MForm(evalKeyOut.Keys[0].Value[i][0][1], evalKeyOut.Keys[0].Value[i][0][1])
		ekg.ringQ.MForm(evalKeyOut.Keys[0].Value[i][1][0], evalKeyOut.Keys[0].Value[i][1][0])
		ekg.ringP.MForm(evalKeyOut.Keys[0].Value[i][1][1], evalKeyOut.Keys[0].Value[i][1][1])
	}
}

// MarshalBinary encodes the target element on a slice of bytes.
func (share *RKGShare) MarshalBinary() ([]byte, error) {
	//we have modulus * bitLog * Len of 1 ring rings
	data := make([]byte, 1+2*(share.Value[0][0][0].GetDataLen(true)+share.Value[0][0][1].GetDataLen(true))*len(share.Value))
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

		if inc, err = elem[0][0].WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc

		if inc, err = elem[0][1].WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc

		if inc, err = elem[1][0].WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc

		if inc, err = elem[1][1].WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc
	}

	return data, nil

}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RKGShare) UnmarshalBinary(data []byte) (err error) {
	share.Value = make([][2][2]*ring.Poly, data[0])
	ptr := 1
	var inc int
	for i := range share.Value {
		share.Value[i][0][0] = new(ring.Poly)
		if inc, err = share.Value[i][0][0].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}

		ptr += inc

		share.Value[i][0][1] = new(ring.Poly)
		if inc, err = share.Value[i][0][1].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}

		ptr += inc

		share.Value[i][1][0] = new(ring.Poly)
		if inc, err = share.Value[i][1][0].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}

		ptr += inc

		share.Value[i][1][1] = new(ring.Poly)
		if inc, err = share.Value[i][1][1].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}

		ptr += inc
	}

	return nil
}
