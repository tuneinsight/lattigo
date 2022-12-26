package drlwe

import (
	"errors"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// RKGProtocol is the structure storing the parameters and and precomputations for the collective relinearization key generation protocol.
type RKGProtocol struct {
	params rlwe.Parameters

	gaussianSamplerQ *ring.GaussianSampler
	ternarySamplerQ  *ring.TernarySampler // sampling in Montgomery form

	tmpPoly1 ringqp.Poly
	tmpPoly2 ringqp.Poly
}

// ShallowCopy creates a shallow copy of RKGProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RKGProtocol can be used concurrently.
func (ekg *RKGProtocol) ShallowCopy() *RKGProtocol {
	var err error
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := ekg.params

	return &RKGProtocol{
		params:           ekg.params,
		gaussianSamplerQ: ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
		ternarySamplerQ:  ring.NewTernarySamplerWithHammingWeight(prng, params.RingQ(), params.HammingWeight(), false),
		tmpPoly1:         params.RingQP().NewPoly(),
		tmpPoly2:         params.RingQP().NewPoly(),
	}
}

// RKGShare is a share in the RKG protocol.
type RKGShare struct {
	Value [][][2]ringqp.Poly
}

// RKGCRP is a type for common reference polynomials in the RKG protocol.
type RKGCRP [][]ringqp.Poly

// NewRKGProtocol creates a new RKG protocol struct.
func NewRKGProtocol(params rlwe.Parameters) *RKGProtocol {
	rkg := new(RKGProtocol)
	rkg.params = params

	var err error
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	rkg.gaussianSamplerQ = ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))
	rkg.ternarySamplerQ = ring.NewTernarySamplerWithHammingWeight(prng, params.RingQ(), params.HammingWeight(), false)
	rkg.tmpPoly1 = params.RingQP().NewPoly()
	rkg.tmpPoly2 = params.RingQP().NewPoly()
	return rkg
}

// AllocateShare allocates the share of the EKG protocol.
func (ekg *RKGProtocol) AllocateShare() (ephSk *rlwe.SecretKey, r1 *RKGShare, r2 *RKGShare) {
	params := ekg.params
	ephSk = rlwe.NewSecretKey(params)
	r1, r2 = new(RKGShare), new(RKGShare)

	decompRNS := params.DecompRNS(params.QCount()-1, params.PCount()-1)
	decompPw2 := params.DecompPw2(params.QCount()-1, params.PCount()-1)

	r1.Value = make([][][2]ringqp.Poly, decompRNS)
	r2.Value = make([][][2]ringqp.Poly, decompRNS)

	for i := 0; i < decompRNS; i++ {
		r1.Value[i] = make([][2]ringqp.Poly, decompPw2)
		r2.Value[i] = make([][2]ringqp.Poly, decompPw2)
		for j := 0; j < decompPw2; j++ {
			r1.Value[i][j][0] = ekg.params.RingQP().NewPoly()
			r1.Value[i][j][1] = ekg.params.RingQP().NewPoly()
			r2.Value[i][j][0] = ekg.params.RingQP().NewPoly()
			r2.Value[i][j][1] = ekg.params.RingQP().NewPoly()
		}
	}
	return
}

// SampleCRP samples a common random polynomial to be used in the RKG protocol from the provided
// common reference string.
func (ekg *RKGProtocol) SampleCRP(crs CRS) RKGCRP {
	params := ekg.params
	decompRNS := params.DecompRNS(params.QCount()-1, params.PCount()-1)
	decompPw2 := params.DecompPw2(params.QCount()-1, params.PCount()-1)

	crp := make([][]ringqp.Poly, decompRNS)
	us := ringqp.NewUniformSampler(crs, *params.RingQP())
	for i := range crp {
		crp[i] = make([]ringqp.Poly, decompPw2)
		for j := range crp[i] {
			crp[i][j] = params.RingQP().NewPoly()
			us.Read(crp[i][j])
		}
	}
	return RKGCRP(crp)
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(sk *rlwe.SecretKey, crp RKGCRP, ephSkOut *rlwe.SecretKey, shareOut *RKGShare) {
	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	levelQ := sk.LevelQ()
	levelP := sk.LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ

	hasModulusP := levelP > -1

	if hasModulusP {
		// Computes P * sk
		ringQ.MulScalarBigint(sk.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], ekg.tmpPoly1.Q)
	} else {
		levelP = 0
		ring.CopyLvl(levelQ, sk.Value.Q, ekg.tmpPoly1.Q)
	}

	ringQ.InvMForm(ekg.tmpPoly1.Q, ekg.tmpPoly1.Q)

	// u
	ekg.ternarySamplerQ.Read(ephSkOut.Value.Q)
	if hasModulusP {
		ringQP.ExtendBasisSmallNormAndCenter(ephSkOut.Value.Q, levelP, nil, ephSkOut.Value.P)
	}
	ringQP.NTT(ephSkOut.Value, ephSkOut.Value)
	ringQP.MForm(ephSkOut.Value, ephSkOut.Value)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {
			// h = e
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j][0].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j][0].Q, levelP, nil, shareOut.Value[i][j][0].P)
			}

			ringQP.NTT(shareOut.Value[i][j][0], shareOut.Value[i][j][0])

			// h = sk*CrtBaseDecompQi + e
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handles the case where nb pj does not divides nb qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.Tables[index].Modulus
				skP := ekg.tmpPoly1.Q.Coeffs[index]
				h := shareOut.Value[i][j][0].Q.Coeffs[index]

				for w := 0; w < ringQ.N(); w++ {
					h[w] = ring.CRed(h[w]+skP[w], qi)
				}
			}

			// h = sk*CrtBaseDecompQi + -u*a + e
			ringQP.MulCoeffsMontgomeryAndSub(ephSkOut.Value, crp[i][j], shareOut.Value[i][j][0])

			// Second Element
			// e_2i
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j][1].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j][1].Q, levelP, nil, shareOut.Value[i][j][1].P)
			}

			ringQP.NTT(shareOut.Value[i][j][1], shareOut.Value[i][j][1])
			// s*a + e_2i
			ringQP.MulCoeffsMontgomeryAndAdd(sk.Value, crp[i][j], shareOut.Value[i][j][1])
		}

		ringQ.MulScalar(ekg.tmpPoly1.Q, 1<<ekg.params.Pow2Base(), ekg.tmpPoly1.Q)
	}
}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Upon receiving the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(ephSk, sk *rlwe.SecretKey, round1 *RKGShare, shareOut *RKGShare) {

	levelQ := sk.LevelQ()
	levelP := sk.LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	// (u_i - s_i)
	ringQP.Sub(ephSk.Value, sk.Value, ekg.tmpPoly1)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i:
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {

			// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

			// (AggregateShareRoundTwo samples) * sk
			ringQP.MulCoeffsMontgomeryConstant(round1.Value[i][j][0], sk.Value, shareOut.Value[i][j][0])

			// (AggregateShareRoundTwo samples) * sk + e_1i
			ekg.gaussianSamplerQ.Read(ekg.tmpPoly2.Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(ekg.tmpPoly2.Q, levelP, nil, ekg.tmpPoly2.P)
			}

			ringQP.NTT(ekg.tmpPoly2, ekg.tmpPoly2)
			ringQP.Add(shareOut.Value[i][j][0], ekg.tmpPoly2, shareOut.Value[i][j][0])

			// second part
			// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j][1].Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j][1].Q, levelP, nil, shareOut.Value[i][j][1].P)
			}

			ringQP.NTT(shareOut.Value[i][j][1], shareOut.Value[i][j][1])
			ringQP.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1, round1.Value[i][j][1], shareOut.Value[i][j][1])
		}
	}
}

// AggregateShares combines two RKG shares into a single one.
func (ekg *RKGProtocol) AggregateShares(share1, share2, shareOut *RKGShare) {

	levelQ := share1.Value[0][0][0].Q.Level()

	var levelP int
	if share1.Value[0][0][0].P != nil {
		levelP = share1.Value[0][0][0].P.Level()
	}

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(share1.Value[i][j][0], share2.Value[i][j][0], shareOut.Value[i][j][0])
			ringQP.Add(share1.Value[i][j][1], share2.Value[i][j][1], shareOut.Value[i][j][1])
		}
	}
}

// GenRelinearizationKey computes the generated RLK from the public shares and write the result in evalKeyOut.
func (ekg *RKGProtocol) GenRelinearizationKey(round1 *RKGShare, round2 *RKGShare, evalKeyOut *rlwe.RelinearizationKey) {

	levelQ := round1.Value[0][0][0].Q.Level()

	var levelP int
	if round1.Value[0][0][0].P != nil {
		levelP = round1.Value[0][0][0].P.Level()
	}

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(round1.Value)
	BITDecomp := len(round1.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(round2.Value[i][j][0], round2.Value[i][j][1], evalKeyOut.Keys[0].Value[i][j].Value[0])
			evalKeyOut.Keys[0].Value[i][j].Value[1].Copy(round1.Value[i][j][1])
			ringQP.MForm(evalKeyOut.Keys[0].Value[i][j].Value[0], evalKeyOut.Keys[0].Value[i][j].Value[0])
			ringQP.MForm(evalKeyOut.Keys[0].Value[i][j].Value[1], evalKeyOut.Keys[0].Value[i][j].Value[1])
		}
	}
}

// MarshalBinary encodes the target element on a slice of bytes.
func (share *RKGShare) MarshalBinary() ([]byte, error) {
	//we have modulus * bitLog * Len of 1 ring rings
	data := make([]byte, 2+2*share.Value[0][0][0].MarshalBinarySize64()*len(share.Value)*len(share.Value[0]))
	if len(share.Value) > 0xFF {
		return []byte{}, errors.New("RKGShare : uint8 overflow on length")
	}

	if len(share.Value[0]) > 0xFF {
		return []byte{}, errors.New("RKGShare : uint8 overflow on length")
	}

	data[0] = uint8(len(share.Value))
	data[1] = uint8(len(share.Value[0]))

	//write all of our rings in the data
	//write all the polys
	ptr := 2
	var inc int
	var err error
	for i := range share.Value {
		for _, el := range share.Value[i] {

			if inc, err = el[0].Encode64(data[ptr:]); err != nil {
				return []byte{}, err
			}
			ptr += inc

			if inc, err = el[1].Encode64(data[ptr:]); err != nil {
				return []byte{}, err
			}
			ptr += inc
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RKGShare) UnmarshalBinary(data []byte) (err error) {
	share.Value = make([][][2]ringqp.Poly, data[0])
	ptr := 2
	var inc int
	for i := range share.Value {
		share.Value[i] = make([][2]ringqp.Poly, data[1])
		for j := range share.Value[i] {

			if inc, err = share.Value[i][j][0].Decode64(data[ptr:]); err != nil {
				return err
			}
			ptr += inc

			if inc, err = share.Value[i][j][1].Decode64(data[ptr:]); err != nil {
				return err
			}
			ptr += inc
		}
	}

	return nil
}
