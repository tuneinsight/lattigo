package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"math"
)

// rkgProtocolState is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type rkgProtocolState struct {
	ringContext     *ring.Context
	ternarySampler  *ring.TernarySampler
	gaussianSampler *ring.KYSampler
	bitDecomp       uint64
	bitLog          uint64
	tmpPoly1        *ring.Poly
	tmpPoly2        *ring.Poly
}

type rkgShareRoundOne [][]*ring.Poly
type rkgShareRoundTwo [][][2]*ring.Poly
type rkgShareRoundThree [][]*ring.Poly

func (ekg *rkgProtocolState) AllocateShares() (r1 rkgShareRoundOne, r2 rkgShareRoundTwo, r3 rkgShareRoundThree) {
	r1 = make([][]*ring.Poly, len(ekg.ringContext.Modulus))
	r2 = make([][][2]*ring.Poly, len(ekg.ringContext.Modulus))
	r3 = make([][]*ring.Poly, len(ekg.ringContext.Modulus))
	for i := range ekg.ringContext.Modulus {
		r1[i] = make([]*ring.Poly, ekg.bitLog)
		r2[i] = make([][2]*ring.Poly, ekg.bitLog)
		r3[i] = make([]*ring.Poly, ekg.bitLog)
		for w := uint64(0); w < ekg.bitLog; w++ {
			r1[i][w] = ekg.ringContext.NewPoly()
			r2[i][w][0] = ekg.ringContext.NewPoly()
			r2[i][w][1] = ekg.ringContext.NewPoly()
			r3[i][w] = ekg.ringContext.NewPoly()
		}
	}
	return
}

// NewEkgProtocol creates a new rkgProtocolState object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewEkgProtocol(context *bfv.BfvContext, bitDecomp uint64) *rkgProtocolState {
	ekg := new(rkgProtocolState)
	ekg.ringContext = context.ContextQ()
	ekg.ternarySampler = context.TernarySampler()
	ekg.gaussianSampler = context.GaussianSampler()
	ekg.bitDecomp = bitDecomp
	ekg.bitLog = uint64(math.Ceil(float64(60) / float64(bitDecomp)))
	ekg.tmpPoly1 = ekg.ringContext.NewPoly()
	ekg.tmpPoly2 = ekg.ringContext.NewPoly()
	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first round).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *rkgProtocolState) NewEphemeralKey(p float64) (ephemeralKey *ring.Poly, err error) {
	if ephemeralKey, err = ekg.ternarySampler.SampleMontgomeryNTTNew(p); err != nil {
		return nil, err
	}
	return
}

// GenShareRoundOne is the first of three rounds of the rkgProtocolState protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *rkgProtocolState) GenShareRoundOne(u, sk *ring.Poly, crp [][]*ring.Poly, shareOut rkgShareRoundOne) {

	mredParams := ekg.ringContext.GetMredParams()

	// Given a base decomposition w (here the CRT decomposition)
	// computes [-u_i*a + s_i*w + e_i]
	// where a = crp
	for i, qi := range ekg.ringContext.Modulus {

		for w := uint64(0); w < ekg.bitLog; w++ {

			// h = e
			ekg.gaussianSampler.SampleNTT(shareOut[i][w])

			// h = sk*CrtBaseDecompQi + e
			for j := uint64(0); j < ekg.ringContext.N; j++ {
				shareOut[i][w].Coeffs[i][j] += ring.PowerOf2(sk.Coeffs[i][j], ekg.bitDecomp*w, qi, mredParams[i])
			}

			// h = sk*CrtBaseDecompQi + -u*a + e
			ekg.ringContext.MulCoeffsMontgomeryAndSub(u, crp[i][w], shareOut[i][w])
		}
	}

	return
}

func (ekg *rkgProtocolState) AggregateShareRoundOne(share1, share2, shareOut rkgShareRoundOne) {

	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			ekg.ringContext.Add(share1[i][w], share2[i][w], shareOut[i][w])
		}
	}
}

// GenShareRoundTwo is the second of three rounds of the rkgProtocolState protocol. Uppon received the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *rkgProtocolState) GenShareRoundTwo(round1 rkgShareRoundOne, sk *ring.Poly, crp [][]*ring.Poly, shareOut rkgShareRoundTwo) {

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {

			// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

			// (AggregateShareRoundTwo samples) * sk
			ekg.ringContext.MulCoeffsMontgomery(round1[i][w], sk, shareOut[i][w][0])

			// (AggregateShareRoundTwo samples) * sk + e_1i
			ekg.gaussianSampler.SampleNTT(ekg.tmpPoly1)
			ekg.ringContext.Add(shareOut[i][w][0], ekg.tmpPoly1, shareOut[i][w][0])

			// Second Element
			// e_2i
			ekg.gaussianSampler.SampleNTT(shareOut[i][w][1])
			// s*a + e_2i
			ekg.ringContext.MulCoeffsMontgomeryAndAdd(sk, crp[i][w], shareOut[i][w][1])
		}
	}
	ekg.tmpPoly1.Zero()
}

// AggregateShareRoundTwo is the first part of the third and last round of the rkgProtocolState protocol. Uppon receiving the j-1 elements, each party
// computues :
//
// [sum(s_j * (-u*a + s*w + e) + e_j1), sum(s_j*a + e_j2)]
//
// = [s * (-u*a + s*w + e) + e_1, s*a + e_2].
func (ekg *rkgProtocolState) AggregateShareRoundTwo(share1, share2, shareOut rkgShareRoundTwo) {

	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			ekg.ringContext.Add(share1[i][w][0], share2[i][w][0], shareOut[i][w][0])
			ekg.ringContext.Add(share1[i][w][1], share2[i][w][1], shareOut[i][w][1])
		}
	}
}

// GenShareRound3 is the second pard of the third and last round of the rkgProtocolState protocol. Each party operates a key-switch on [s*a + e_2],
// by computing :
//
// [(u_i - s_i)*(s*a + e_2)]
//
// and broadcasts the result the other j-1 parties.
func (ekg *rkgProtocolState) GenShareRound3(round2 rkgShareRoundTwo, u, sk *ring.Poly, shareOut rkgShareRoundThree) {

	// (u_i - s_i)
	ekg.ringContext.Sub(u, sk, ekg.tmpPoly1)

	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
			ekg.gaussianSampler.SampleNTT(shareOut[i][w])
			ekg.ringContext.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1, round2[i][w][1], shareOut[i][w])
		}
	}
}


func (ekg *rkgProtocolState) AggregateShareRound3(share1, share2, share3 rkgShareRoundThree) {
	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			ekg.ringContext.Add(share1[i][w], share2[i][w], share3[i][w])
		}
	}
}


func (ekg *rkgProtocolState) GenRelinearizationKey(round2 rkgShareRoundTwo, round3 rkgShareRoundThree, evalKeyOut *bfv.EvaluationKey) {

	evk := [][][][2]*ring.Poly{make([][][2]*ring.Poly, len(ekg.ringContext.Modulus))}
	key := evk[0]
	for i := range ekg.ringContext.Modulus {
		key[i] = make([][2]*ring.Poly, ekg.bitLog)
		for w := uint64(0); w < ekg.bitLog; w++ {
			key[i][w][0] = ekg.ringContext.NewPoly()
			ekg.ringContext.Add(round2[i][w][0], round3[i][w], key[i][w][0])
			key[i][w][1] = round2[i][w][1] // A propper copy will be done by evk.SetRelinKeys


			ekg.ringContext.MForm(key[i][w][0], key[i][w][0])
			ekg.ringContext.MForm(key[i][w][1], key[i][w][1])
		}
	}

	evalKeyOut.SetRelinKeys(evk, ekg.bitLog)
}