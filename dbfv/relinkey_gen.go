package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"math"
)

// RKGProtocol is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type RKGProtocol struct {
	ringContext     *ring.Context
	keyswitchprimes []uint64
	alpha           uint64
	beta            uint64
	gaussianSampler *ring.KYSampler
	tmpPoly1        *ring.Poly
	tmpPoly2        *ring.Poly
	polypool        *ring.Poly
}

type RKGShareRoundOne []*ring.Poly
type RKGShareRoundTwo [][2]*ring.Poly
type RKGShareRoundThree []*ring.Poly

func (ekg *RKGProtocol) AllocateShares() (r1 RKGShareRoundOne, r2 RKGShareRoundTwo, r3 RKGShareRoundThree) {
	r1 = make([]*ring.Poly, ekg.beta)
	r2 = make([][2]*ring.Poly, ekg.beta)
	r3 = make([]*ring.Poly, ekg.beta)
	for i := uint64(0); i < ekg.beta; i++ {
		r1[i] = ekg.ringContext.NewPoly()
		r2[i][0] = ekg.ringContext.NewPoly()
		r2[i][1] = ekg.ringContext.NewPoly()
		r3[i] = ekg.ringContext.NewPoly()
	}
	return
}

// NewEkgProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewEkgProtocol(context *bfv.BfvContext) *RKGProtocol {

	ekg := new(RKGProtocol)
	ekg.ringContext = context.ContextKeys()

	ekg.keyswitchprimes = make([]uint64, len(context.KeySwitchPrimes()))
	for i, pi := range context.KeySwitchPrimes() {
		ekg.keyswitchprimes[i] = pi
	}

	ekg.alpha = uint64(len(ekg.keyswitchprimes))
	ekg.beta = uint64(math.Ceil(float64(len(ekg.ringContext.Modulus)-len(ekg.keyswitchprimes)) / float64(ekg.alpha)))

	ekg.gaussianSampler = context.GaussianSampler()

	ekg.tmpPoly1 = ekg.ringContext.NewPoly()
	ekg.tmpPoly2 = ekg.ringContext.NewPoly()
	ekg.polypool = ekg.ringContext.NewPoly()

	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first round).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *RKGProtocol) NewEphemeralKey(p float64) (ephemeralKey *ring.Poly, err error) {
	if ephemeralKey, err = ekg.ringContext.SampleTernaryMontgomeryNTTNew(p); err != nil {
		return nil, err
	}
	return
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

	for _, pj := range ekg.keyswitchprimes {
		ekg.ringContext.MulScalar(ekg.polypool, pj, ekg.polypool)
	}

	ekg.ringContext.InvMForm(ekg.polypool, ekg.polypool)

	for i := uint64(0); i < ekg.beta; i++ {

		// h = e
		ekg.gaussianSampler.SampleNTT(shareOut[i])

		// h = sk*CrtBaseDecompQi + e
		for j := uint64(0); j < ekg.alpha; j++ {

			index = i*ekg.alpha + j

			for w := uint64(0); w < ekg.ringContext.N; w++ {
				shareOut[i].Coeffs[index][w] = ring.CRed(shareOut[i].Coeffs[index][w]+ekg.polypool.Coeffs[index][w], ekg.ringContext.Modulus[index])
			}

			// Handles the case where nb pj does not divides nb qi
			if index == uint64(len(ekg.ringContext.Modulus)-len(ekg.keyswitchprimes)-1) {
				break
			}
		}

		// h = sk*CrtBaseDecompQi + -u*a + e
		ekg.ringContext.MulCoeffsMontgomeryAndSub(u, crp[i], shareOut[i])
	}

	ekg.polypool.Zero()

	return
}

func (ekg *RKGProtocol) AggregateShareRoundOne(share1, share2, shareOut RKGShareRoundOne) {

	for i := uint64(0); i < ekg.beta; i++ {
		ekg.ringContext.Add(share1[i], share2[i], shareOut[i])
	}

}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Uppon received the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(round1 RKGShareRoundOne, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShareRoundTwo) {

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := uint64(0); i < ekg.beta; i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// (AggregateShareRoundTwo samples) * sk
		ekg.ringContext.MulCoeffsMontgomery(round1[i], sk, shareOut[i][0])

		// (AggregateShareRoundTwo samples) * sk + e_1i
		ekg.gaussianSampler.SampleNTT(ekg.tmpPoly1)
		ekg.ringContext.Add(shareOut[i][0], ekg.tmpPoly1, shareOut[i][0])

		// Second Element
		// e_2i
		ekg.gaussianSampler.SampleNTT(shareOut[i][1])
		// s*a + e_2i
		ekg.ringContext.MulCoeffsMontgomeryAndAdd(sk, crp[i], shareOut[i][1])
	}

}

// AggregateShareRoundTwo is the first part of the third and last round of the RKGProtocol protocol. Uppon receiving the j-1 elements, each party
// computues :
//
// [sum(s_j * (-u*a + s*w + e) + e_j1), sum(s_j*a + e_j2)]
//
// = [s * (-u*a + s*w + e) + e_1, s*a + e_2].
func (ekg *RKGProtocol) AggregateShareRoundTwo(share1, share2, shareOut RKGShareRoundTwo) {

	for i := uint64(0); i < ekg.beta; i++ {
		ekg.ringContext.Add(share1[i][0], share2[i][0], shareOut[i][0])
		ekg.ringContext.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}

}

// GenShareRoundThree is the second pard of the third and last round of the RKGProtocol protocol. Each party operates a key-switch on [s*a + e_2],
// by computing :
//
// [(u_i - s_i)*(s*a + e_2)]
//
// and broadcasts the result the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundThree(round2 RKGShareRoundTwo, u, sk *ring.Poly, shareOut RKGShareRoundThree) {

	// (u_i - s_i)
	ekg.ringContext.Sub(u, sk, ekg.tmpPoly1)

	for i := uint64(0); i < ekg.beta; i++ {

		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		ekg.gaussianSampler.SampleNTT(shareOut[i])
		ekg.ringContext.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1, round2[i][1], shareOut[i])
	}
}

func (ekg *RKGProtocol) AggregateShareRoundThree(share1, share2, shareOut RKGShareRoundThree) {
	for i := uint64(0); i < ekg.beta; i++ {
		ekg.ringContext.Add(share1[i], share2[i], shareOut[i])
	}
}

func (ekg *RKGProtocol) GenRelinearizationKey(round2 RKGShareRoundTwo, round3 RKGShareRoundThree, evalKeyOut *bfv.EvaluationKey) {

	key := evalKeyOut.Get()[0].Get()
	for i := uint64(0); i < ekg.beta; i++ {

		ekg.ringContext.Add(round2[i][0], round3[i], key[i][0])
		key[i][1].Copy(round2[i][1])

		ekg.ringContext.MForm(key[i][0], key[i][0])
		ekg.ringContext.MForm(key[i][1], key[i][1])

	}
}
