package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// RKGProtocol is a structure storing the parameters for the collective evaluation-key generation.
type RKGProtocol struct {
	dckksContext    *dckksContext
	polypool        *ring.Poly
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
}

// RKGShareRoundOne is a struct storing the round one share of the RKG protocol.
type RKGShareRoundOne []*ring.Poly

// RKGShareRoundTwo is a struct storing the round two share of the RKG protocol.
type RKGShareRoundTwo [][2]*ring.Poly

// RKGShareRoundThree is a struct storing the round three share of the RKG protocol.
type RKGShareRoundThree []*ring.Poly

// AllocateShares allocates the shares of the RKG protocol.
func (ekg *RKGProtocol) AllocateShares() (r1 RKGShareRoundOne, r2 RKGShareRoundTwo, r3 RKGShareRoundThree) {

	contextQP := ekg.dckksContext.contextQP

	r1 = make([]*ring.Poly, ekg.dckksContext.beta)
	r2 = make([][2]*ring.Poly, ekg.dckksContext.beta)
	r3 = make([]*ring.Poly, ekg.dckksContext.beta)
	for i := uint64(0); i < ekg.dckksContext.beta; i++ {
		r1[i] = contextQP.NewPoly()
		r2[i][0] = contextQP.NewPoly()
		r2[i][1] = contextQP.NewPoly()
		r3[i] = contextQP.NewPoly()
	}
	return
}

// NewEkgProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key.
func NewEkgProtocol(params *ckks.Parameters) *RKGProtocol {

	if !params.IsValid() {
		panic("cannot NewEkgProtocol : params not valid (check if they where generated properly)")
	}

	ekg := new(RKGProtocol)
	dckksContext := newDckksContext(params)
	ekg.dckksContext = dckksContext
	ekg.polypool = ekg.dckksContext.contextQP.NewPoly()
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ekg.gaussianSampler = ring.NewGaussianSampler(prng, dckksContext.contextQP)
	ekg.ternarySampler = ring.NewTernarySampler(prng, dckksContext.contextQP, 0.5, true)
	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first rounds).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *RKGProtocol) NewEphemeralKey() (ephemeralKey *ring.Poly) {
	return ekg.ternarySampler.ReadNewNTT()
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + P*s_i + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(u, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShareRoundOne) {

	contextQP := ekg.dckksContext.contextQP

	var index uint64

	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	ekg.polypool.Copy(sk)

	contextQP.MulScalarBigint(ekg.polypool, ekg.dckksContext.contextP.ModulusBigint, ekg.polypool)

	contextQP.InvMForm(ekg.polypool, ekg.polypool)

	for i := uint64(0); i < ekg.dckksContext.beta; i++ {

		// h = e
		ekg.gaussianSampler.SampleNTTLvl(uint64(len(contextQP.Modulus)-1), shareOut[i], ekg.dckksContext.params.Sigma, uint64(6*ekg.dckksContext.params.Sigma))

		// h = sk*CrtBaseDecompQi + e
		for j := uint64(0); j < ekg.dckksContext.alpha; j++ {

			index = i*ekg.dckksContext.alpha + j

			qi := contextQP.Modulus[index]
			tmp0 := ekg.polypool.Coeffs[index]
			tmp1 := shareOut[i].Coeffs[index]

			for w := uint64(0); w < contextQP.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= uint64(len(ekg.dckksContext.params.Qi)-1) {
				break
			}
		}

		// h = sk*CrtBaseDecompQi + -u*a + e
		contextQP.MulCoeffsMontgomeryAndSub(u, crp[i], shareOut[i])
	}

	ekg.polypool.Zero()

	return
}

// AggregateShareRoundOne sums share1 with share2 on shareOut.
func (ekg *RKGProtocol) AggregateShareRoundOne(share1, share2, shareOut RKGShareRoundOne) {

	contextQP := ekg.dckksContext.contextQP

	for i := uint64(0); i < ekg.dckksContext.beta; i++ {
		contextQP.Add(share1[i], share2[i], shareOut[i])
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

	contextQP := ekg.dckksContext.contextQP

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := uint64(0); i < ekg.dckksContext.beta; i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// (AggregateShareRoundTwo samples) * sk
		contextQP.MulCoeffsMontgomery(round1[i], sk, shareOut[i][0])

		// (AggregateShareRoundTwo samples) * sk + e_1i
		ekg.gaussianSampler.SampleNTTLvl(uint64(len(contextQP.Modulus)-1), ekg.polypool, ekg.dckksContext.params.Sigma, uint64(6*ekg.dckksContext.params.Sigma))
		contextQP.Add(shareOut[i][0], ekg.polypool, shareOut[i][0])

		// Second Element
		// e_2i
		ekg.gaussianSampler.SampleNTTLvl(uint64(len(contextQP.Modulus)-1), shareOut[i][1], ekg.dckksContext.params.Sigma, uint64(6*ekg.dckksContext.params.Sigma))
		// s*a + e_2i
		contextQP.MulCoeffsMontgomeryAndAdd(sk, crp[i], shareOut[i][1])
	}

	ekg.polypool.Zero()

}

// AggregateShareRoundTwo is the first part of the third and last round of the RKGProtocol protocol. Upon receiving the j-1 elements, each party
// computes :
//
// [sum(s_j * (-u*a + s*w + e) + e_j1), sum(s_j*a + e_j2)]
//
// = [s * (-u*a + s*w + e) + e_1, s*a + e_2].
func (ekg *RKGProtocol) AggregateShareRoundTwo(share1, share2, shareOut RKGShareRoundTwo) {

	contextQP := ekg.dckksContext.contextQP

	for i := uint64(0); i < ekg.dckksContext.beta; i++ {
		contextQP.Add(share1[i][0], share2[i][0], shareOut[i][0])
		contextQP.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}

}

// GenShareRoundThree is the second pard of the third and last round of the RKGProtocol protocol. Each party operates a key-switch on [s*a + e_2],
// by computing :
//
// [(u_i - s_i)*(s*a + e_2)]
//
// and broadcasts the result to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundThree(round2 RKGShareRoundTwo, u, sk *ring.Poly, shareOut RKGShareRoundThree) {

	contextQP := ekg.dckksContext.contextQP

	// (u_i - s_i)
	contextQP.Sub(u, sk, ekg.polypool)

	for i := uint64(0); i < ekg.dckksContext.beta; i++ {

		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		ekg.gaussianSampler.SampleNTTLvl(uint64(len(contextQP.Modulus)-1), shareOut[i], ekg.dckksContext.params.Sigma, uint64(6*ekg.dckksContext.params.Sigma))
		contextQP.MulCoeffsMontgomeryAndAdd(ekg.polypool, round2[i][1], shareOut[i])
	}
}

// AggregateShareRoundThree sums share1 with share2 on shareOut.
func (ekg *RKGProtocol) AggregateShareRoundThree(share1, share2, shareOut RKGShareRoundThree) {

	contextQP := ekg.dckksContext.contextQP

	for i := uint64(0); i < ekg.dckksContext.beta; i++ {
		contextQP.Add(share1[i], share2[i], shareOut[i])
	}
}

// GenRelinearizationKey finalizes the protocol and returns the collective EvalutionKey.
func (ekg *RKGProtocol) GenRelinearizationKey(round2 RKGShareRoundTwo, round3 RKGShareRoundThree, evalKeyOut *ckks.EvaluationKey) {

	contextQP := ekg.dckksContext.contextQP

	key := evalKeyOut.Get().Get()
	for i := uint64(0); i < ekg.dckksContext.beta; i++ {

		contextQP.Add(round2[i][0], round3[i], key[i][0])
		key[i][1].Copy(round2[i][1])

		contextQP.MForm(key[i][0], key[i][0])
		contextQP.MForm(key[i][1], key[i][1])

	}
}
