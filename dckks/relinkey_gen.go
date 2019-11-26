package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

// RKGProtocol is a structure storing the parameters for the collective evaluation-key generation.
type RKGProtocol struct {
	ckksContext *ckks.Context
	polypool    *ring.Poly
}

// RKGShareRoundOne is a struct storing the round one share of the RKG protocol.
type RKGShareRoundOne []*ring.Poly

// RKGShareRoundTwo is a struct storing the round two share of the RKG protocol.
type RKGShareRoundTwo [][2]*ring.Poly

// RKGShareRoundThree is a struct storing the round three share of the RKG protocol.
type RKGShareRoundThree []*ring.Poly

// AllocateShares allocates the shares of the RKG protocol.
func (ekg *RKGProtocol) AllocateShares() (r1 RKGShareRoundOne, r2 RKGShareRoundTwo, r3 RKGShareRoundThree) {

	contextKeys := ekg.ckksContext.ContextKeys()

	r1 = make([]*ring.Poly, ekg.ckksContext.Beta())
	r2 = make([][2]*ring.Poly, ekg.ckksContext.Beta())
	r3 = make([]*ring.Poly, ekg.ckksContext.Beta())
	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {
		r1[i] = contextKeys.NewPoly()
		r2[i][0] = contextKeys.NewPoly()
		r2[i][1] = contextKeys.NewPoly()
		r3[i] = contextKeys.NewPoly()
	}
	return
}

// NewEkgProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key.
func NewEkgProtocol(ckksContext *ckks.Context) *RKGProtocol {

	ekg := new(RKGProtocol)
	ekg.ckksContext = ckksContext
	ekg.polypool = ekg.ckksContext.ContextKeys().NewPoly()

	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first rounds).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *RKGProtocol) NewEphemeralKey(p float64) (ephemeralKey *ring.Poly) {
	return ekg.ckksContext.ContextKeys().SampleTernaryMontgomeryNTTNew(p)
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + P*s_i + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(u, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShareRoundOne) {

	contextKeys := ekg.ckksContext.ContextKeys()

	var index uint64

	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	ekg.polypool.Copy(sk)

	contextKeys.MulScalarBigint(ekg.polypool, ekg.ckksContext.ContextP().ModulusBigint, ekg.polypool)

	contextKeys.InvMForm(ekg.polypool, ekg.polypool)

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {

		// h = e
		ekg.ckksContext.GaussianSampler().SampleNTT(shareOut[i])

		// h = sk*CrtBaseDecompQi + e
		for j := uint64(0); j < ekg.ckksContext.Alpha(); j++ {

			index = i*ekg.ckksContext.Alpha() + j

			qi := contextKeys.Modulus[index]
			tmp0 := ekg.polypool.Coeffs[index]
			tmp1 := shareOut[i].Coeffs[index]

			for w := uint64(0); w < contextKeys.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= uint64(len(ekg.ckksContext.ContextQ().Modulus)-1) {
				break
			}
		}

		// h = sk*CrtBaseDecompQi + -u*a + e
		contextKeys.MulCoeffsMontgomeryAndSub(u, crp[i], shareOut[i])
	}

	ekg.polypool.Zero()

	return
}

// AggregateShareRoundOne sums share1 with share2 on shareOut.
func (ekg *RKGProtocol) AggregateShareRoundOne(share1, share2, shareOut RKGShareRoundOne) {

	contextKeys := ekg.ckksContext.ContextKeys()

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {
		contextKeys.Add(share1[i], share2[i], shareOut[i])
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

	contextKeys := ekg.ckksContext.ContextKeys()

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// (AggregateShareRoundTwo samples) * sk
		contextKeys.MulCoeffsMontgomery(round1[i], sk, shareOut[i][0])

		// (AggregateShareRoundTwo samples) * sk + e_1i
		ekg.ckksContext.GaussianSampler().SampleNTT(ekg.polypool)
		contextKeys.Add(shareOut[i][0], ekg.polypool, shareOut[i][0])

		// Second Element
		// e_2i
		ekg.ckksContext.GaussianSampler().SampleNTT(shareOut[i][1])
		// s*a + e_2i
		contextKeys.MulCoeffsMontgomeryAndAdd(sk, crp[i], shareOut[i][1])
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

	contextKeys := ekg.ckksContext.ContextKeys()

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {
		contextKeys.Add(share1[i][0], share2[i][0], shareOut[i][0])
		contextKeys.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}

}

// GenShareRoundThree is the second pard of the third and last round of the RKGProtocol protocol. Each party operates a key-switch on [s*a + e_2],
// by computing :
//
// [(u_i - s_i)*(s*a + e_2)]
//
// and broadcasts the result to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundThree(round2 RKGShareRoundTwo, u, sk *ring.Poly, shareOut RKGShareRoundThree) {

	contextKeys := ekg.ckksContext.ContextKeys()

	// (u_i - s_i)
	contextKeys.Sub(u, sk, ekg.polypool)

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {

		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		ekg.ckksContext.GaussianSampler().SampleNTT(shareOut[i])
		contextKeys.MulCoeffsMontgomeryAndAdd(ekg.polypool, round2[i][1], shareOut[i])
	}
}

// AggregateShareRoundThree sums share1 with share2 on shareOut.
func (ekg *RKGProtocol) AggregateShareRoundThree(share1, share2, shareOut RKGShareRoundThree) {

	contextKeys := ekg.ckksContext.ContextKeys()

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {
		contextKeys.Add(share1[i], share2[i], shareOut[i])
	}
}

// GenRelinearizationKey finalizes the protocol and returns the collective EvalutionKey.
func (ekg *RKGProtocol) GenRelinearizationKey(round2 RKGShareRoundTwo, round3 RKGShareRoundThree, evalKeyOut *ckks.EvaluationKey) {

	contextKeys := ekg.ckksContext.ContextKeys()

	key := evalKeyOut.Get().Get()
	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {

		contextKeys.Add(round2[i][0], round3[i], key[i][0])
		key[i][1].Copy(round2[i][1])

		contextKeys.MForm(key[i][0], key[i][0])
		contextKeys.MForm(key[i][1], key[i][1])

	}
}
