package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// RKGProtocolNaive is a structure storing the parameters for the naive EKG protocol.
type RKGProtocolNaive struct {
	dckksContext             *dckksContext
	polypool                 *ring.Poly
	gaussianSampler          *ring.GaussianSampler
	ternarySamplerMontgomery *ring.TernarySampler
}

// NewRKGProtocolNaive creates a new RKGProtocolNaive object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewRKGProtocolNaive(params *ckks.Parameters) (rkg *RKGProtocolNaive) {

	rkg = new(RKGProtocolNaive)
	dckksContext := newDckksContext(params)
	rkg.dckksContext = dckksContext
	rkg.polypool = dckksContext.ringQP.NewPoly()
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	rkg.gaussianSampler = ring.NewGaussianSampler(prng, dckksContext.ringQP, params.Sigma(), uint64(6*params.Sigma()))
	rkg.ternarySamplerMontgomery = ring.NewTernarySampler(prng, dckksContext.ringQP, 0.5, true)

	return
}

// RKGNaiveShareRoundOne is a struct storing the round one share of the RKG naive protocol.
type RKGNaiveShareRoundOne [][2]*ring.Poly

// RKGNaiveShareRoundTwo is a struct storing the round two share of the RKG naive protocol.
type RKGNaiveShareRoundTwo [][2]*ring.Poly

// AllocateShares allocates the share of the RKG naive protocol.
func (rkg *RKGProtocolNaive) AllocateShares() (r1 RKGNaiveShareRoundOne, r2 RKGNaiveShareRoundTwo) {
	contextQP := rkg.dckksContext.ringQP

	r1 = make([][2]*ring.Poly, rkg.dckksContext.beta)
	r2 = make([][2]*ring.Poly, rkg.dckksContext.beta)

	for i := uint64(0); i < rkg.dckksContext.beta; i++ {
		r1[i][0] = contextQP.NewPoly()
		r1[i][1] = contextQP.NewPoly()

		r2[i][0] = contextQP.NewPoly()
		r2[i][1] = contextQP.NewPoly()
	}

	return
}

// GenShareRoundOne is the first of two rounds of the naive EKG protocol. Using the shared public key "cpk",
// each party generates a pseudo-encryption of s*w of the form :
//
// [cpk[0] * u_i + P * s_i + e_0i, cpk[1] * u_i + e_1i]
//
// and broadcasts it to all other j-1 parties.
func (rkg *RKGProtocolNaive) GenShareRoundOne(sk *ring.Poly, pk [2]*ring.Poly, shareOut RKGNaiveShareRoundOne) {

	contextQP := rkg.dckksContext.ringQP

	rkg.polypool.Copy(sk)

	contextQP.MulScalarBigint(rkg.polypool, rkg.dckksContext.ringP.ModulusBigint, rkg.polypool)

	contextQP.InvMForm(rkg.polypool, rkg.polypool)

	var index uint64

	for i := uint64(0); i < rkg.dckksContext.beta; i++ {

		// h_0 = e0
		rkg.gaussianSampler.Read(shareOut[i][0])
		contextQP.NTT(shareOut[i][0], shareOut[i][0])

		// h_1 = e1
		rkg.gaussianSampler.Read(shareOut[i][1])
		contextQP.NTT(shareOut[i][1], shareOut[i][1])

		// h_0 = e0 + [sk*P*(qiBarre*qiStar)%qi = sk*P, else 0]

		for j := uint64(0); j < rkg.dckksContext.alpha; j++ {

			index = i*rkg.dckksContext.alpha + j

			qi := contextQP.Modulus[index]
			tmp0 := rkg.polypool.Coeffs[index]
			tmp1 := shareOut[i][0].Coeffs[index]

			for w := uint64(0); w < contextQP.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= rkg.dckksContext.params.QiCount() {
				break
			}
		}
	}

	for i := uint64(0); i < rkg.dckksContext.beta; i++ {
		// u
		rkg.ternarySamplerMontgomery.Read(rkg.polypool)
		contextQP.NTT(rkg.polypool, rkg.polypool)
		// h_0 = pk_0 * u + e0 + P * sk * (qiBarre*qiStar)%qi
		contextQP.MulCoeffsMontgomeryAndAdd(pk[0], rkg.polypool, shareOut[i][0])
		// h_1 = pk_1 * u + e1 + P * sk * (qiBarre*qiStar)%qi
		contextQP.MulCoeffsMontgomeryAndAdd(pk[1], rkg.polypool, shareOut[i][1])
	}

	rkg.polypool.Zero()
}

// AggregateShareRoundOne is the second part of the first round of the naive EKG protocol. Upon receiving the j-1 elements, each party computes :
//
// [sum(cpk[0] * u_j + P * s_j + e_0j), sum(cpk[1] * u_j + e_1j)]
//
// = [cpk[0] * u + P * s + e_0, cpk[1]*u + e_1]
func (rkg *RKGProtocolNaive) AggregateShareRoundOne(share1, share2, shareOut RKGNaiveShareRoundOne) {

	contextQP := rkg.dckksContext.ringQP

	for i := uint64(0); i < rkg.dckksContext.beta; i++ {
		contextQP.Add(share1[i][0], share2[i][0], shareOut[i][0])
		contextQP.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}

}

// GenShareRoundTwo is the first part of the second round, each party computes :
//
// [s_i * (cpk[0] * u + P * s + e_0) + v_i * cpk[0] + e_2i, s_i * (cpk[1] * u + e_1) + cpk[1] * v_i + e_3i]
//
// = [ cpk[0] * (u * s_i) + P * (s * s_i) + (s_i * e_0) + v_i*cpk[0] + e_2i, cpk[1] * u * s_i + (s_i * e_1) + cpk[1] * v_i + e_3i]
//
// And each party broadcasts this last result to the other j-1 parties.
func (rkg *RKGProtocolNaive) GenShareRoundTwo(round1 RKGNaiveShareRoundOne, sk *ring.Poly, pk [2]*ring.Poly, shareOut RKGNaiveShareRoundTwo) {

	contextQP := rkg.dckksContext.ringQP

	for i := uint64(0); i < rkg.dckksContext.beta; i++ {

		// h_0 = sum(samples[0]) * sk
		// h_1 = sum(samples[1]) * sk
		contextQP.MulCoeffsMontgomery(round1[i][0], sk, shareOut[i][0])
		contextQP.MulCoeffsMontgomery(round1[i][1], sk, shareOut[i][1])

		// v
		rkg.ternarySamplerMontgomery.Read(rkg.polypool)
		contextQP.NTT(rkg.polypool, rkg.polypool)

		// h_0 = sum(samples[0]) * sk + pk0 * v
		contextQP.MulCoeffsMontgomeryAndAdd(pk[0], rkg.polypool, shareOut[i][0])

		// h_1 = sum(samples[1]) * sk + pk1 * v
		contextQP.MulCoeffsMontgomeryAndAdd(pk[1], rkg.polypool, shareOut[i][1])

		// h_0 = sum(samples[0]) * sk + pk0 * v + e2
		rkg.gaussianSampler.Read(rkg.polypool)
		contextQP.NTT(rkg.polypool, rkg.polypool)
		contextQP.Add(shareOut[i][0], rkg.polypool, shareOut[i][0])

		// h_1 = sum(samples[1]) * sk + pk1 * v + e3
		rkg.gaussianSampler.Read(rkg.polypool)
		contextQP.NTT(rkg.polypool, rkg.polypool)
		contextQP.Add(shareOut[i][1], rkg.polypool, shareOut[i][1])

	}

	rkg.polypool.Zero()
}

// AggregateShareRoundTwo is the second part of the second and last round of the naive EKG protocol. Upon receiving the j-1 elements,
// each party computes :
//
// [ sum(cpk[0] * (u*s_i) + P * (s*s_i) + (s_i*e_0) + v_i*cpk[0] + e_2i), sum(cpk[1]*u*s_i + (s_i*e_1) + cpk[1] * v_i + e_3i)]
//
// = [cpk[0] * (s*u + v) + (P * s^2 ) + s*e_0 + e_2, ckp[1] * (s*u + v) + s*e_1 + e_3]
//
// = [-s*b + P * s^2 - (s*u + b) * e_cpk + s*e_0 + e_2, b + s*e_1 + e_3]
func (rkg *RKGProtocolNaive) AggregateShareRoundTwo(share1, share2, shareOut RKGNaiveShareRoundTwo) {

	contextQP := rkg.dckksContext.ringQP

	for i := uint64(0); i < rkg.dckksContext.beta; i++ {
		contextQP.Add(share1[i][0], share2[i][0], shareOut[i][0])
		contextQP.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}
}

// GenRelinearizationKey finalizes the protocol and returns the collective EvalutionKey.
func (rkg *RKGProtocolNaive) GenRelinearizationKey(round2 RKGNaiveShareRoundTwo, evalKeyOut *ckks.EvaluationKey) {

	contextQP := rkg.dckksContext.ringQP

	key := evalKeyOut.Get().Get()
	for i := uint64(0); i < rkg.dckksContext.beta; i++ {

		key[i][0].Copy(round2[i][0])
		key[i][1].Copy(round2[i][1])

		contextQP.MForm(key[i][0], key[i][0])
		contextQP.MForm(key[i][1], key[i][1])

	}
}
