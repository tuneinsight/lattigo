package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

// EkgProtocolNaive is a structure storing the parameters for the naive EKG protocol.
type EkgProtocolNaive struct {
	ckksContext     *ckks.CkksContext
	gaussianSampler *ring.KYSampler
	ternarySampler  *ring.TernarySampler
	polypool        *ring.Poly
}

// NewEkgProtocolNaive creates a new EkgProtocolNaive object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewEkgProtocolNaive(ckksContext *ckks.CkksContext) (ekg *EkgProtocolNaive) {
	ekg = new(EkgProtocolNaive)
	ekg.ckksContext = ckksContext
	ekg.gaussianSampler = ckksContext.GaussianSampler()
	ekg.ternarySampler = ckksContext.TernarySampler()
	ekg.polypool = ckksContext.ContextKeys().NewPoly()
	return ekg
}

// GenSamples is the first of two rounds of the naive EKG protocol. Using the shared public key "cpk",
// each party generates a pseudo-encryption of s*w of the form :
//
// [cpk[0] * u_i + P * s_i + e_0i, cpk[1] * u_i + e_1i]
//
// and broadcasts it to all other j-1 parties.
func (ekg *EkgProtocolNaive) GenSamples(sk *ring.Poly, pk [2]*ring.Poly) (h [][2]*ring.Poly) {

	contextKeys := ekg.ckksContext.ContextKeys()

	h = make([][2]*ring.Poly, ekg.ckksContext.Beta())

	ekg.polypool.Copy(sk)

	for _, pj := range ekg.ckksContext.KeySwitchPrimes() {
		contextKeys.MulScalar(ekg.polypool, pj, ekg.polypool)
	}

	contextKeys.InvMForm(ekg.polypool, ekg.polypool)

	var index uint64

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {
		// h_0 = e0
		h[i][0] = ekg.gaussianSampler.SampleNTTNew()
		// h_1 = e1
		h[i][1] = ekg.gaussianSampler.SampleNTTNew()

		// h_0 = e0 + [sk*P*(qiBarre*qiStar)%qi = sk*P, else 0]

		for j := uint64(0); j < ekg.ckksContext.Alpha(); j++ {

			index = i*ekg.ckksContext.Alpha() + j

			for w := uint64(0); w < contextKeys.N; w++ {
				h[i][0].Coeffs[index][w] = ring.CRed(h[i][0].Coeffs[index][w]+ekg.polypool.Coeffs[index][w], contextKeys.Modulus[index])
			}

			// Handles the case where nb pj does not divides nb qi
			if index == uint64(len(ekg.ckksContext.ContextQ().Modulus)) {
				break
			}
		}
	}

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {
		// u
		ekg.ternarySampler.SampleMontgomeryNTT(0.5, ekg.polypool)
		// h_0 = pk_0 * u + e0 + P * sk * (qiBarre*qiStar)%qi
		contextKeys.MulCoeffsMontgomeryAndAdd(pk[0], ekg.polypool, h[i][0])
		// h_1 = pk_1 * u + e1 + P * sk * (qiBarre*qiStar)%qi
		contextKeys.MulCoeffsMontgomeryAndAdd(pk[1], ekg.polypool, h[i][1])
	}

	ekg.polypool.Zero()

	return h
}

// Aggregate is the first part of the second and last round of the naive EKG protocol. Uppon receiving the j-1 elements, each party computes :
//
// [sum(cpk[0] * u_j + P * s_j + e_0j), sum(cpk[1] * u_j + e_1j)]
//
// = [cpk[0] * u + P * s + e_0, cpk[1]*u + e_1]
//
// Using this intermediate result, each party computes :
//
// [s_i * (cpk[0] * u + P * s + e_0) + v_i * cpk[0] + e_2i, s_i * (cpk[1] * u + e_1) + cpk[1] * v_i + e_3i]
//
// = [ cpk[0] * (u * s_i) + P * (s * s_i) + (s_i * e_0) + v_i*cpk[0] + e_2i, cpk[1] * u * s_i + (s_i * e_1) + cpk[1] * v_i + e_3i]
//
// And party broadcast this last result to the other j-1 parties.
func (ekg *EkgProtocolNaive) Aggregate(sk *ring.Poly, pk [2]*ring.Poly, samples [][][2]*ring.Poly) (h [][2]*ring.Poly) {

	contextKeys := ekg.ckksContext.ContextKeys()

	h = make([][2]*ring.Poly, ekg.ckksContext.Beta())

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {

		h[i][0] = samples[0][i][0].CopyNew()
		h[i][1] = samples[0][i][1].CopyNew()

		// h_0 = sum(samples[0])
		// h_1 = sum(samples[1])
		for j := 1; j < len(samples); j++ {
			contextKeys.AddNoMod(h[i][0], samples[j][i][0], h[i][0])
			contextKeys.AddNoMod(h[i][1], samples[j][i][1], h[i][1])

			if j&7 == 7 {
				contextKeys.Reduce(h[i][0], h[i][0])
				contextKeys.Reduce(h[i][1], h[i][1])
			}
		}

		if (len(samples)-1)&7 == 7 {
			contextKeys.Reduce(h[i][0], h[i][0])
			contextKeys.Reduce(h[i][1], h[i][1])
		}

		// h_0 = sum(samples[0]) * sk
		// h_1 = sum(samples[1]) * sk
		contextKeys.MulCoeffsMontgomery(h[i][0], sk, h[i][0])
		contextKeys.MulCoeffsMontgomery(h[i][1], sk, h[i][1])

		// v
		ekg.ternarySampler.SampleMontgomeryNTT(0.5, ekg.polypool)

		// h_0 = sum(samples[0]) * sk + pk0 * v
		contextKeys.MulCoeffsMontgomeryAndAdd(pk[0], ekg.polypool, h[i][0])

		// h_1 = sum(samples[1]) * sk + pk1 * v
		contextKeys.MulCoeffsMontgomeryAndAdd(pk[1], ekg.polypool, h[i][1])

		// h_0 = sum(samples[0]) * sk + pk0 * v + e2
		ekg.gaussianSampler.SampleNTT(ekg.polypool)
		contextKeys.Add(h[i][0], ekg.polypool, h[i][0])

		// h_1 = sum(samples[1]) * sk + pk1 * v + e3
		ekg.gaussianSampler.SampleNTT(ekg.polypool)
		contextKeys.Add(h[i][1], ekg.polypool, h[i][1])

	}

	ekg.polypool.Zero()

	return h
}

// Finalize is the second part of the second and last round of the naive EKG protocol. Uppon receiving the j-1 elements,
// each party computes :
//
// [ sum(cpk[0] * (u*s_i) + P * (s*s_i) + (s_i*e_0) + v_i*cpk[0] + e_2i), sum(cpk[1]*u*s_i + (s_i*e_1) + cpk[1] * v_i + e_3i)]
//
// = [cpk[0] * (s*u + v) + (P * s^2 ) + s*e_0 + e_2, ckp[1] * (s*u + v) + s*e_1 + e_3]
//
// = [-s*b + P * s^2 - (s*u + b) * e_cpk + s*e_0 + e_2, b + s*e_1 + e_3]
func (ekg *EkgProtocolNaive) Finalize(h [][][2]*ring.Poly) (evaluationKey [][2]*ring.Poly) {

	contextKeys := ekg.ckksContext.ContextKeys()

	evaluationKey = make([][2]*ring.Poly, ekg.ckksContext.Beta())

	for i := uint64(0); i < ekg.ckksContext.Beta(); i++ {

		evaluationKey[i][0] = h[0][i][0].CopyNew()
		evaluationKey[i][1] = h[0][i][1].CopyNew()

		for j := 1; j < len(h); j++ {
			contextKeys.AddNoMod(evaluationKey[i][0], h[j][i][0], evaluationKey[i][0])
			contextKeys.AddNoMod(evaluationKey[i][1], h[j][i][1], evaluationKey[i][1])

			if j&7 == 7 {
				contextKeys.Reduce(evaluationKey[i][0], evaluationKey[i][0])
				contextKeys.Reduce(evaluationKey[i][1], evaluationKey[i][1])
			}
		}

		if (len(h)-1)&7 == 7 {
			contextKeys.Reduce(evaluationKey[i][0], evaluationKey[i][0])
			contextKeys.Reduce(evaluationKey[i][1], evaluationKey[i][1])
		}

		contextKeys.MForm(evaluationKey[i][0], evaluationKey[i][0])
		contextKeys.MForm(evaluationKey[i][1], evaluationKey[i][1])

	}

	return evaluationKey
}
