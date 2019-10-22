package dckks

import (
	"github.com/ldsec/lattigo/ring"
	"math"
)

// EkgProtocolNaive is a structure storing the parameters for the naive EKG protocol.
type EkgProtocolNaive struct {
	context         *ring.Context
	keyswitchprimes []uint64
	alpha           uint64
	beta            uint64
	gaussianSampler *ring.KYSampler
	ternarySampler  *ring.TernarySampler
	polypool        *ring.Poly
}

// NewEkgProtocolNaive creates a new EkgProtocolNaive object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewEkgProtocolNaive(context *ring.Context, keyswitchprimes []uint64) (ekg *EkgProtocolNaive) {
	ekg = new(EkgProtocolNaive)
	ekg.context = context

	ekg.keyswitchprimes = make([]uint64, len(keyswitchprimes))
	for i := range keyswitchprimes {
		ekg.keyswitchprimes[i] = keyswitchprimes[i]
	}

	ekg.alpha = uint64(len(keyswitchprimes))
	ekg.beta = uint64(math.Ceil(float64(len(context.Modulus)-len(keyswitchprimes)) / float64(ekg.alpha)))

	ekg.gaussianSampler = context.NewKYSampler(3.19, 19)
	ekg.ternarySampler = context.NewTernarySampler()
	ekg.polypool = context.NewPoly()
	return ekg
}

// GenSamples is the first of two rounds of the naive EKG protocol. Using the shared public key "cpk",
// each party generates a pseudo-encryption of s*w of the form :
//
// [cpk[0] * u_i + P * s_i + e_0i, cpk[1] * u_i + e_1i]
//
// and broadcasts it to all other j-1 parties.
func (ekg *EkgProtocolNaive) GenSamples(sk *ring.Poly, pk [2]*ring.Poly) (h [][2]*ring.Poly) {

	h = make([][2]*ring.Poly, ekg.beta)

	ekg.polypool.Copy(sk)

	for _, pj := range ekg.keyswitchprimes {
		ekg.context.MulScalar(ekg.polypool, pj, ekg.polypool)
	}

	ekg.context.InvMForm(ekg.polypool, ekg.polypool)

	for i := uint64(0); i < ekg.beta; i++ {
		// h_0 = e0
		h[i][0] = ekg.gaussianSampler.SampleNTTNew()
		// h_1 = e1
		h[i][1] = ekg.gaussianSampler.SampleNTTNew()

		// h_0 = e0 + [sk*P*(qiBarre*qiStar)%qi = sk*P, else 0]

		for j := uint64(0); j < ekg.alpha; j++ {

			for w := uint64(0); w < ekg.context.N; w++ {
				h[i][0].Coeffs[i*ekg.alpha+j][w] = ring.CRed(h[i][0].Coeffs[i*ekg.alpha+j][w]+ekg.polypool.Coeffs[i*ekg.alpha+j][w], ekg.context.Modulus[i*ekg.alpha+j])
			}

			// Handles the case where nb pj does not divides nb qi
			if i*ekg.alpha+j == uint64(len(ekg.context.Modulus)-len(ekg.keyswitchprimes)-1) {
				break
			}
		}
	}

	for i := uint64(0); i < ekg.beta; i++ {
		// u
		ekg.ternarySampler.SampleMontgomeryNTT(0.5, ekg.polypool)
		// h_0 = pk_0 * u + e0 + P * sk * (qiBarre*qiStar)%qi
		ekg.context.MulCoeffsMontgomeryAndAdd(pk[0], ekg.polypool, h[i][0])
		// h_1 = pk_1 * u + e1 + P * sk * (qiBarre*qiStar)%qi
		ekg.context.MulCoeffsMontgomeryAndAdd(pk[1], ekg.polypool, h[i][1])
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

	h = make([][2]*ring.Poly, ekg.beta)

	for i := uint64(0); i < ekg.beta; i++ {

		h[i][0] = samples[0][i][0].CopyNew()
		h[i][1] = samples[0][i][1].CopyNew()

		// h_0 = sum(samples[0])
		// h_1 = sum(samples[1])
		for j := 1; j < len(samples); j++ {
			ekg.context.AddNoMod(h[i][0], samples[j][i][0], h[i][0])
			ekg.context.AddNoMod(h[i][1], samples[j][i][1], h[i][1])

			if j&7 == 7 {
				ekg.context.Reduce(h[i][0], h[i][0])
				ekg.context.Reduce(h[i][1], h[i][1])
			}
		}

		if (len(samples)-1)&7 == 7 {
			ekg.context.Reduce(h[i][0], h[i][0])
			ekg.context.Reduce(h[i][1], h[i][1])
		}

		// h_0 = sum(samples[0]) * sk
		// h_1 = sum(samples[1]) * sk
		ekg.context.MulCoeffsMontgomery(h[i][0], sk, h[i][0])
		ekg.context.MulCoeffsMontgomery(h[i][1], sk, h[i][1])

		// v
		ekg.ternarySampler.SampleMontgomeryNTT(0.5, ekg.polypool)

		// h_0 = sum(samples[0]) * sk + pk0 * v
		ekg.context.MulCoeffsMontgomeryAndAdd(pk[0], ekg.polypool, h[i][0])

		// h_1 = sum(samples[1]) * sk + pk1 * v
		ekg.context.MulCoeffsMontgomeryAndAdd(pk[1], ekg.polypool, h[i][1])

		// h_0 = sum(samples[0]) * sk + pk0 * v + e2
		ekg.gaussianSampler.SampleNTT(ekg.polypool)
		ekg.context.Add(h[i][0], ekg.polypool, h[i][0])

		// h_1 = sum(samples[1]) * sk + pk1 * v + e3
		ekg.gaussianSampler.SampleNTT(ekg.polypool)
		ekg.context.Add(h[i][1], ekg.polypool, h[i][1])

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

	evaluationKey = make([][2]*ring.Poly, ekg.beta)

	for i := uint64(0); i < ekg.beta; i++ {

		evaluationKey[i][0] = h[0][i][0].CopyNew()
		evaluationKey[i][1] = h[0][i][1].CopyNew()

		for j := 1; j < len(h); j++ {
			ekg.context.AddNoMod(evaluationKey[i][0], h[j][i][0], evaluationKey[i][0])
			ekg.context.AddNoMod(evaluationKey[i][1], h[j][i][1], evaluationKey[i][1])

			if j&7 == 7 {
				ekg.context.Reduce(evaluationKey[i][0], evaluationKey[i][0])
				ekg.context.Reduce(evaluationKey[i][1], evaluationKey[i][1])
			}
		}

		if (len(h)-1)&7 == 7 {
			ekg.context.Reduce(evaluationKey[i][0], evaluationKey[i][0])
			ekg.context.Reduce(evaluationKey[i][1], evaluationKey[i][1])
		}

		ekg.context.MForm(evaluationKey[i][0], evaluationKey[i][0])
		ekg.context.MForm(evaluationKey[i][1], evaluationKey[i][1])

	}

	return evaluationKey
}
