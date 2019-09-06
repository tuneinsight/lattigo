package dbfv

import (
	"github.com/lca1/lattigo/ring"
	"math"
)

// EkgProtocolNaive is a structure storing the parameters for the naive EKG protocol.
type EkgProtocolNaive struct {
	context         *ring.Context
	gaussianSampler *ring.KYSampler
	ternarySampler  *ring.TernarySampler
	bitDecomp       uint64
	bitLog          uint64
	polypool        *ring.Poly
}

// NewEkgProtocolNaive creates a new EkgProtocolNaive object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewEkgProtocolNaive(context *ring.Context, bitDecomp uint64) *EkgProtocolNaive {
	ekg := new(EkgProtocolNaive)
	ekg.context = context
	ekg.gaussianSampler = context.NewKYSampler(3.19, 19)
	ekg.ternarySampler = context.NewTernarySampler()
	ekg.bitDecomp = bitDecomp
	ekg.bitLog = uint64(math.Ceil(float64(60) / float64(bitDecomp)))
	ekg.polypool = context.NewPoly()
	return ekg
}

// GenSamples is the first of two rounds of the naive EKG protocol. Using the shared public key "cpk",
// each party generates a pseudo-encryption of s*w of the form :
//
// [cpk[0]*u_i + s_i * w + e_0i, cpk[1]*u_i + e_1i]
//
// and broadcasts it to all other j-1 parties.
func (ekg *EkgProtocolNaive) GenSamples(sk *ring.Poly, pk [2]*ring.Poly) (h [][][2]*ring.Poly) {

	h = make([][][2]*ring.Poly, len(ekg.context.Modulus))

	mredParams := ekg.context.GetMredParams()

	for i, qi := range ekg.context.Modulus {

		h[i] = make([][2]*ring.Poly, ekg.bitLog)

		for w := uint64(0); w < ekg.bitLog; w++ {

			// h_0 = e0
			h[i][w][0] = ekg.gaussianSampler.SampleNTTNew()
			// h_1 = e1
			h[i][w][1] = ekg.gaussianSampler.SampleNTTNew()

			// h_0 = e0 + [sk*w*(qiBarre*qiStar)%qi = 1<<w, else 0]
			for j := uint64(0); j < ekg.context.N; j++ {
				h[i][w][0].Coeffs[i][j] += ring.PowerOf2(sk.Coeffs[i][j], ekg.bitDecomp*w, qi, mredParams[i])
			}

			// u
			ekg.ternarySampler.SampleMontgomeryNTT(0.5, ekg.polypool)

			// h_0 = pk_0 * u + e0 + sk*w*(qiBarre*qiStar)%qi
			ekg.context.MulCoeffsMontgomeryAndAdd(pk[0], ekg.polypool, h[i][w][0])
			// h_1 = pk_1 * u + e1 + sk*w*(qiBarre*qiStar)%qi
			ekg.context.MulCoeffsMontgomeryAndAdd(pk[1], ekg.polypool, h[i][w][1])
		}
	}

	ekg.polypool.Zero()

	return

}

// Aggregate is the first part of the second and last round of the naive EKG protocol. Uppon receiving the j-1 elements, each party computes :
//
// [sum(cpk[0]*u_j + s_j * w + e_0j), sum(cpk[1]*u_j + e_1j)]
//
// = [cpk[0]*u + s * w + e_0, cpk[1]*u + e_1]
//
// Using this intermediate result, each party computes :
//
// [s_i * (cpk[0]*u + s * w + e_0) + v_i*cpk[0] + e_2i, s_i*(cpk[1]*u + e_1) + cpk[1] * v_i + e_3i]
//
// = [ cpk[0] * (u*s_i) + (s*s_i) * w + (s_i*e_0) + v_i*cpk[0] + e_2i, cpk[1]*u*s_i + (s_i*e_1) + cpk[1] * v_i + e_3i]
//
// And party broadcast this last result to the other j-1 parties.
func (ekg *EkgProtocolNaive) Aggregate(sk *ring.Poly, pk [2]*ring.Poly, samples [][][][2]*ring.Poly) (h [][][2]*ring.Poly) {

	h = make([][][2]*ring.Poly, len(ekg.context.Modulus))

	for i := range h {

		h[i] = make([][2]*ring.Poly, ekg.bitLog)

		for w := uint64(0); w < ekg.bitLog; w++ {

			h[i][w][0] = samples[0][i][w][0].CopyNew()
			h[i][w][1] = samples[0][i][w][1].CopyNew()

			// h_0 = sum(samples[0])
			// h_1 = sum(samples[1])
			for j := 1; j < len(samples); j++ {
				ekg.context.AddNoMod(h[i][w][0], samples[j][i][w][0], h[i][w][0])
				ekg.context.AddNoMod(h[i][w][1], samples[j][i][w][1], h[i][w][1])

				if j&7 == 7 {
					ekg.context.Reduce(h[i][w][0], h[i][w][0])
					ekg.context.Reduce(h[i][w][1], h[i][w][1])
				}
			}

			if (len(samples)-1)&7 == 7 {
				ekg.context.Reduce(h[i][w][0], h[i][w][0])
				ekg.context.Reduce(h[i][w][1], h[i][w][1])
			}

			// h_0 = sum(samples[0]) * sk
			// h_1 = sum(samples[1]) * sk
			ekg.context.MulCoeffsMontgomery(h[i][w][0], sk, h[i][w][0])
			ekg.context.MulCoeffsMontgomery(h[i][w][1], sk, h[i][w][1])

			// v
			ekg.ternarySampler.SampleMontgomeryNTT(0.5, ekg.polypool)

			// h_0 = sum(samples[0]) * sk + pk0 * v
			ekg.context.MulCoeffsMontgomeryAndAdd(pk[0], ekg.polypool, h[i][w][0])

			// h_1 = sum(samples[1]) * sk + pk1 * v
			ekg.context.MulCoeffsMontgomeryAndAdd(pk[1], ekg.polypool, h[i][w][1])

			// h_0 = sum(samples[0]) * sk + pk0 * v + e2
			ekg.gaussianSampler.SampleNTT(ekg.polypool)
			ekg.context.Add(h[i][w][0], ekg.polypool, h[i][w][0])

			// h_1 = sum(samples[1]) * sk + pk1 * v + e3
			ekg.gaussianSampler.SampleNTT(ekg.polypool)
			ekg.context.Add(h[i][w][1], ekg.polypool, h[i][w][1])
		}
	}

	ekg.polypool.Zero()

	return
}

// Finalize is the second part of the second and last round of the naive EKG protocol. Uppon receiving the j-1 elements,
// each party computes :
//
// [ sum(cpk[0] * (u*s_i) + (s*s_i) * w + (s_i*e_0) + v_i*cpk[0] + e_2i), sum(cpk[1]*u*s_i + (s_i*e_1) + cpk[1] * v_i + e_3i)]
//
// = [cpk[0] * (s*u + v) + (s^2 * w) + s*e_0 + e_2, ckp[1] * (s*u + v) + s*e_1 + e_3]
//
// = [-s*b + s^2 * w - (s*u + b) * e_cpk + s*e_0 + e_2, b + s*e_1 + e_3]
func (ekg *EkgProtocolNaive) Finalize(h [][][][2]*ring.Poly) (evaluationKey [][][2]*ring.Poly) {

	evaluationKey = make([][][2]*ring.Poly, len(ekg.context.Modulus))

	for i := range ekg.context.Modulus {

		evaluationKey[i] = make([][2]*ring.Poly, ekg.bitLog)

		for w := uint64(0); w < ekg.bitLog; w++ {

			evaluationKey[i][w][0] = h[0][i][w][0].CopyNew()
			evaluationKey[i][w][1] = h[0][i][w][1].CopyNew()

			for j := 1; j < len(h); j++ {
				ekg.context.AddNoMod(evaluationKey[i][w][0], h[j][i][w][0], evaluationKey[i][w][0])
				ekg.context.AddNoMod(evaluationKey[i][w][1], h[j][i][w][1], evaluationKey[i][w][1])

				if j&7 == 7 {
					ekg.context.Reduce(evaluationKey[i][w][0], evaluationKey[i][w][0])
					ekg.context.Reduce(evaluationKey[i][w][1], evaluationKey[i][w][1])
				}
			}

			if (len(h)-1)&7 == 7 {
				ekg.context.Reduce(evaluationKey[i][w][0], evaluationKey[i][w][0])
				ekg.context.Reduce(evaluationKey[i][w][1], evaluationKey[i][w][1])
			}

			ekg.context.MForm(evaluationKey[i][w][0], evaluationKey[i][w][0])
			ekg.context.MForm(evaluationKey[i][w][1], evaluationKey[i][w][1])
		}
	}

	return
}
