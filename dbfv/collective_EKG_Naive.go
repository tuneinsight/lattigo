package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"math"
)

type EkgProtocolNaive struct {
	context         *ring.Context
	gaussianSampler *ring.KYSampler
	ternarySampler  *ring.TernarySampler
	bitDecomp       uint64
	bitLog          uint64
	polypool        *ring.Poly
}

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
				h[i][w][0].Coeffs[i][j] += bfv.PowerOf2(sk.Coeffs[i][j], ekg.bitDecomp*w, qi, mredParams[i])
			}

			// u
			ekg.ternarySampler.SampleMontgomeryNTT(ekg.polypool)

			// h_0 = pk_0 * u + e0 + sk*w*(qiBarre*qiStar)%qi
			ekg.context.MulCoeffsMontgomeryAndAdd(pk[0], ekg.polypool, h[i][w][0])
			// h_1 = pk_1 * u + e1 + sk*w*(qiBarre*qiStar)%qi
			ekg.context.MulCoeffsMontgomeryAndAdd(pk[1], ekg.polypool, h[i][w][1])
		}
	}

	ekg.polypool.Zero()

	return

}

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
			ekg.ternarySampler.SampleMontgomeryNTT(ekg.polypool)

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
