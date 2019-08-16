package dbfv

import (
	"github.com/lca1/lattigo/bfv"
	"github.com/lca1/lattigo/ring"
	"math"
)

type EkgProtocol struct {
	context         *ring.Context
	gaussianSampler *ring.KYSampler
	bitDecomp       uint64
	bitLog          uint64
	polypool        *ring.Poly
}

func NewEkgProtocol(context *ring.Context, bitDecomp uint64) *EkgProtocol {
	ekg := new(EkgProtocol)
	ekg.context = context
	ekg.gaussianSampler = context.NewKYSampler(3.19, 19)
	ekg.bitDecomp = bitDecomp
	ekg.bitLog = uint64(math.Ceil(float64(60) / float64(bitDecomp)))
	ekg.polypool = context.NewPoly()
	return ekg
}

// Ephemeral Key u (needs to be stored among the 3 first round)
func (ekg *EkgProtocol) NewEphemeralKey() *ring.Poly {
	ephemeralKey := ekg.context.NewTernaryPoly()
	ekg.context.NTT(ephemeralKey, ephemeralKey)
	ekg.context.MForm(ephemeralKey, ephemeralKey)
	return ephemeralKey
}

// Ephemeral Key u (needs to be stored among the 3 first round)
func (ekg *EkgProtocol) GenSamples(u, sk *ring.Poly, crp [][]*ring.Poly) (h [][]*ring.Poly) {

	h = make([][]*ring.Poly, len(ekg.context.Modulus))

	mredParams := ekg.context.GetMredParams()

	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + s*w_i + e_i]
	// where a_i = crp_i
	for i, qi := range ekg.context.Modulus {

		h[i] = make([]*ring.Poly, ekg.bitLog)

		for w := uint64(0); w < ekg.bitLog; w++ {

			// h = e
			h[i][w] = ekg.gaussianSampler.SampleNTTNew()

			// h = sk*CrtBaseDecompQi + e
			for j := uint64(0); j < ekg.context.N; j++ {
				h[i][w].Coeffs[i][j] += bfv.PowerOf2(sk.Coeffs[i][j], ekg.bitDecomp*w, qi, mredParams[i])
			}

			// h = sk*CrtBaseDecompQi + -u*a + e
			ekg.context.MulCoeffsMontgomeryAndSub(u, crp[i][w], h[i][w])
		}
	}

	return
}

// Round 2
func (ekg *EkgProtocol) Aggregate(sk *ring.Poly, samples [][][]*ring.Poly, crp [][]*ring.Poly) (h [][][2]*ring.Poly) {

	h = make([][][2]*ring.Poly, len(ekg.context.Modulus))

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := range ekg.context.Modulus {

		h[i] = make([][2]*ring.Poly, ekg.bitLog)

		for w := uint64(0); w < ekg.bitLog; w++ {

			// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

			// First Element
			h[i][w][0] = samples[0][i][w].CopyNew()

			// Continues with the sum samples
			for j := 1; j < len(samples); j++ {
				ekg.context.AddNoMod(h[i][w][0], samples[j][i][w], h[i][w][0])

				if j&7 == 7 {
					ekg.context.Reduce(h[i][w][0], h[i][w][0])
				}
			}

			if (len(samples)-1)&7 != 7 {
				ekg.context.Reduce(h[i][w][0], h[i][w][0])
			}

			// (Sum samples) * sk
			ekg.context.MulCoeffsMontgomery(h[i][w][0], sk, h[i][w][0])

			// (Sum samples) * sk + e_1i
			ekg.gaussianSampler.SampleNTT(ekg.polypool)
			ekg.context.Add(h[i][w][0], ekg.polypool, h[i][w][0])

			// Second Element

			// e_2i
			h[i][w][1] = ekg.gaussianSampler.SampleNTTNew()
			// s*a + e_2i
			ekg.context.MulCoeffsMontgomeryAndAdd(sk, crp[i][w], h[i][w][1])

		}

	}

	ekg.polypool.Zero()

	return
}

// Round 3

//Part 1
func (ekg *EkgProtocol) Sum(samples [][][][2]*ring.Poly) (h [][][2]*ring.Poly) {

	h = make([][][2]*ring.Poly, len(ekg.context.Modulus))

	for i := range ekg.context.Modulus {

		h[i] = make([][2]*ring.Poly, ekg.bitLog)

		for w := uint64(0); w < ekg.bitLog; w++ {

			h[i][w][0] = samples[0][i][w][0].CopyNew()
			h[i][w][1] = samples[0][i][w][1].CopyNew()

			for j := 1; j < len(samples); j++ {
				ekg.context.AddNoMod(h[i][w][0], samples[j][i][w][0], h[i][w][0])
				ekg.context.AddNoMod(h[i][w][1], samples[j][i][w][1], h[i][w][1])

				if j&7 == 7 {
					ekg.context.Reduce(h[i][w][0], h[i][w][0])
					ekg.context.Reduce(h[i][w][1], h[i][w][1])
				}
			}
			if (len(samples)-1)&7 != 7 {
				ekg.context.Reduce(h[i][w][0], h[i][w][0])
				ekg.context.Reduce(h[i][w][1], h[i][w][1])
			}
		}
	}

	return
}

// Part 2
func (ekg *EkgProtocol) KeySwitch(u, sk *ring.Poly, samples [][][2]*ring.Poly) (h1 [][]*ring.Poly) {

	h1 = make([][]*ring.Poly, len(ekg.context.Modulus))

	// (u_i - s_i)
	mask := ekg.context.NewPoly()
	ekg.context.Sub(u, sk, mask)

	for i := range ekg.context.Modulus {

		h1[i] = make([]*ring.Poly, ekg.bitLog)

		for w := uint64(0); w < ekg.bitLog; w++ {

			// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
			h1[i][w] = ekg.gaussianSampler.SampleNTTNew()
			ekg.context.MulCoeffsMontgomeryAndAdd(mask, samples[i][w][1], h1[i][w])
		}
	}

	return h1
}

// Round 4

func (ekg *EkgProtocol) ComputeEVK(h1 [][][]*ring.Poly, h [][][2]*ring.Poly) (collectiveEVK [][][2]*ring.Poly) {

	collectiveEVK = make([][][2]*ring.Poly, len(ekg.context.Modulus))

	// collectiveEVK[i][0] = h[i][0] + sum(h1[i])
	// collectiveEVK[i][1] = h[i][1]
	for i := range ekg.context.Modulus {

		collectiveEVK[i] = make([][2]*ring.Poly, ekg.bitLog)

		for w := uint64(0); w < ekg.bitLog; w++ {

			collectiveEVK[i][w][0] = h[i][w][0].CopyNew()
			collectiveEVK[i][w][1] = h[i][w][1].CopyNew()

			for j := range h1 {
				ekg.context.AddNoMod(collectiveEVK[i][w][0], h1[j][i][w], collectiveEVK[i][w][0])

				if j&7 == 7 {
					ekg.context.Reduce(collectiveEVK[i][w][0], collectiveEVK[i][w][0])
				}
			}

			if (len(h1)-1)&7 != 7 {
				ekg.context.Reduce(collectiveEVK[i][w][0], collectiveEVK[i][w][0])
			}

			ekg.context.MForm(collectiveEVK[i][w][0], collectiveEVK[i][w][0])
			ekg.context.MForm(collectiveEVK[i][w][1], collectiveEVK[i][w][1])
		}
	}

	return
}
