package dckks

import (
	"github.com/ldsec/lattigo/ring"
	"math"
)

// EkgProtocol is a structure storing the parameters for the collective evaluation-key generation.
type EkgProtocol struct {
	context         *ring.Context
	keyswitchprimes []uint64
	alpha           uint64
	beta            uint64
	ternarySampler  *ring.TernarySampler
	gaussianSampler *ring.KYSampler
	polypool        *ring.Poly
}

// NewEkgProtocol creates a new EkgProtocol object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewEkgProtocol(context *ring.Context, keyswitchprimes []uint64) *EkgProtocol {
	ekg := new(EkgProtocol)
	ekg.context = context

	ekg.keyswitchprimes = make([]uint64, len(keyswitchprimes))
	for i := range keyswitchprimes {
		ekg.keyswitchprimes[i] = keyswitchprimes[i]
	}

	ekg.alpha = uint64(len(keyswitchprimes))
	ekg.beta = uint64(math.Ceil(float64(len(context.Modulus)-len(keyswitchprimes)) / float64(ekg.alpha)))

	ekg.ternarySampler = context.NewTernarySampler()
	ekg.gaussianSampler = context.NewKYSampler(3.19, 19)
	ekg.polypool = context.NewPoly()
	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first round).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *EkgProtocol) NewEphemeralKey(p float64) (ephemeralKey *ring.Poly, err error) {
	if ephemeralKey, err = ekg.ternarySampler.SampleMontgomeryNTTNew(p); err != nil {
		return nil, err
	}
	return
}

// GenSamples is the first of three rounds of the EkgProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + P*s_i + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *EkgProtocol) GenSamples(u, sk *ring.Poly, crp []*ring.Poly) (h []*ring.Poly) {

	h = make([]*ring.Poly, ekg.beta)

	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	ekg.polypool.Copy(sk)

	for _, pj := range ekg.keyswitchprimes {
		ekg.context.MulScalar(ekg.polypool, pj, ekg.polypool)
	}

	ekg.context.InvMForm(ekg.polypool, ekg.polypool)

	for i := uint64(0); i < ekg.beta; i++ {

		// h = e
		h[i] = ekg.gaussianSampler.SampleNTTNew()

		// h = sk*CrtBaseDecompQi + e
		for j := uint64(0); j < ekg.alpha; j++ {

			for w := uint64(0); w < ekg.context.N; w++ {
				h[i].Coeffs[i*ekg.alpha+j][w] = ring.CRed(h[i].Coeffs[i*ekg.alpha+j][w]+ekg.polypool.Coeffs[i*ekg.alpha+j][w], ekg.context.Modulus[i*ekg.alpha+j])
			}

			// Handles the case where nb pj does not divides nb qi
			if i*ekg.alpha+j == uint64(len(ekg.context.Modulus)-len(ekg.keyswitchprimes)-1) {
				break
			}
		}

		// h = sk*CrtBaseDecompQi + -u*a + e
		ekg.context.MulCoeffsMontgomeryAndSub(u, crp[i], h[i])
	}

	ekg.polypool.Zero()

	return
}

// Aggregate is the second of three rounds of the EkgProtocol protocol. Uppon received the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + P*s_j + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + P*s + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *EkgProtocol) Aggregate(sk *ring.Poly, samples [][]*ring.Poly, crp []*ring.Poly) (h [][2]*ring.Poly) {

	h = make([][2]*ring.Poly, ekg.beta)

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := uint64(0); i < ekg.beta; i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// First Element
		h[i][0] = samples[0][i].CopyNew()

		// Continues with the sum samples
		for j := 1; j < len(samples); j++ {
			ekg.context.AddNoMod(h[i][0], samples[j][i], h[i][0])

			if j&7 == 7 {
				ekg.context.Reduce(h[i][0], h[i][0])
			}
		}

		if (len(samples)-1)&7 != 7 {
			ekg.context.Reduce(h[i][0], h[i][0])
		}

		// (Sum samples) * sk
		ekg.context.MulCoeffsMontgomery(h[i][0], sk, h[i][0])

		// (Sum samples) * sk + e_1i
		ekg.gaussianSampler.SampleNTT(ekg.polypool)
		ekg.context.Add(h[i][0], ekg.polypool, h[i][0])

		// Second Element

		// e_2i
		h[i][1] = ekg.gaussianSampler.SampleNTTNew()
		// s*a + e_2i
		ekg.context.MulCoeffsMontgomeryAndAdd(sk, crp[i], h[i][1])

	}

	ekg.polypool.Zero()

	return
}

// Sum is the first part of the third and last round of the EkgProtocol protocol. Uppon receiving the j-1 elements, each party
// computues :
//
// [sum(s_j * (-u*a + P*s + e) + e_j1), sum(s_j*a + e_j2)]
//
// = [s * (-u*a + P*s + e) + e_1, s*a + e_2].
func (ekg *EkgProtocol) Sum(samples [][][2]*ring.Poly) (h [][2]*ring.Poly) {

	h = make([][2]*ring.Poly, ekg.beta)

	for i := uint64(0); i < ekg.beta; i++ {

		h[i][0] = samples[0][i][0].CopyNew()
		h[i][1] = samples[0][i][1].CopyNew()

		for j := 1; j < len(samples); j++ {
			ekg.context.AddNoMod(h[i][0], samples[j][i][0], h[i][0])
			ekg.context.AddNoMod(h[i][1], samples[j][i][1], h[i][1])

			if j&7 == 7 {
				ekg.context.Reduce(h[i][0], h[i][0])
				ekg.context.Reduce(h[i][1], h[i][1])
			}
		}
		if (len(samples)-1)&7 != 7 {
			ekg.context.Reduce(h[i][0], h[i][0])
			ekg.context.Reduce(h[i][1], h[i][1])
		}

	}

	return
}

// KeySwitch is the second pard of the third and last round of the EkgProtocol protocol. Each party operates a key-switch on [s*a + e_2],
// by computing :
//
// [(u_i - s_i)*(s*a + e_2)]
//
// and broadcasts the result the other j-1 parties.
func (ekg *EkgProtocol) KeySwitch(u, sk *ring.Poly, samples [][2]*ring.Poly) (h1 []*ring.Poly) {

	h1 = make([]*ring.Poly, ekg.beta)

	// (u_i - s_i)
	mask := ekg.context.NewPoly()
	ekg.context.Sub(u, sk, mask)

	for i := uint64(0); i < ekg.beta; i++ {

		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		h1[i] = ekg.gaussianSampler.SampleNTTNew()
		ekg.context.MulCoeffsMontgomeryAndAdd(mask, samples[i][1], h1[i])

	}

	return
}

// ComputeEVK is third part ot the third and last round of the EkgProtocol protocol. Uppon receiving the other j-1 elements, each party computes :
//
// [s * (-u*a + P*s + e) + e_1 + sum([(u_j - s_j)*(s*a + e_2)])]
//
// = [s * (-u*a + P*s + e) + e_1 + (u - s)*(s*a + e_2)]
//
// = [-s*u*a + P*s^2 + s*e + e_1 + s*u*a -s^2*a + (u - s)*e_2]
//
// = [-s^2*a + P*s^2 + e_1 + (u - s)*e_2]
//
// = [-s^2*a + P*s^2 + e]
//
// The evaluation key is therefor : [-s*b + P*s^2 + e, s*b]
func (ekg *EkgProtocol) ComputeEVK(h1 [][]*ring.Poly, h [][2]*ring.Poly) (collectiveEVK [][2]*ring.Poly) {

	collectiveEVK = make([][2]*ring.Poly, ekg.beta)

	// collectiveEVK[i][0] = h[i][0] + sum(h1[i])
	// collectiveEVK[i][1] = h[i][1]
	for i := uint64(0); i < ekg.beta; i++ {

		collectiveEVK[i][0] = h[i][0].CopyNew()
		collectiveEVK[i][1] = h[i][1].CopyNew()

		for j := range h1 {
			ekg.context.AddNoMod(collectiveEVK[i][0], h1[j][i], collectiveEVK[i][0])

			if j&7 == 7 {
				ekg.context.Reduce(collectiveEVK[i][0], collectiveEVK[i][0])
			}
		}

		if (len(h1)-1)&7 != 7 {
			ekg.context.Reduce(collectiveEVK[i][0], collectiveEVK[i][0])
		}

		ekg.context.MForm(collectiveEVK[i][0], collectiveEVK[i][0])
		ekg.context.MForm(collectiveEVK[i][1], collectiveEVK[i][1])
	}

	return
}
