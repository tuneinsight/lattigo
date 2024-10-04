package ring

import (
	"encoding/binary"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

const (
	rn = 3.442619855899
)

// GaussianSampler keeps the state of a truncated Gaussian polynomial sampler.
type GaussianSampler struct {
	*baseSampler
	*randomBuffer
	xe         DiscreteGaussian
	montgomery bool
}

// NewGaussianSampler creates a new instance of GaussianSampler from a PRNG, a ring definition and the truncated
// Gaussian distribution parameters. Sigma is the desired standard deviation and bound is the maximum coefficient norm in absolute
// value.
func NewGaussianSampler(prng sampling.PRNG, baseRing *Ring, X DiscreteGaussian, montgomery bool) (g *GaussianSampler) {
	g = new(GaussianSampler)
	g.baseSampler = &baseSampler{}
	g.prng = prng
	g.randomBuffer = newRandomBuffer()
	g.baseRing = baseRing
	g.xe = X
	g.montgomery = montgomery
	return
}

// AtLevel returns an instance of the target GaussianSampler that operates at the target level.
// This instance is not thread safe and cannot be used concurrently to the base instance.
func (g *GaussianSampler) AtLevel(level int) Sampler {
	return &GaussianSampler{
		baseSampler:  g.baseSampler.AtLevel(level),
		randomBuffer: g.randomBuffer,
		xe:           g.xe,
		montgomery:   g.montgomery,
	}
}

// Read samples a truncated Gaussian polynomial on "pol" at the maximum level in the default ring, standard deviation and bound.
func (g *GaussianSampler) Read(pol Poly) {
	g.read(pol, func(a, b, c uint64) uint64 {
		return b
	})
}

// ReadNew samples a new truncated Gaussian polynomial at the maximum level in the default ring, standard deviation and bound.
func (g *GaussianSampler) ReadNew() (pol Poly) {
	pol = g.baseRing.NewPoly()
	g.Read(pol)
	return pol
}

// ReadAndAdd samples a truncated Gaussian polynomial at the given level for the receiver's default standard deviation and bound and adds it on "pol".
func (g *GaussianSampler) ReadAndAdd(pol Poly) {
	g.read(pol, func(a, b, c uint64) uint64 {
		return CRed(a+b, c)
	})
}

func (g *GaussianSampler) read(pol Poly, f func(a, b, c uint64) uint64) {
	var norm float64

	var sign uint64

	r := g.baseRing

	level := r.level

	if _, err := g.prng.Read(g.randomBufferN); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}

	moduli := r.ModuliChain()[:level+1]

	bound := g.xe.Bound
	sigma := float64(g.xe.Sigma)

	N := r.N()

	coeffs := pol.Coeffs

	// If the standard deviation is greater than float64 precision
	// and the bound is greater than uint64, we switch to an approximation
	// using arbitrary precision.
	//
	// The approximation of the large norm sampling is done by sampling
	// a uniform value [0, sigma] * ceil(norm) * sign.
	if sigma > 0x20000000000000 && bound > 0xffffffffffffffff {

		sigmaInt := new(big.Int)
		new(big.Float).SetFloat64(sigma).Int(sigmaInt)

		Qi := make([]*big.Int, len(moduli))

		for i, qi := range moduli {
			Qi[i] = bignum.NewInt(qi)
		}

		boundInt := new(big.Int)
		new(big.Float).SetFloat64(bound).Int(boundInt)

		coeff := new(big.Int)

		normInt := new(big.Int)
		normFlo := new(big.Float)
		normIntLowBits := new(big.Int)

		for i := 0; i < N; i++ {

			for {
				// Sample norm with sigma = 1 and sign
				norm, sign = g.normFloat64()

				// Sets normFlo = norm * sigma with precision 53 bits
				// and 0.5 for rounding discretization
				normFlo.SetFloat64(norm*sigma + 0.5)

				// Discretizes to an integer
				normFlo.Int(normInt)

				// Derive the number of zero bits: normInt>>53
				normIntLowBits.Rsh(normInt, 53)

				// Sample in the size of the number of zero bits and adds
				// (normInt + rand(normInt>>53)) * sign
				// This might not be constant time
				if normIntLowBits.Cmp(new(big.Int)) > 0 {
					normInt.Add(normInt, bignum.RandInt(g.prng, normIntLowBits))
				}

				normInt.Mul(normInt, bignum.NewInt(2*int64(sign)-1))

				if normInt.Cmp(boundInt) < 1 {
					break
				}
			}

			for j, qi := range moduli {
				coeffs[j][i] = f(coeffs[j][i], coeff.Mod(normInt, Qi[j]).Uint64(), qi)
			}
		}

	} else {

		var coeffInt uint64

		for i := 0; i < N; i++ {

			for {
				norm, sign = g.normFloat64()

				if v := norm * sigma; v <= bound {
					coeffInt = uint64(v + 0.5) // rounding
					break
				}
			}

			for j, qi := range moduli {
				coeffs[j][i] = f(coeffs[j][i], (coeffInt*sign)|(qi-coeffInt)*(sign^1), qi)
			}
		}
	}

	if g.montgomery {
		g.baseRing.MForm(pol, pol)
	}
}

// NormFloat64 returns a normally distributed float64 in
// the range [-math.MaxFloat64, +math.MaxFloat64], bounds included,
// with standard normal distribution (mean = 0, stddev = 1).
// To produce a different normal distribution, callers can
// adjust the output using:
//
//	sample = NormFloat64() * desiredStdDev + desiredMean
//
// Algorithm adapted from https://golang.org/src/math/rand/normal.go
// to use a secure PRNG instead of math/rand.
func (g *GaussianSampler) normFloat64() (float64, uint64) {

	ptr := g.ptr
	buff := g.randomBufferN
	prng := g.prng
	buffLen := len(buff)

	read := func() {
		if ptr == buffLen {
			if _, err := prng.Read(buff); err != nil {
				// Sanity check, this error should not happen.
				panic(err)
			}
			ptr = 0
		}
	}

	randU32 := func() (x uint32) {
		read()
		x = binary.LittleEndian.Uint32(buff[ptr : ptr+4])
		ptr += 8 // Avoids buffer misalignment
		return
	}

	randF64 := func() (x float64) {
		read()
		x = float64(binary.LittleEndian.Uint64(buff[ptr:ptr+8])&0x1fffffffffffff) / float64(0x1fffffffffffff)
		ptr += 8
		return
	}

	for {

		juint32 := randU32()

		j := int32(juint32 & 0x7fffffff)
		sign := uint64(juint32 >> 31)

		i := j & 0x7F

		x := float64(j) * float64(wn[i])

		// 1 (>99%)
		if uint32(j) < kn[i] {
			g.ptr = ptr
			return x, sign
		}

		// 2 (<1%)
		if i == 0 {

			// This extra work is only required for the base strip.
			for {

				x = -math.Log(randF64()) * (1.0 / rn)
				y := -math.Log(randF64())

				if y+y >= x*x {
					break
				}
			}

			g.ptr = ptr
			return x + rn, sign
		}

		// 3
		if fn[i]+float32(randF64())*(fn[i-1]-fn[i]) < float32(math.Exp(-0.5*x*x)) {
			g.ptr = ptr
			return x, sign
		}
	}
}

var kn = [128]uint32{
	0x76ad2212, 0x0, 0x600f1b53, 0x6ce447a6, 0x725b46a2,
	0x7560051d, 0x774921eb, 0x789a25bd, 0x799045c3, 0x7a4bce5d,
	0x7adf629f, 0x7b5682a6, 0x7bb8a8c6, 0x7c0ae722, 0x7c50cce7,
	0x7c8cec5b, 0x7cc12cd6, 0x7ceefed2, 0x7d177e0b, 0x7d3b8883,
	0x7d5bce6c, 0x7d78dd64, 0x7d932886, 0x7dab0e57, 0x7dc0dd30,
	0x7dd4d688, 0x7de73185, 0x7df81cea, 0x7e07c0a3, 0x7e163efa,
	0x7e23b587, 0x7e303dfd, 0x7e3beec2, 0x7e46db77, 0x7e51155d,
	0x7e5aabb3, 0x7e63abf7, 0x7e6c222c, 0x7e741906, 0x7e7b9a18,
	0x7e82adfa, 0x7e895c63, 0x7e8fac4b, 0x7e95a3fb, 0x7e9b4924,
	0x7ea0a0ef, 0x7ea5b00d, 0x7eaa7ac3, 0x7eaf04f3, 0x7eb3522a,
	0x7eb765a5, 0x7ebb4259, 0x7ebeeafd, 0x7ec2620a, 0x7ec5a9c4,
	0x7ec8c441, 0x7ecbb365, 0x7ece78ed, 0x7ed11671, 0x7ed38d62,
	0x7ed5df12, 0x7ed80cb4, 0x7eda175c, 0x7edc0005, 0x7eddc78e,
	0x7edf6ebf, 0x7ee0f647, 0x7ee25ebe, 0x7ee3a8a9, 0x7ee4d473,
	0x7ee5e276, 0x7ee6d2f5, 0x7ee7a620, 0x7ee85c10, 0x7ee8f4cd,
	0x7ee97047, 0x7ee9ce59, 0x7eea0eca, 0x7eea3147, 0x7eea3568,
	0x7eea1aab, 0x7ee9e071, 0x7ee98602, 0x7ee90a88, 0x7ee86d08,
	0x7ee7ac6a, 0x7ee6c769, 0x7ee5bc9c, 0x7ee48a67, 0x7ee32efc,
	0x7ee1a857, 0x7edff42f, 0x7ede0ffa, 0x7edbf8d9, 0x7ed9ab94,
	0x7ed7248d, 0x7ed45fae, 0x7ed1585c, 0x7ece095f, 0x7eca6ccb,
	0x7ec67be2, 0x7ec22eee, 0x7ebd7d1a, 0x7eb85c35, 0x7eb2c075,
	0x7eac9c20, 0x7ea5df27, 0x7e9e769f, 0x7e964c16, 0x7e8d44ba,
	0x7e834033, 0x7e781728, 0x7e6b9933, 0x7e5d8a1a, 0x7e4d9ded,
	0x7e3b737a, 0x7e268c2f, 0x7e0e3ff5, 0x7df1aa5d, 0x7dcf8c72,
	0x7da61a1e, 0x7d72a0fb, 0x7d30e097, 0x7cd9b4ab, 0x7c600f1a,
	0x7ba90bdc, 0x7a722176, 0x77d664e5,
}

var wn = [128]float32{
	1.7290405e-09, 1.2680929e-10, 1.6897518e-10, 1.9862688e-10,
	2.2232431e-10, 2.4244937e-10, 2.601613e-10, 2.7611988e-10,
	2.9073963e-10, 3.042997e-10, 3.1699796e-10, 3.289802e-10,
	3.4035738e-10, 3.5121603e-10, 3.616251e-10, 3.7164058e-10,
	3.8130857e-10, 3.9066758e-10, 3.9975012e-10, 4.08584e-10,
	4.1719309e-10, 4.2559822e-10, 4.338176e-10, 4.418672e-10,
	4.497613e-10, 4.5751258e-10, 4.651324e-10, 4.7263105e-10,
	4.8001775e-10, 4.87301e-10, 4.944885e-10, 5.015873e-10,
	5.0860405e-10, 5.155446e-10, 5.2241467e-10, 5.2921934e-10,
	5.359635e-10, 5.426517e-10, 5.4928817e-10, 5.5587696e-10,
	5.624219e-10, 5.6892646e-10, 5.753941e-10, 5.818282e-10,
	5.882317e-10, 5.946077e-10, 6.00959e-10, 6.072884e-10,
	6.135985e-10, 6.19892e-10, 6.2617134e-10, 6.3243905e-10,
	6.386974e-10, 6.449488e-10, 6.511956e-10, 6.5744005e-10,
	6.6368433e-10, 6.699307e-10, 6.7618144e-10, 6.824387e-10,
	6.8870465e-10, 6.949815e-10, 7.012715e-10, 7.075768e-10,
	7.1389966e-10, 7.202424e-10, 7.266073e-10, 7.329966e-10,
	7.394128e-10, 7.4585826e-10, 7.5233547e-10, 7.58847e-10,
	7.653954e-10, 7.719835e-10, 7.7861395e-10, 7.852897e-10,
	7.920138e-10, 7.987892e-10, 8.0561924e-10, 8.125073e-10,
	8.194569e-10, 8.2647167e-10, 8.3355556e-10, 8.407127e-10,
	8.479473e-10, 8.55264e-10, 8.6266755e-10, 8.7016316e-10,
	8.777562e-10, 8.8545243e-10, 8.932582e-10, 9.0117996e-10,
	9.09225e-10, 9.174008e-10, 9.2571584e-10, 9.341788e-10,
	9.427997e-10, 9.515889e-10, 9.605579e-10, 9.697193e-10,
	9.790869e-10, 9.88676e-10, 9.985036e-10, 1.0085882e-09,
	1.0189509e-09, 1.0296151e-09, 1.0406069e-09, 1.0519566e-09,
	1.063698e-09, 1.0758702e-09, 1.0885183e-09, 1.1016947e-09,
	1.1154611e-09, 1.1298902e-09, 1.1450696e-09, 1.1611052e-09,
	1.1781276e-09, 1.1962995e-09, 1.2158287e-09, 1.2369856e-09,
	1.2601323e-09, 1.2857697e-09, 1.3146202e-09, 1.347784e-09,
	1.3870636e-09, 1.4357403e-09, 1.5008659e-09, 1.6030948e-09,
}

var fn = [128]float32{
	1, 0.9635997, 0.9362827, 0.9130436, 0.89228165, 0.87324303,
	0.8555006, 0.8387836, 0.8229072, 0.8077383, 0.793177,
	0.7791461, 0.7655842, 0.7524416, 0.73967725, 0.7272569,
	0.7151515, 0.7033361, 0.69178915, 0.68049186, 0.6694277,
	0.658582, 0.6479418, 0.63749546, 0.6272325, 0.6171434,
	0.6072195, 0.5974532, 0.58783704, 0.5783647, 0.56903,
	0.5598274, 0.5507518, 0.54179835, 0.5329627, 0.52424055,
	0.5156282, 0.50712204, 0.49871865, 0.49041483, 0.48220766,
	0.4740943, 0.46607214, 0.4581387, 0.45029163, 0.44252872,
	0.43484783, 0.427247, 0.41972435, 0.41227803, 0.40490642,
	0.39760786, 0.3903808, 0.3832238, 0.37613547, 0.36911446,
	0.3621595, 0.35526937, 0.34844297, 0.34167916, 0.33497685,
	0.3283351, 0.3217529, 0.3152294, 0.30876362, 0.30235484,
	0.29600215, 0.28970486, 0.2834622, 0.2772735, 0.27113807,
	0.2650553, 0.25902456, 0.2530453, 0.24711695, 0.241239,
	0.23541094, 0.22963232, 0.2239027, 0.21822165, 0.21258877,
	0.20700371, 0.20146611, 0.19597565, 0.19053204, 0.18513499,
	0.17978427, 0.17447963, 0.1692209, 0.16400786, 0.15884037,
	0.15371831, 0.14864157, 0.14361008, 0.13862377, 0.13368265,
	0.12878671, 0.12393598, 0.119130544, 0.11437051, 0.10965602,
	0.104987256, 0.10036444, 0.095787846, 0.0912578, 0.08677467,
	0.0823389, 0.077950984, 0.073611505, 0.06932112, 0.06508058,
	0.06089077, 0.056752663, 0.0526674, 0.048636295, 0.044660863,
	0.040742867, 0.03688439, 0.033087887, 0.029356318,
	0.025693292, 0.022103304, 0.018592102, 0.015167298,
	0.011839478, 0.008624485, 0.005548995, 0.0026696292,
}
