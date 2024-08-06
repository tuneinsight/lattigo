package ckks

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

var (
	// ExampleParameters128BitLogN14LogQP438 is an example parameters set with logN=14, logQP=435
	// offering 128-bit of security.
	ExampleParameters128BitLogN14LogQP438 = ParametersLiteral{
		// LogN is the log2 of the ring of the ring degree, i.e. log2(2^{14})
		LogN: 14,

		// Q is the ciphertext modulus, which is the product of pair-wise co-prime primes congruent to 1 modulo 2N.
		// Each prime after the first one adds one `level` (i.e. one to the total depth that the parameters can support)
		// and should be usually as close as possible to 2^{LogDefaultScale}.
		// The first prime must be large enough to store the result of the computation. In this example parameters, the
		// first prime is close to 2^{55} and the expected final scaling factor is 2^{LogDefaultScale}, thus the gap
		// between Q[0] and 2^{LogDefaultScale} is 2^{10}, leaving room for a plaintext message whose magnitude can
		// be up to 2^{9} (one bit is reserved for the sign).
		Q: []uint64{
			0x80000000080001, // 55
			0x2000000a0001,   // 45
			0x2000000e0001,   // 45
			0x2000001d0001,   // 45
			0x1fffffcf0001,   // 45
			0x1fffffc20001,   // 45
			0x200000440001,   // 45
		},

		// LogQ allows to specify the primes of Q by their bit size instead.
		//LogQ: []int{}

		// LogDefaultScale is the log2 of the initial scaling factor (i.e the scaling factor that ciphertext and plaintext
		// have by default when allocated).
		// This value is usually the same
		LogDefaultScale: 45,

		// Optional parameters:

		// P is an optional auxiliary modulus added on top of Q for the evaluation keys. This modulus does not contribute
		// to the homomorphic capacity, but has to be taken into account, along with Q, when estimating the the security
		// of a parameter set.
		//
		// The RNS decomposition during the key-switching operation uses by default the primes of Q as decomposition basis,
		// i.e. [q0, q1, q2, ...] and this introduces a noise proportional to the sum of the primes composing Q. To mitigate
		// this noise, we can add an auxiliary prime P that satisfies |P| >= max(|qi|) and that will divide final noise by
		// that same amount, leaving a residual rounding noise (i.e. negligible).
		//
		// Using a single P is practical only up to a certain point because the complexity of the key-switching (and size of the
		// evaluation keys) is proportional to the number of primes in Q times the number of elements in the decomposition basis.
		// Thus by default it is quadratic in the number of primes in Q. To reduce the size of the evaluation keys and the complexity
		// of the key-switching operation, the user can give more than one prime for P. The number of primes in P will determine
		// the size of each element of the RNS decomposition basis. For example, if 3 primes for P are given, then the decomposition
		// basis will be triplets of primes of Q, i.e. [(q0 * q1 * q2), (q3 * q4 * q5), ...]. The number of elements in the decomposition
		// basis is therefor reduced by a factor of 3, and so are the size of the keys and the complexity of the key-switching.
		// As a rule of thumb, allocating sqrt(#qi) primes to P is a good starting point.
		// The drawback of adding more primes to P is these primes contribute to the total modulus used to estimate the security
		// but not to the total homomorphic capacity.
		// For additional information about this hybrid key-switching, see Section 3 of Better Bootstrapping for Approximate Homomorphic
		// Encryption (https://eprint.iacr.org/2019/688.pdf).
		//
		// It is also possible to not allocate any prime to P but then the default RNS decomposition (modulo each prime of Q) will add
		// a substantial error (proportional to the sum of the primes of Q). However, it is still possible to to mitigate this noise by
		// adding an extra power of two decomposition on the evaluation keys. However such decomposition is independent from the parameters
		// and therefor not specified here. See the description of `rlwe.EvaluationKey` and `rlwe.EvaluationKeyParametersLiteral`
		// for additional information.
		P: []uint64{
			0x80000000130001, // 55
			0x7fffffffe90001, // 55
		},

		// LogP allows to specify the primes of P by their bit size instead.
		//LogP: []int{}

		// RingType denotes the type of ring in which we work. By default (ring.Standard) it is Z[X]/(X^{2^{LogN}}+1), and this
		// will instantiates the regular CKKS scheme, i.e. with N/2 slots over the complex domain. However another ring type
		// ring.ConjugateInvariant can be given, which is defined as Z[X + X^-1]/(X^{2^{LogN}}+1). This special ring will instantiate
		// a variant of the CKKS scheme with instead N slots over the reals. See Approximate Homomorphic Encryption over the Conjugate-
		// invariant Ring  for additional details (https://eprint.iacr.org/2018/952).
		// Note that a homomorphic bridge between the two rings exists (see ckks/bridge.go).
		RingType: ring.Standard,

		// Xs is the secret distribution. The default value is a ternary secret with density 2/3 (i.e. each coefficient as an equal chance
		// of being -1, 0 or 1). Other distributions are supported, such as ternary secret with fixed Hamming weight or an error distribution.
		// See lattigo/ring/sampler.go for the available distributions and their parameters.
		Xs: rlwe.DefaultXs,

		// Xe is the error distribution. The default value is a discrete Gaussian with standard deviation 3.2 and bounded by 19.
		// Other distributions are supported, see lattigo/ring/sampler.go for the available distributions and their parameters.
		Xe: rlwe.DefaultXe,
	}
)
