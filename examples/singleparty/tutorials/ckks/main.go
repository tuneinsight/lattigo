package main

import (
	"fmt"
	"math/cmplx"
	"math/rand"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func main() {
	// ============
	// Introduction
	// ============
	//
	// This example showcase the capabilities of encrypted fixed-point approximate arithmetic over the reals/complexes
	// using the implementation of the CKKS scheme as part of the `schemes/ckks` package.
	//
	// The Lattigo library is a library designed around layers.
	// Each layer is a package that provides functionalities for the layers above it.
	//
	// The `ckks` package relies on the `rlwe` package, which themselves relies on the `ring` package:  `ring` -> `rlwe` -> `ckks`.
	//
	// The lowest layer is the `ring` package.
	// The `ring` package provides optimized arithmetic in rings `Z_{Q}[X]/(X^{N}+1)` for `N` a power of two and
	// `QL` the product of `L+1` pairwise NTT friendly primes.
	// It is generic and can be used to implement any scheme based on such rings.
	//
	// The middle layer is the `rlwe` package.
	// This package implements RLWE functionalities that are common to all RLWE based schemes.
	// All objects that are not specific to the CKKS scheme will be imported from the `rlwe` package.
	// Such objects notably `rlwe.Plaintext`, `rlwe.Ciphertext`, `rlwe.SecretKey`, `rlwe.PublicKey` and `rlwe.EvaluationKey`.
	// But also an `rlwe.Evaluator` for all operations that are not scheme specific, such as gadget product and automorphisms,
	// but also more advanced  operations such as the `Trace`.
	//
	// The top layer is the `ckks` package.
	// This package implements the CKKS scheme, and mostly consist in defining the encoding and scheme specific homomorphic operations.
	//

	// =======================================================
	// `rlwe.Ciphertext`, `rlwe.Plaintext` and `rlwe.MetaData`
	// =======================================================
	//
	// Before talking about the capabilities of the `ckks` package, we have to give some information about the `rlwe.Ciphertext` and `rlwe.Plaintext` objects.
	//
	// Both contain the `rlwe.MetaData` struct, which notably holds the following fields:
	//    - `Scale`: the scaling factor. This field is updated dynamically during computations.
	//    - `EncodingDomain`:
	//        - `SlotsDomain`: the usual encoding that provides SIMD operations over the slots.
	//        - `CoefficientDomain`: plain encoding in the RING. Addition behaves as usual, but multiplication will result in negacyclic convolution over the slots.
	//    - `LogSlots`: the log2 of the number of slots. Note that if a ciphertext with n slots is multiplied with a ciphertext of 2n slots, the resulting ciphertext
	//                  will have 2n slots. Because a message `m` of n slots is identical to the message `m|m` of 2n slots.
	//
	// These are all public fields which can be manually edited by advanced users if needed.
	//
	// ======================================================
	// Capabilities of the CKKS Package in the Lattigo Library
	// ======================================================
	//
	// The current capabilities of the `schemes/ckks` package are the following:
	//
	//    - Encoding: encode vectors of type `[]complex128`, `[]float64`, `[]*big.Float` or `[]*bignum.Complex` on `rlwe.Plaintext`
	//
	//    - Addition:
	//        - `rlwe.Ciphertext` + `rlwe.Ciphertext`
	//        - `rlwe.Ciphertext` + `rlwe.Plaintext`
	//        - `rlwe.Ciphertext` + `scalar` of type `complex128`, `float64`, `int`, `int64`, `uint`, `uint64`, `*big.Int`, `*big.Float` or `*bignum.Complex`
	//        - `rlwe.Ciphertext` + `vector` of type `[]complex128`, `[]float64`, `[]*big.Float` or `[]*bignum.Complex`
	//
	//    - Multiplication:
	//        - `rlwe.Ciphertext` * `rlwe.Ciphertext`
	//        - `rlwe.Ciphertext` * `rlwe.Plaintext`
	//        - `rlwe.Ciphertext` * `scalar` of type `complex128`, `float64`, `int`, `int64`, `uint`, `uint64`, `*big.Int`, `*big.Float` or `*bignum.Complex`
	//        - `rlwe.Ciphertext` * `vector` of type `[]complex128`, `[]float64`, `[]*big.Float` or `[]*bignum.Complex`
	//
	//    - Multiplication Fused with Addition (c = c + a*b)
	//        - `rlwe.Ciphertext` + `rlwe.Ciphertext` * `rlwe.Ciphertext`
	//        - `rlwe.Ciphertext` + `rlwe.Ciphertext` * `rlwe.Plaintext`
	//        - `rlwe.Ciphertext` + `rlwe.Ciphertext` * `scalar` of type `complex128`, `float64`, `int`, `int64`, `uint`, `uint64`, `*big.Int`, `*big.Float` or `*bignum.Complex`
	//        - `rlwe.Ciphertext` + `rlwe.Ciphertext` * `vector` of type `[]complex128`, `[]float64`, `[]*big.Float` or `[]*bignum.Complex`
	//
	//    - Rotations & Conjugation
	//
	// On top of the `ckks` package resides the `circuits` module in which specific circuits based on the schemes in `schemes` are implemented.
	// The most salient for the CKKS cryptosystem being:
	//
	//    - Polynomial Evaluation:
	//        - `Single polynomial`: evaluate the same polynomial on all slots of a `rlwe.Ciphertext`
	//        - `Vector polynomial`: evaluate different polynomials on different slots of a `rlwe.Ciphertext`
	//
	//    - Linear Transformations:
	//        - `InnerSum`: aggregate slots inside a `rlwe.Ciphertext`
	//        - `Replicate`: replicate slots inside a `rlwe.Ciphertext`
	//        - `Average`: average the slots inside a `rlwe.Ciphertext`
	//        - `Trace`: evaluate the trace on the slots of `rlwe.Ciphertext`, this
	//        - `LinearTransform`: evaluate a plaintext matrix of type `[][]complex128`, `[][]float64`, `[][]*big.Float` or `[][]*bignum.Complex` on a `rlwe.Ciphertext`
	//
	//    - All methods of the `rlwe.Evaluator`, which are not described here.
	//
	// The `circuits/bootstrapping` package also contains the sub-packages: `bootstrapper` which implements bootstrapping to refresh ciphertexts, enabling arbitrary depth circuits.
	//
	// Note that the package `he/float` also supports a real variant, i.e. plaintext vector of R^{N} (instead of complex vectors C^{N/2}).
	// A homomorphic bridge between the two schemes is also available.
	// This variant can be activated by specifying the `ring.Type` to `ring.ConjugateInvariant` (i.e the ring Z[X + X^{-1}]/(X^{N}+1)) in the `ckks.Parameters` struct.

	// =================================
	// Instantiating the ckks.Parameters
	// =================================
	//
	// We will instantiate a `ckks.Parameters` struct.
	// Unlike other libraries, `Lattigo` doesn't have, yet, a quick constructor.
	// Users must specify all parameters, up to each individual prime size.
	//
	// We will create parameters that are 128-bit secure and allow a depth 7 computation with a scaling factor of 2^{45}.

	var err error
	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(
		ckks.ParametersLiteral{
			LogN:            14,                                    // A ring degree of 2^{14}
			LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45}, // An initial prime of 55 bits and 7 primes of 45 bits
			LogP:            []int{61},                             // The log2 size of the key-switching prime
			LogDefaultScale: 45,                                    // The default log2 of the scaling factor
		}); err != nil {
		panic(err)
	}

	// The ratio between the first prime of size ~2^{55} and the scaling factor 2^{45} is ~2^{10}.
	// This means that these parameter can accommodate for values as large as 2^{9} (signed values).
	// To be able to store larger values, either the scale has to be reduced or the first prime increased.
	// Because the maximum size for the primes of the modulus Q is 60, if we want to store larger values
	// with precision, we will need to reserve the first two primes.

	// We get the encoding precision of the parameters in bits, which is min(53, log2(DefaultScale)).
	// It is always at least 53 (double float precision).
	// This precision is notably the precision used by the encoder to encode/decode values.
	prec := params.EncodingPrecision() // we will need this value later

	// Note that the following fields in the `ckks.ParametersLiteral`are optional, but can be manually specified by advanced users:
	//   - `Xs`: the secret distribution (default uniform ternary)
	//   - `Xe`: the error distribution (default discrete Gaussian with standard deviation of 3.2 and truncated to 19)
	//   - `PowBase`: the log2 of the binary decomposition (default 0, i.e. infinity, i.e. no decomposition)
	//   - `RingType`: the ring to be used, (default Z[X]/(X^{N}+1))
	//
	// We can check the total logQP of the parameters with `params.LogQP()`.
	// For a ring degree 2^{14}, we must ensure that LogQP <= 438 to ensure at least 128 bits of security.

	// ==============
	// Key Generation
	// ==============
	//
	// To generate any key, be it the secret key, the public key or evaluation keys, we first need to instantiate the key generator.
	kgen := rlwe.NewKeyGenerator(params)

	// For now we will generate the following keys:
	//   - SecretKey: the secret from which all other keys are derived
	//   - PublicKey: an encryption of zero, which can be shared and enable anyone to encrypt plaintexts.
	//   - RelinearizationKey: an evaluation key which is used during ciphertext x ciphertext multiplication to ensure ciphertext compactness.
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk) // Note that we can generate any number of public keys associated to the same Secret Key.
	rlk := kgen.GenRelinearizationKeyNew(sk)

	// To store and manage the loading of evaluation keys, we instantiate a struct that complies to the `rlwe.EvaluationKeySetInterface` Interface.
	// The package `rlwe` provides a simple struct that complies to this interface, but a user can design its own struct compliant to the `rlwe.EvaluationKeySetInterface`
	// for example to manage the loading/saving/persistence of the keys in the memory.
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	// ====================
	// Plaintext Generation
	// ====================
	//
	// We use the default number of slots, which is N/2.
	// It is possible to use less slots, however it most situations, there is no reason to do so.
	LogSlots := params.LogMaxSlots()
	Slots := 1 << LogSlots

	// We generate a vector of `[]complex128` with both the real and imaginary part uniformly distributed in [-1, 1]
	/* #nosec G404 -- this is a plaintext vector  */
	r := rand.New(rand.NewSource(0))
	values1 := make([]complex128, Slots)
	for i := 0; i < Slots; i++ {
		values1[i] = complex(2*r.Float64()-1, 2*r.Float64()-1)
	}

	// We allocate a new plaintext, at the maximum level.
	// We can allocate plaintexts at lower levels to optimize memory consumption for operations that we know will happen at a lower level.
	// Plaintexts (and ciphertexts) are by default created with the following metadata:
	//   - `Scale`: `params.DefaultScale()` (which is 2^{45} in this example)
	//   - `EncodingDomain`: `rlwe.SlotsDomain` (this is the default value)
	//   - `LogSlots`: `params.MaxLogSlots` (which is LogN-1=13 in this example)
	// We can check that the plaintext was created at the maximum level with pt1.Level().
	pt1 := ckks.NewPlaintext(params, params.MaxLevel())

	// Then we need to instantiate the encoder, which will enable us to embed our `values` of type `[]complex128` on a `rlwe.Plaintext`.
	// By default the encoder will use the params.DefaultPrecision(), but a user can specify a custom precision as an optional argument,
	// for example `ckks.NewEncoder(params, 256)`.
	ecd := ckks.NewEncoder(params)

	ecd2 := ckks.NewEncoder(ckks.Parameters(params))

	// And we encode our `values` on the plaintext.
	// Note that the encoder will check the metadata of the plaintext and adapt the encoding accordingly.
	// For example, one can modify the `Scale`, `EncodingDomain` or `LogSlots` fields change the way the encoding behaves.
	if err = ecd2.Encode(values1, pt1); err != nil {
		panic(err)
	}

	// =====================
	// Ciphertext Generation
	// =====================
	//
	// To generate ciphertexts we need an encryptor.
	// An encryptor will accept both a secret key or a public key,
	// in this example we will use the public key.
	enc := rlwe.NewEncryptor(params, pk)

	// And we create the ciphertext.
	// Note that the metadata of the plaintext will be copied on the resulting ciphertext.
	ct1, err := enc.EncryptNew(pt1)
	if err != nil {
		panic(err)
	}

	// It is also possible to first allocate the ciphertext the same way it was done
	// for the plaintext with with `ct := ckks.NewCiphertext(params, 1, pt.Level())`,
	// enabling allocation free encryptions (for example if the ciphertext has to be
	// serialized right away).

	// =========
	// Decryptor
	// =========
	//
	// We are able to generate ciphertext from plaintext using the encryptor.
	// To do the converse, generate plaintexts from ciphertexts, we need to instantiate a decryptor.
	// Obviously, the decryptor will only accept the secret key.
	dec := rlwe.NewDecryptor(params, sk)

	// ================
	// Evaluator Basics
	// ================
	//
	// Before anything, we must instantiate the evaluator, and we provide the evaluation key struct.
	eval := ckks.NewEvaluator(params, evk)

	// For the purpose of the example, we will create a second vector of random values.
	values2 := make([]complex128, Slots)
	for i := 0; i < Slots; i++ {
		values2[i] = complex(2*r.Float64()-1, 2*r.Float64()-1)
	}

	pt2 := ckks.NewPlaintext(params, params.MaxLevel())

	// ===========================
	// Managing the Scaling Factor
	// ===========================
	//
	// Before going further and showcasing the capabilities of the evaluator, we must talk
	// about the maintenance of the scaling factor.
	// This is a very central topic, especially for the full-RNS variant of fixed-point
	// approximate homomorphic encryption over the reals/complexes.
	// Messages are encoded on integer polynomials, and thus to keep the precision real
	// coefficients need to be scaled before being discretized to integers.
	// When two messages are multiplied together, the scaling factor of the resulting message
	// is the product of the two initial scaling factors.
	//
	// For example, let D0 * m0 and D1 * m1, be two messages scaled by D0 and D1 respectively.
	// Their multiplication will result in a new messages D0 * D1 * m0 * m1.
	// This means that without any maintenance, the scaling factor will grow exponentially.
	//
	// To control the growth of the scaling factor, we have the rescaling operation.
	// The rescaling operation divides a ciphertext by the prime of its current level and
	// returns a new ciphertext with one less level and scaling factor divided by this prime.
	//
	// The main  difficulty arises from the primes used for the rescaling, since they do not
	// divide the scaling factor.
	//
	// Throughout this example we will show ways to properly manage this scaling factor to both
	// keep it as close as possible to the default scaling factor (in this example 2^{45}) and
	// minimizing the error.
	// In fact we will show that it is usually possible to keep the scaling factor always at 2^{45},
	// even though the primes are not powers of two.

	fmt.Printf("========\n")
	fmt.Printf("ADDITION\n")
	fmt.Printf("========\n")
	fmt.Printf("\n")
	// Additions are often seen as a trivial operation.
	// However in the case of the full-RNS implementation we have to be careful.
	// Indeed, we must ensure that when adding two ciphertexts, those ciphertexts have the same exact scale,
	// else an error proportional to the difference of the scale will be introduced.
	//
	// The evaluator will try to compensate if the ciphertexts do not have the same scale,
	// but only up to an integer multiplication (which is "free").
	// This means that if one scale is an integer multiple of the other (e.g. 2^{45} * q0 and 2^{45}),
	// then the evaluator will take that into account and properly operate the addition.
	//
	// However, if one of the scales is a fraction of the other (e.g. 2^{45} * q0 and 2^{45} * q1),
	// the evaluator isn't able to reconciliate the scales and will treat the ciphertext with the
	// smallest scale as being at the scale of the largest one.
	//
	// This will introduce an approximation error proportional to q^{45} * q0 / 2^{45} * q1 = q0/q1 in the addition.
	//
	// Thus, when users are manually calling the addition between ciphertexts and/or plaintexts,
	// they must ensure that both operands have scales that are an integer multiple of the other.

	// ciphertext + ciphertext
	if err = ecd.Encode(values2, pt2); err != nil {
		panic(err)
	}

	ct2, err := enc.EncryptNew(pt2)
	if err != nil {
		panic(err)
	}

	want := make([]complex128, Slots)
	for i := 0; i < Slots; i++ {
		want[i] = values1[i] + values2[i]
	}

	// A small comment about the precision stats.
	// Theses stats show the -log2 of the matching bits on the right side of the decimal point.
	// Because values are not normalized, large values will show as having a low precision, even if left side of of the decimal point (integer part) is correct.
	// Eventually this will be fixed, by normalizing with the maximum value decrypted.
	ct3, err := eval.AddNew(ct1, ct2)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Addition - ct + ct%s", ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	// ciphertext + plaintext
	ct3, err = eval.AddNew(ct1, pt2)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Addition - ct + pt%s", ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	// ciphertext + vector
	// Note that the evaluator will encode this vector at the scale of the input ciphertext to ensure a noiseless addition.
	ct3, err = eval.AddNew(ct1, values2)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Addition - ct + vector%s", ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	// ciphertext + scalar
	scalar := 3.141592653589793 + 1.4142135623730951i
	for i := 0; i < Slots; i++ {
		want[i] = values1[i] + scalar
	}

	// Similarly, if we give a scalar, it will be scaled by the scale of the input ciphertext to ensure a noiseless addition.
	ct3, err = eval.AddNew(ct1, scalar)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Addition - ct + scalar%s", ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	fmt.Printf("==============\n")
	fmt.Printf("MULTIPLICATION\n")
	fmt.Printf("==============\n")
	fmt.Printf("\n")

	for i := 0; i < Slots; i++ {
		want[i] = values1[i] * values2[i]
	}

	// We could simply call the multiplication on ct1 and ct2, however since a rescaling is needed afterward,
	// we also want to properly control the scale of the result.
	// Our goal is to keep the scale to the default one, i.e. 2^{45} in this example.
	// However, the rescaling operation divides by one (or multiple) primes qi,
	// with the shape 2^{s} +/- k*2N + 1, which are obviously not powers of two.
	// The best way to achieve this goal is to ensure that the scale before the rescaling is 2^{45} * prime_to_rescale.
	// This way the division is exact and we fall back on the default scaling factor.
	//
	// Given a ciphertext of scale 2^{45}, the easiest way to achieve this result is to scale ct2
	// by the prime that will be used by the rescaling, which params.Q()[min(ct1.Level(), ct2.Level())].
	//
	// So, for this example, we will show how to create a new ciphertext at the correct scale.
	//
	// To do so, we manually specify the scaling factor of the plaintext:
	pt2.Scale = rlwe.NewScale(params.Q()[ct1.Level()])

	// Then we encode the values (recall that the encoding is done according to the metadata of the plaintext)
	if err = ecd.Encode(values2, pt2); err != nil {
		panic(err)
	}

	// and we encrypt (recall that the metadata of the plaintext are copied on the created ciphertext)
	if err := enc.Encrypt(pt2, ct2); err != nil {
		panic(err)
	}

	res, err := eval.MulRelinNew(ct1, ct2)
	if err != nil {
		panic(err)
	}

	// The scaling factor of res should be equal to ct1.Scale * ct2.Scale
	ctScale := &res.Scale.Value // We need to access the pointer to have it display correctly in the command line
	fmt.Printf("Scale before rescaling: %f\n", ctScale)

	// To control the growth of the scaling factor, we call the rescaling operation.
	// Such rescaling operation should be called at the latest before the next multiplication.
	// Each rescaling operation consumes a level, reducing the homomorphic capacity of the ciphertext.
	// If a ciphertext reaches the level 0, it can no longer be rescaled and any further multiplication
	// risks inducing a plaintext overflow.
	if err = eval.Rescale(res, res); err != nil {
		panic(err)
	}

	Scale := params.DefaultScale().Value

	// And we check that we are back on our feet with a scale of 2^{45} but with one less level
	fmt.Printf("Scale after rescaling: %f == %f: %t and %d == %d+1: %t\n", ctScale, &Scale, ctScale.Cmp(&Scale) == 0, ct1.Level(), res.Level(), ct1.Level() == res.Level()+1)
	fmt.Printf("\n")

	// For the sake of conciseness, we will not rescale the output for the other multiplication example.
	// But this maintenance operation should usually be called (either before of after the multiplication depending on the choice of noise management)
	// to control the magnitude of the plaintext scale.
	fmt.Printf("Multiplication - ct * ct%s", ckks.GetPrecisionStats(params, ecd, dec, want, res, 0, false).String())

	// ciphertext + plaintext
	ct3, err = eval.MulRelinNew(ct1, pt2)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Multiplication - ct * pt%s", ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	// ciphertext + vector
	// Note that when giving non-encoded vectors, the evaluator will internally encode this vector with the appropriate scale that ensure that
	// the following rescaling operation will make the resulting ciphertext fall back on it's previous scale.
	ct3, err = eval.MulRelinNew(ct1, values2)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Multiplication - ct * vector%s", ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	// ciphertext + scalar (scalar = pi + sqrt(2) * i)
	for i := 0; i < Slots; i++ {
		want[i] = values1[i] * scalar
	}

	// Similarly, when giving a scalar, the scalar is encoded with the appropriate scale to get back to the original ciphertext scale after the rescaling.
	// Additionally, the multiplication with a Gaussian integer does not increase the scale of the ciphertext, thus does not require rescaling and does not consume a level.
	// For example, multiplication/division by the imaginary unit `i` is free in term of level consumption and can be used without moderation.
	ct3, err = eval.MulRelinNew(ct1, scalar)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Multiplication - ct * scalar%s", ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	fmt.Printf("======================\n")
	fmt.Printf("ROTATION & CONJUGATION\n")
	fmt.Printf("======================\n")
	fmt.Printf("\n")

	// Before being able to do any rotations, we need to generate the corresponding Galois keys.
	// A Galois key is a special type of `rlwe.EvaluationKey` that enables automorphisms X^{i} -> X^{i*k mod 2N} mod X^{N} + 1 on ciphertext
	// Some of these automorphisms act like cyclic rotations on plaintext encoded in the `SlotsDomain`.
	//
	// Galois keys can be large depending on the parameters, and one Galois key is needed per automorphism.
	// Therefore it is important to design circuits that minimize the numbers of these keys.
	//
	// In this example we will rotate a ciphertext by 5 positions to the left, as well as get the complex conjugate.
	// This corresponds to the following values for k which we call "galois elements":
	rot := 5
	galEls := []uint64{
		// The galois element for the cyclic rotations by 5 positions to the left.
		params.GaloisElement(rot),
		// The galois element for the complex conjugatation.
		params.GaloisElementForComplexConjugation(),
	}

	// We then generate the `rlwe.GaloisKey`s element that corresponds to these galois elements.
	// And we update the evaluator's `rlwe.EvaluationKeySet` with the new keys.
	eval = eval.WithKey(rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(galEls, sk)...))

	// Rotation by 5 positions to the left
	for i := 0; i < Slots; i++ {
		want[i] = values1[(i+5)%Slots]
	}

	ct3, err = eval.RotateNew(ct1, rot)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Rotation by k=%d %s", rot, ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	// Conjugation
	for i := 0; i < Slots; i++ {
		want[i] = complex(real(values1[i]), -imag(values1[i]))
	}

	ct3, err = eval.ConjugateNew(ct1)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Conjugation %s", ckks.GetPrecisionStats(params, ecd, dec, want, ct3, 0, false).String())

	// Note that rotations and conjugation only add a fixed additive noise independent of the ciphertext noise.
	// If the parameters are set correctly, this noise can be rounding error (thus negligible).
	// It is recommended apply the rescaling operation after such operations rather than before.
	// This way, the noise is added in the lower bits of the ciphertext and gets erased by the rescaling.

	fmt.Printf("=====================\n")
	fmt.Printf("POLYNOMIAL EVALUATION\n")
	fmt.Printf("=====================\n")
	fmt.Printf("\n")

	// The evaluator can evaluate polynomials in standard and Chebyshev basis.
	// The evaluation is optimal in depth consumption and ensures that all additions are noiseless.
	// The package `utils/bignum` also provide a way to approximate smooth functions with a Chebyshev interpolation.
	// Eventually, we will also add the multi-interval minimax approximation.
	//
	// Let define a function, for example, the SiLU.
	// The signature needed is `func(x *bignum.Complex) (y *bignum.Complex)` so we must accommodate for it first:

	// Yes SiLU over the complex!
	SiLU := func(x complex128) (y complex128) {
		return x / (cmplx.Exp(-x) + 1)
	}

	// We must also give an interval [a, b], for example [-8, 8], in which we approximate SiLU, as well as the degree of approximation.
	// With 7 levels, we can evaluate a polynomial of degree up to 127.
	// However, since we will be in the Chebyshev basis, we must also take into consideration the change of basis
	// y = (2*x - a - b)/(b-a), which will usually consume a level.
	// Often it is however possible to include this linear transformation in previous step of a circuit, to save a level.
	// Since we do not have any previous operation in this example, we will have to operate the change of basis, thus
	// the maximum polynomial degree for depth 6 is 63.

	interval := bignum.Interval{
		Nodes: 63,
		A:     *bignum.NewFloat(-8, prec),
		B:     *bignum.NewFloat(8, prec),
	}

	// We generate the `bignum.Polynomial` which stores the degree 63 Chevyshev approximation of the SiLU function in the interval [-8, 8]
	poly := bignum.ChebyshevApproximation(SiLU, interval)

	// The struct `bignum.Polynomial` comes with an handy evaluation method
	tmp := bignum.NewComplex().SetPrec(prec)
	for i := 0; i < Slots; i++ {
		want[i] = poly.Evaluate(tmp.SetComplex128(values1[i])).Complex128()
	}

	// First, we must operate the change of basis for the Chebyshev evaluation y = (2*x-a-b)/(b-a) = scalarmul * x + scalaradd
	scalarmul, scalaradd := poly.ChangeOfBasis()

	res, err = eval.MulNew(ct1, scalarmul)
	if err != nil {
		panic(err)
	}

	if err = eval.Add(res, scalaradd, res); err != nil {
		panic(err)
	}

	if err = eval.Rescale(res, res); err != nil {
		panic(err)
	}

	polyEval := polynomial.NewEvaluator(params, eval)

	// And we evaluate this polynomial on the ciphertext
	// The last argument, `params.DefaultScale()` is the scale that we want the ciphertext
	// to have after the evaluation, which is usually the default scale, 2^{45} in this example.
	// Other values can be specified, but they should be close to the default scale, else the
	// depth consumption will not be optimal.
	if res, err = polyEval.Evaluate(res, poly, params.DefaultScale()); err != nil {
		panic(err)
	}

	fmt.Printf("Polynomial Evaluation %s", ckks.GetPrecisionStats(params, ecd, dec, want, res, 0, false).String())

	// =============================
	// Vector Polynomials Evaluation
	// =============================
	//

	fmt.Printf("======================\n")
	fmt.Printf("LINEAR TRANSFORMATIONS\n")
	fmt.Printf("======================\n")
	fmt.Printf("\n")

	// The `circuits/lintrans` package provides a multiple handy linear transformations.
	// We will start with the inner sum.
	// Thus method allows to aggregate `n` sub-vectors of size `batch`.
	// For example given a vector [x0, x1, x2, x3, x4, x5, x6, x7], batch = 2 and n = 3
	// it will return the vector [x0+x2+x4, x1+x3+x5, x2+x4+x6, x3+x5+x7, x4+x6+x0, x5+x7+x1, x6+x0+x2, x7+x1+x3]
	// Observe that the inner sum wraps around the vector, this behavior must be taken into account.

	batch := 37
	n := 127

	// The innersum operations is carried out with log2(n) + HW(n) automorphisms and we need to
	// generate the corresponding Galois keys and provide them to the `Evaluator`.
	eval = eval.WithKey(rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(params.GaloisElementsForInnerSum(batch, n), sk)...))

	// Plaintext circuit
	copy(want, values1)
	for i := 1; i < n; i++ {
		for j, vi := range utils.RotateSlice(values1, i*batch) {
			want[j] += vi
		}
	}

	if err := eval.InnerSum(ct1, batch, n, res); err != nil {
		panic(err)
	}

	// Note that this method can obviously be used to average values.
	// For a good noise management, it is recommended to first multiply the values by 1/n, then
	// apply the innersum and then only apply the rescaling.
	fmt.Printf("Innersum %s", ckks.GetPrecisionStats(params, ecd, dec, want, res, 0, false).String())

	// The replicate operation is exactly the same as the innersum operation, but in reverse
	eval = eval.WithKey(rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(params.GaloisElementsForReplicate(batch, n), sk)...))

	// Plaintext circuit
	copy(want, values1)
	for i := 1; i < n; i++ {
		for j, vi := range utils.RotateSlice(values1, -i*batch) { //Note the minus sign
			want[j] += vi
		}
	}

	if err := eval.Replicate(ct1, batch, n, res); err != nil {
		panic(err)
	}

	fmt.Printf("Replicate %s", ckks.GetPrecisionStats(params, ecd, dec, want, res, 0, false).String())

	// And we arrive to the linear transformation.
	// This method enables to evaluate arbitrary Slots x Slots matrices on a ciphertext.
	// What matters is not the size of the matrix, but the number of non-zero diagonals, as
	// the complexity of this operation is 2sqrt(#non-zero-diags).
	//
	// First lets explain what we mean by non-zero diagonal.
	// As an example, lets take the following 4x4 matrix:
	//   0 1 2 3 (diagonal index)
	// | 1 2 3 0 |
	// | 0 1 2 3 |
	// | 3 0 1 2 |
	// | 2 3 0 1 |
	//
	// This matrix has 3 non zero diagonals at indexes [0, 1, 2]:
	//   - 0: [1, 1, 1, 1]
	//   - 1: [2, 2, 2, 2]
	//   - 2: [3, 3, 3, 3]
	//

	nonZeroDiagonals := []int{-15, -4, -1, 0, 1, 2, 3, 4, 15}

	// We allocate the non-zero diagonals and populate them
	diagonals := make(lintrans.Diagonals[complex128])

	for _, i := range nonZeroDiagonals {
		tmp := make([]complex128, Slots)

		for j := range tmp {
			tmp[j] = complex(2*r.Float64()-1, 2*r.Float64()-1)
		}

		diagonals[i] = tmp
	}

	// We create the linear transformation of type complex128 (float64, *big.Float and *bignum.Complex are also possible)
	// Here we use the default structs of the rlwe package, which is compliant to the rlwe.LinearTransformationParameters interface
	// But a user is free to use any struct compliant to this interface.
	// See the definition of the interface for more information about the parameters.
	ltparams := lintrans.Parameters{
		DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
		LevelQ:                    ct1.Level(),
		LevelP:                    params.MaxLevelP(),
		Scale:                     rlwe.NewScale(params.Q()[ct1.Level()]),
		LogDimensions:             ct1.LogDimensions,
		LogBabyStepGiantStepRatio: 1,
	}

	// We allocated the rlwe.LinearTransformation.
	// The allocation takes into account the parameters of the linear transformation.
	lt := lintrans.NewTransformation(params, ltparams)

	// We encode our linear transformation on the allocated rlwe.LinearTransformation.
	// Not that trying to encode a linear transformation with different non-zero diagonals,
	// plaintext dimensions or baby-step giant-step ratio than the one used to allocate the
	// rlwe.LinearTransformation will return an error.
	if err := lintrans.Encode(ecd, diagonals, lt); err != nil {
		panic(err)
	}

	// Then we generate the corresponding Galois keys.
	// The list of Galois elements can also be obtained with `lt.GaloisElements`
	// but this requires to have it pre-allocated, which is not always desirable.
	galEls = lintrans.GaloisElements(params, ltparams)

	ltEval := lintrans.NewEvaluator(eval.WithKey(rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(galEls, sk)...)))

	// And we valuate the linear transform
	if err := ltEval.Evaluate(ct1, lt, res); err != nil {
		panic(err)
	}

	// Result is not returned rescaled
	if err = eval.Rescale(res, res); err != nil {
		panic(err)
	}

	// We evaluate the same circuit in plaintext
	want = diagonals.Evaluate(values1, newVec, add, muladd)

	fmt.Printf("vector x matrix %s", ckks.GetPrecisionStats(params, ecd, dec, want, res, 0, false).String())

	// =============================
	// Homomorphic Encoding/Decoding
	// =============================
	//

	// ============
	// Bootstrapping
	// ============
	//

	// ==========
	// CONCURRENCY
	// ==========
	//
	// Lattigo does not implement low level concurrency yet.
	// Currently concurrency must be done at the circuit level.
	//
	// By design, structs outside of the parameters are not thread safe.
	// For example, one cannot use an encoder to encode concurrently on different plaintexts.
	// However, all structs (for which it makes sens) have the method `ShallowCopy`, which creates
	// a copy of the original struct with new internal buffers, that is safe to use concurrently.

}

func newVec(size int) []complex128 {
	return make([]complex128, size)
}

func add(a, b, c []complex128) {
	for i := range a {
		c[i] = a[i] + b[i]
	}
}

func muladd(a, b, c []complex128) {
	for i := range a {
		c[i] = c[i] + a[i]*b[i]
	}
}
