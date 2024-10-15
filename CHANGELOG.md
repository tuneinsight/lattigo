# Changelog
All notable changes to this library are documented in this file.

## [6.1.0] - 04.10.2024
- Update of `PrecisionStats` in `ckks/precision.go`:  
  - The precision is now computed as the min/max/average/... of the log of the error (instead of the log of the min/max/average/... of the error).
  - fields renamed (`MinPrecision` -> `MINLog2Prec`, `MaxPrecision` -> `MAXLog2Prec`, ...)
  - `rlwe.Scale` has a `.Log2()` method
- Update of `mod1.Parameters` fields (made public, some removed)
- Improvement of the relinearization key-generation protocol (reduce the degree of the shares)
- Serialization of bootstrapping keys
- Lower noise incurred by `ModUp`
- Evaluation keys can be compressed (public element `a` can be generated from a seed)
- More doc formatting 
- Fix various bugs: 
  - `ShallowCopy` of the CKKS bootstrapping evaluator and BFV evaluator not deep enough.
  - PSI example failing 
  - Incorrect reset of pointer in uniform sampler
  - Error when doing inverse NTT with small degree 
  - Mod1Evaluator changes the input ciphertext

## [6.0.0] - 06.08.2024
- Deprecated Go versions `1.18`, `1.19` and `1.20`. The minimum version is now `1.21`, due to the use of the slices of std library package.
  - Removal of all slice utility functions in `utils/slices.go` that are provided in the standard library.
- Golang test cache is deleted before invoking the unit test suite via `make test`.
- Deletion of the `he` package and its abstraction layers for the BGV and CKKS, refocus of the library onto the scheme level, i.e., the `schemes` package.
- Extraction of the homomorphic circuits from the `he` into a new `circuits` package.
  - Simplification of the API for several circuits:
	- Removal of the circuit-specific evaluator interfaces, e.g., `EvaluatorForLinearTransformation`. These interfaces are replaced with a scheme-agnostic evaluator in `schemes/schemes.go` due to the refocus of the Lattigo towards the individual cryptosystems.
	- The individual homomorphic circuits are organized by schemes, in other words in the packages `circuits/bgv`, `circuits/ckks` and `circuits/common` where the latter bundles scheme-generic functionalities of circuits common to all deriving schemes.
- Absorb the `bfv` package into the `bgv` package.
- Rename `mhe` into `multiparty`.
- Slot-wise permutations as part of the `LinearTransformation` circuits:
  - Permutation are handled through `lintrans.Permutation` and `lintrans.PermutationMapping` that induce `lintrans.Diagonals` from which a regular linear transformation can be bootstrapped.
- New ring packing API:
  - New packing evaluator: `rlwe.RingPackingEvaluator`:
	- `NewRingPackingEvaluator(evk *RingPackingEvaluationKey)`
	- `Extract(ct *rlwe.Ciphertext, idx []int, naive bool) (cts map[int]*rlwe.Ciphertext, err error)`
	- `Repack(cts map[int]*rlwe.Ciphertext, naive bool) (ct *rlwe.Ciphertext, err error)`
	- `Split(ctN, ctEvenNHalf, ctOddNHalf *rlwe.Ciphertext) (err error)`
	- `Merge(ctEvenNHalf, ctOddNHalf, ctN *rlwe.Ciphertext) (err error)`
	- `ShallowCopy() *RingPackingEvaluator`
  - New packing evaluation key:
	- `rlwe.RingPackingEvaluationKey`
- Streamlined unit test context generation to reduce boilerplate code duplication:
  - `schemes/*/test_utils.go` to be reused in packages `schemes` and packages that depend on `schemes`.
- Introduction of the `lintrans.Parameters.LevelQ` and `lintrans.Parameters.LevelP` fields and the removal of the `lintrans.Parameters.Level` field.
- Optimizations:
  - Reduce the number of masking operations in `rlwe.Evaluator.gadgetProductSinglePAndBitDecompLazy`.
  - Reduce degree of relinearization key shares in `multiparty/keygen_relin.go`.
- Add `gosec` exception directive to weak randomness in unit tests.
- Fix various linter warnings.
- Various docstring formatting fixes and the addition of `godoc` links through the `[]` operator.
- Updated `README.md` with new package hierarchy figures `lattigo-hierachy.svg` and new issue policy.

## [5.0.0] - 15.11.2023
- Deprecated Go versions `1.14`, `1.15`, `1.16`, and `1.17`. The minimum version is now `1.18`, due to the required use of generics.
- Golang Security Checker pass.
- Dereferenced most inputs and pointers methods whenever possible. Pointers methods/inputs are now mostly used when the struct implementing the method and/or the input is intended to be modified.
- Improved serialization interface:
    - Low-entropy structs (such as parameters or rings) have been updated to use more compatible `json.Marshal` as underlying marshaller.
    - High-entropy structs, such as structs storing keys or encrypted values now all satisfy the following interface:
        - `WriteTo(io.Writer) (int64, error)`: writes the object to a standard `io.Writer` interface. The method is optimized and most efficient when writing on writers that expose their own internal buffer (see the `buffer.Writer` interface).
        - `ReadFrom(io.Reader) (int64, error)`: reads an object from a standard `io.Reader` interface. The method is optimized and most efficient when reading from readers that expose their own internal buffers (see the `buffer.Writer` interface).
        - `MarshalBinary() ([]byte, error)`: the previously available, standard `encoding.BinaryMarshaler` interface.
        - `UnmarshalBinary([]byte) (error)`: the previously available, standard `encoding.BinaryUnmarshaler` interface.
        - `BinarySize() int`: size in bytes when written to an `io.Writer` or when marshalled.
    - Streamlined and simplified all tests related to serialization. They can now be implemented with a single line of code with `RequireSerializerCorrect` that checks the correctness of the above interface as well as equality between bites written using `WriteTo` and bytes generated using `MarshalBinary`.
- Improved consistency across method names and across packages/schemes:
    - All sub-strings `NoMod`, `NoModDown` and `Constant` in method names have been replaced by the sub-string `Lazy`. For example `AddNoMod` and `MulCoeffsMontgomeryConstant` become `AddLazy` and `MulCoeffsMontgomeryLazy` respectively.
    - All sub-strings `And` in methods names have been replaced by the sub-string `Then`. For example `MulAndAdd` becomes `MulThenAdd`.
    - All sub-strings `Inv` have been replaced by `I` for consistency. For example `InvNTT` becomes `INTT`.
    - All sub-strings `Params` and equivalent, referring to pre-computed constants, have been replaced by `Constant`. For example `ModUpParams` becomes `ModUpConstants`.
-  New top-level packages that provide a more convenient and streamlined user-interface to HE:
    - `he`: Package `he` defines common high-level interfaces and implements common high-level operations in a scheme-agnostic way.
        - The common operations in Linear Transformations
        - The common operations in Polynomial Evaluation
    - `he/hefloat`: Package `hefloat` implements fixed-point approximate encrypted arithmetic over real/complex numbers.
      This package provides all the functionalities of the `schemes/ckks` package, as well as additional more advanced circuits, such as:
        - Linear Transformations
        - Homomorphic encoding/decoding
        - Polynomial Evaluation
        - Composite Minimax Polynomial Evaluation
        - Homomorphic modular reduction (x mod 1)
        - GoldschmidtDivision (x in [0, 2])
        - Full domain division (x in [-max, -min] U [min, max])
        - Sign and Step piece-wise functions (x in [-1, 1] and [0, 1] respectively)
        - Min/Max between values in [-0.5, 0.5]
    - `he/hefloat/bootstrapper`: Package `bootstrapper` implements bootstrapping for fixed-point approximate homomorphic encryption over the real/complex numbers.
    It improves on the original implementation with the following features:
        - Bootstrapping batches of ciphertexts of smaller dimension and/or with sparse packing with automatic ring-degree switching and $0$-depth packing/unpacking.
        - Bootstrapping for the Conjugate Invariant CKKS with optimal throughput.
        - Decorrelation between the bootstrapping parameters and residual parameters: the user doesn't need to manage two sets of parameters anymore and the user 
          only needs to provide the residual parameters (what should remain after the evaluation of the bootstrapping circuit)
        - Out-of-the-box usability with default parameterization independent of the residual parameters.
        - In-depth parameterization for advanced users with 16 tunable parameters.
        - Improved implementation of META-BTS, providing arbitrary precision bootstrapping from only one additional small prime.
    - `he/heint`: Package `heint` implements encrypted modular arithmetic over the integers.
        - Linear Transformations
        - Polynomial Evaluation 
    - `he/hebin`: Package`hebin` implements blind rotations evaluation for R-LWE schemes.
- Moved the default parameters of all schemes to the `examples` package, where they are now referred to as **example** parameter sets to better convey the idea that they should not be used as such in real applications.
- BFV: 
    - The code of the package `bfv` has been replaced by a wrapper of the package `bgv` and moved to the package `schemes/bfv`.
- BGV:
    - The code the `bgv` package has been moved to the package `schemes/bfv`
    - The package `bgv` has been rewritten to implement a unification of the textbook BFV and BGV schemes under a single scheme. This unification offers all the functionalities of the BFV and BGV schemes under a single scheme.
    - Changes to the `Encoder`:
        - `NewEncoder` now returns an `*Encoder` instead of an interface.
        - Updated and uniformized the `Encoder` API. It now satisfies the generic `he.Encoder` interface.
        - The encoding will be performed according to the plaintext `MetaData`.
    - Changes to the `Evaluator`:
        - `NewEvaluator` now returns an `*Evaluator` instead of an interface.
        - Updated and uniformized the `Evaluator` API. It now satisfies the generic `he.Evaluator` interface.
    - Changes to the `Parameters`:
        - Enabled plaintext moduli with a smaller 2N-th root of unity than the ring degree.
        - Replaced the default parameters by a single example parameter.
        - Added a test parameter set with small plaintext modulus.
- CKKS:
    - The code of the `ckks` package has been moved to the package `schemes/ckks`.
    - Changes to the `Encoder`:
        - Enabled the encoding of plaintexts of any sparsity (previously hard-capped at a minimum of 8 slots).
        - Unified `encoderComplex128` and `encoderBigComplex`.
        - Updated and uniformized the `Encoder`API. It now satisfies the generic `he.Encoder` interface.
        - The encoding will be performed according to the plaintext `MetaData`.
    - Changes to the `Evaluator`: 
        - `NewEvaluator` now returns an `*Evaluator` instead of an interface.
        - Updated and uniformized the `Evaluator` API. It now satisfies the generic `he.Evaluator` interface.
        - Improved and generalized the internal implementation of the `Evaluator` to enable arbitrary precision encrypted arithmetic.
    - Changes to the `Parameters`:
        - Replaced the default parameters by a single example parameter.
        - Renamed the field `LogScale` of the `ParametersLiteralStruct` to `LogPlaintextScale`.
    - Changes to the tests:
        - Tests do not use the default parameters anymore but specific and optimized test parameters.
        - Added two test parameters `TESTPREC45` for 45-bit precision and `TESTPREC90` for 90-bit precision.
    - Others:
        - Updated the Chebyshev interpolation with arbitrary precision arithmetic and moved the code to `utils/bignum/approximation`.
- RLWE:
    - The package `rlwe` has been moved to `core/rlwe`.
    - The package `ringqp` has been moved to `ring/ringqp`.
    - Changes to the `Parameters`:
        - It is now possible to specify both the secret and error distributions via the `Xs` and `Xe` fields of the `ParameterLiteral` struct.
        - Removed the concept of rotation, everything is now defined in terms of Galois elements.
        - Renamed methods to better reflect their purpose and to generalize them.
        - Added methods related to plaintext parameters and noise.
        - Removed the field `Pow2Base` which is now a parameter of the struct `EvaluationKey`.
    - Changes to the `Encryptor`:
        - `EncryptorPublicKey` and `EncryptorSecretKey` are now public.
        - Encryptors instantiated with a `rlwe.PublicKey` can now encrypt over `rlwe.ElementInterface[ringqp.Poly]` (i.e. generating of `rlwe.GadgetCiphertext` encryptions of zero with `rlwe.PublicKey`).
    - Changes to the `Decryptor`:
        - `NewDecryptor` returns a `*Decryptor` instead of an interface.
    - Changes to the `Evaluator`:
        - Updated all methods of the `Evaluator` to work with operands in and out of the NTT domain.
        - Renamed `SwitchKeys` to `ApplyEvaluationKey`.
        - Renamed `Evaluator.Merge` to `Evaluator.Pack` and generalized `Evaluator.Pack` to be able to take into account the packing `X^{N/n}` of the ciphertext.
        - `Evaluator.Pack` is not recursive anymore and gives the option to zero (or not) slots which are not multiples of `X^{N/n}`.
        - Added the methods `CheckAndGetGaloisKey` and `CheckAndGetRelinearizationKey` to safely check and get the corresponding `EvaluationKeys`.
        - Added the method `InnerFunction`, which applies a user-defined bi-operand function on the Ciphertext with a tree-like combination.
    - Changes to the Keys structs:
        - Added `EvaluationKeySet`, which enables users to provide custom loading/saving/persistence policies and implementation for the `EvaluationKeys`.
        - `SwitchingKey` has been renamed `EvaluationKey` to better convey that these are public keys used during the evaluation phase of a circuit. All methods and variable names have been renamed accordingly.
        - The struct `RotationKeySet` holding a map of `SwitchingKeys` has been replaced by the struct `GaloisKey` holding a single `EvaluationKey`.
        - The `RelinearizationKey` type now stores a single GSW-like encryption of `s^2`, which is what the schemes' relinearization methods currently support.
    - Changes to the `KeyGenerator`:
        - The `NewKeyGenerator` returns a `*KeyGenerator` instead of an interface.
        - Simplified the `KeyGenerator`: methods to generate specific sets of `rlwe.GaloisKey` have been removed. Instead, the corresponding method on `rlwe.Parameters` allows to get the appropriate `GaloisElement`s.
        - Improved the API consistency of the `rlwe.KeyGenerator`. Methods that allocate elements have the suffix `New`. Added corresponding in-place methods.
        - It is now possible to generate `rlwe.EvaluationKey`, `rlwe.GaloisKey` and `rlwe.RelinearizationKey` at specific levels (for both `Q` and `P`) and with a specific `BaseTwoDecomposition` by passing the corresponding pre-allocated key.
    - Changes to the `MetaData`:
        - Content of the `MetaData` struct is now divided into `PlaintextMetaData` and `CiphertextMetaData`.
        - `PlaintextMetaData` contains the fields:
            - `Scale`
            - `LogDimensions`: represents the concept of plaintext algebra dimensions (e.g. BGV/BFV = [2, n] and CKKS = [1, n/2])
            - `IsBatched`: Boolean indicating if the plaintext is batched or not.
        - `CiphertextMetaData` contains the fields:
            - `IsNTT`: Boolean indicating whether the ciphertext is in the NTT domain.
            - `IsMontgomery`: Boolean indicating whether the ciphertext is in the Montgomery domain.
    - Changes to the tests:
        - Added accurate noise bounds for the tests.
        - Substantially increased the test coverage of `rlwe` (for both the amount of operations and parameters).
        - Substantially increased the number of benchmarked operations in `rlwe`.
    - Other changes:
        - Added generic `Element[T]` which serves as a common underlying type for ciphertext types.
        - The argument `level` is now optional for `NewCiphertext` and `NewPlaintext`.
        - `EvaluationKey` (and all parent structs) and `GadgetCiphertext` now take an optional argument `rlwe.EvaluationKeyParameters` that allows to specify the level `Q` and `P` and the `BaseTwoDecomposition`.
        - Allocating zero `rlwe.EvaluationKey`, `rlwe.GaloisKey` and `rlwe.RelinearizationKey` now takes an optional struct `rlwe.EvaluationKeyParameters` specifying the levels `Q` and `P` and the `BaseTwoDecomposition` of the key.
        - Changed `[]*ring.Poly` to `structs.Vector[ring.Poly]` and `[]ringqp.Poly` to `structs.Vector[ringqp.Poly]`.
        - Replaced the struct `CiphertextQP` by `Element[ringqp.Poly]`.
        - Added basic interface description for `Parameters`, `Encryptor`, `PRNGEncryptor`, `Decryptor`, `Evaluator` and `PolynomialEvaluator`.
        - All structs that can be serialized now implement the method V Equal(V) bool.
        - Setting to negative values the Hamming weight of the secret or the standard deviation of the error through `NewParameters` will instantiate these fields as zero values and return a warning (as an error).
- DRLWE:
    - The package `drlwe` has been renamed `mhe`.
    - Renamed:
            - `NewCKGProtocol` to `NewPublicKeyGenProtocol`.
            - `NewRKGProtocol` to `NewRelinKeyGenProtocol`.
            - `NewCKSProtocol` to `NewGaloisKeyGenProtocol`.
            - `NewRTGProtocol` to `NewKeySwitchProtocol`.
            - `NewPCKSProtocol` to `NewPublicKeySwitchProtocol`.
    - Replaced `[dbfv/dbfv/dckks].MaskedTransformShare` by `drlwe.RefreshShare`.
    - Added `EvaluationKeyGenProtocol` to enable users to generate generic `rlwe.EvaluationKey` (previously only the `GaloisKey`).
    - It is now possible to specify the levels of the modulus `Q` and `P`, as well as the `BaseTwoDecomposition` via the optional struct `rlwe.EvaluationKeyParameters`, when generating `rlwe.EvaluationKey`, `rlwe.GaloisKey` and `rlwe.RelinearizationKey`.
    - Arbitrarily large smudging noise is now supported.
    - Fixed `CollectiveKeySwitching` and `PublicCollectiveKeySwitching` smudging noise to not be rescaled by `P`.
    - Tests and benchmarks in package other than the `RLWE` and `DRLWE` packages that were merely wrapper of methods of the `RLWE` or `DRLWE` have been removed and/or moved to the `RLWE` and `DRLWE` packages.
    - Improved the GoDoc of the protocols.
    - Added accurate noise bounds for the tests.
- DBFV:
    - The package `dbfv`, which was merely a wrapper of the package `dbgv`, has been removed.
- DBGV:
    - The package `dbgv` has been renamed `mheint` and moved to `mhe/mheint`.
- DCKKS:
    - The package `dckks` has been renamed `mhefloat` and moved to `mhe/mhefloat`.
- RGSW:
    - The package `rgsw` has been moved to `core/rgsw`.
    - Expanded the encryptor to be able encrypt from an `rlwe.PublicKey`.
    - Added tests for encryption and external product.
- RING: 
    - Changes to sampling:
        - Updated Gaussian sampling to work with arbitrary size standard deviation and bounds.
        - Added a generic `Sampler` interface.
    - Added finite field polynomial interpolation.
    - Re-enabled NTT for ring degree smaller than 16.
    - Replaced  `Log2OfInnerSum` by `Log2OfStandardDeviation` in the `ring` package, which returns the log2 of the standard deviation of the coefficients of a polynomial.
    - Renamed `Permute[...]` by `Automorphism[...]` in the `ring` package.
    - Added non-NTT `Automorphism` support for the `ConjugateInvariant` ring.
    - Replaced all prime generation methods by `NTTFriendlyPrimesGenerator` which provides a more user friendly API and better functionality.
    - Added large standard deviation sampling.
    - Refactoring of the `ring.Ring` object:
        - The `ring.Ring` object is now composed of a slice of `ring.SubRings` structs, which store the pre-computations for modular arithmetic and NTT for their respective prime.
        - The methods `ModuliChain`, `ModuliChainLength`, `MaxLevel`, `Level` have been added to the `ring.Ring` type. 
        - Added the `BinaryMarshaller` interface implementation for `ring.Ring` types. It marshals the factors and the primitive roots, removing the need for factorization and enabling a deterministic ring reconstruction.
        - Removed all methods with the API `[...]Lvl(level, ...)`. Instead, to perform operations at a specific level, a lower-level `ring.Ring` type can be obtained using `ring.Ring.AtLevel(level)` (which is allocation-free).
        - Subring-level methods such as `NTTSingle` or `AddVec` are now accessible via `ring.Ring.SubRing[level].Method(*)`. Note that the consistency changes across method names also apply to these methods. For example, `NTTSingle` and `AddVec` are now simply `NTT` and `Add` when called via a `SubRing` object.
        - Updated `ModDownQPtoQNTT` to round the RNS division (instead of flooring).
        - The `NumberTheoreticTransformer` interface no longer has to be implemented for arbitrary `*SubRing` and it abstracts this parameterization as its instantiation.
        - The core NTT method now takes `N` as an input, enabling NTT of different dimensions without having to modify the internal value of the ring degree in the `ring.Ring` object.
- UTILS: 
    - Updated methods with generics when applicable.
    - Added public factorization methods `GetFactors`, `GetFactorPollardRho` and `GetFactorECM`.
    - Added subpackage `sampling` which regroups the various random bytes and number generator that were previously present in the package `utils`.
    - Added the package `utils/bignum` which provides arbitrary precision arithmetic, tools to create and evaluate polynomials, and tools to perform polynomial approximations of functions, notably Chebyshev and Multi-Interval Minimax approximations.
    - Added subpackage `buffer` which implements custom methods to efficiently write and read slices on any writer or reader implementing a subset interface of the `bufio.Writer` and `bufio.Reader`.
        - Added `Writer` interface and methods to write specific objects on a `Writer`.
        - Added `Reader` interface and methods to read specific objects from a `Reader`.
        - Added `RequireSerializerCorrect` which checks that an object satisfies `io.WriterTo`, `io.ReaderFrom`, `encoding.BinaryMarshaler` and `encoding.BinaryUnmarshaler`, and that these interfaces are correctly implemented.
    - Added subpackage `structs`:
        - New structs:
            - `Map[K constraints.Integer, T any] map[K]*T`.
            - `Matrix[T any] [][]T`.
            - `Vector[T any] []T`.
        - All the above structs satisfy the following interfaces:
            - `(T) CopyNew() *T`.
            - `(T) BinarySize() (int)`.
            - `(T) WriteTo(io.Writer) (int64, error)`.
            - `(T) ReadFrom(io.Reader) (int64, error)`.
            - `(T) MarshalBinary() ([]byte, error)`.
            - `(T) UnmarshalBinary([]byte) (error)`.
            - `(T) Equal(T) bool`.
    
## [4.1.0] - 2022-11-22 
- Further improved the generalization of the code across schemes through the `rlwe` package and the introduction of a generic scale management interface.
- All: uniformized the `prec` type to `uint` for `*big.Float` types.
- All: renamed `WriteTo<32/64>` to `Encode<32/64>` and `DecodePoly<32/64>` to `Decode<32/64>`, added similar method to `rlwe.Ciphertext`.
- RLWE: added the type `rlwe.Scale`, which is now a field in the `rlwe.Parameters`.
- RLWE: added the struct `MedaData` which stores the `Scale`, and boolean flags `IsNTT` and `IsMontgomery`. 
- RLWE: added the field `MetaData` to the `rlwe.Plaintext`, `rlwe.Ciphertext`, `rlwe.CiphertextQP`.
- RLWE: added `DefaultScale` and `DefaultNTTFlag` to the `rlwe.ParametersLiteral` struct. These are optional fields which are automatically set by the respective schemes.
- RLWE: elements from `rlwe.NewPlaintext(*)` and `rlwe.NewCiphertext(*)` are given default `IsNTT` and `Scale` values taken from the `rlwe.Parameters`, which depend on the scheme used. These values can be overwritten/modified manually.
- RLWE: added `logGap` parameter to `Evaluator.Expand`, which enables to extract only coefficients whose degree is a multiple of `2^logGap`.
- BFV: the level of the plaintext and ciphertext must now be specified when creating them.
- CKKS: significantly reduced the pre-computation time of the roots, especially for the arbitrary precision encoder.
- CKKS/BGV: abstracted the scaling factor, using `rlwe.Scale`. See the description of the struct for more information.
- BFV/BGV: added the flag `-print-noise` to print the residual noise, after decryption, during the tests.
- BFV/BGV/CKKS: added scheme specific global constant `DefaultNTTFlag`.
- BFV/BGV/CKKS: removed scheme-specific ciphertexts and plaintexts types. They are replaced by generic `rlwe.Ciphertext` and `rlwe.Plaintext`.
- BFV/BGV/CKKS: removed scheme-specific `KeyGenerator`, `Encryptor` and `Decryptor`. They have been replaced by `rlwe.KeyGenerator`, `rlwe.Encryptor` and `rlwe.Decryptor`. The API go instantiate those struct from the scheme specific API, e.g. `bgv.NewEncryptor`, is still available but will return its corresponding `rlwe` struct.
- BFV/BGV/CKKS: removed the following deprecated methods, when applicable
    - `AddNoMod`, `AddNoModNew`, `SubNoMod`, `SubNoModNew`, `Reduce` and `ReduceNew`
    - `PowerOf2`, `Power` and `PowerNew` which are replaced by `PolynomialBasis` and `GenPower`.
- BFV/BGV/CKKS: the naive method algorithms for `InnerSum` and `Replicate` have been removed. The method names `InnerSumLog` and `ReplicateLog` have been replaced by `InnerSum` and `Replicate` respectively.

## [4.0.0] - 2022-10-04
- Added BGV/DBGV schemes.
- ALL: added default parameters for LogN=11 and LogN=10.
- RING: prime generation no longer skips the first candidate.
- RING: reworked marshalling of `ring.Poly` object. The new available methods are:
    - `ring.Poly` now has a `.Buff` 1-dimensional slice which is the only heavy allocation of a `ring.Poly`. The `.Coeffs` 2-dimensional slice is a re-slicing of `.Buff`.
    - `GetDataLen64` and `GetDataLen32`: gets the length in bytes of an encoded `ring.Poly` object.
    - `WriteTo64` and `WriteTo32`: encodes a `ring.Poly` object on a pre-allocated slice of bytes.
    - `WriteCoeffsTo64` and `WriteCoeffsTo32`: encodes a slice of coefficients on a pre-allocated slice of bytes.
    - `DecodeCoeffs64` and `DecodeCoeffs32`: decodes a slice of bytes on a slice of coefficients.
    - `DecodePoly64` and `DecodePoly32`: decodes a slice of bytes on a pre-allocated `ring.Poly` object.
- RING: renamed `ring.Poly.Degree()` to `ring.Poly.N()` for consistency.
- RING: removed `ring.Poly.LenModuli()` deprecated method.
- RING: changed `ring.NewPoly` to take the `level` as argument instead of the number of moduli, for consistency.
- RLWE: added several types of ciphertexts:
    - `rlwe.CiphertextQP` represents a ciphertext that is encrypted in the extended ring R_QP.
    - `rlwe.GadgetCiphertext` represents an encryption in the extended ring R_QP of a plaintext that is decomposed in the CRT and power-of-two basis (e.g., public switching keys).
- RLWE: changed representation of `rlwe.PublicKey` types which are now stored in Montgomery form, consistently with all other key types.
- RLWE: changed `rlwe.SwitchingKey` type to use `rlwe.GadgetCiphertext` internally.
- RLWE: generalized `rlwe.KeySwitcher` into `rlwe.Evaluator`, which provides new functionalities:
    - `DecomposeNTT`: decomposes a polynomial modulo the special RNS basis and extends its basis from Q to QP.
    - `DecomposeSingleNTT`: decomposes a polynomial modulo a single power of the special RNS basis and extends its basis from Q to QP.
    - `ExpandRLWE`: extracts each coefficient of a RLWE sample to the degree-0 coefficient of multiple RLWE samples.
    - `MergeRLWE`: merges the degree-0 coefficient of multiple RLWE samples into a single RLWE sample.
    - `GadgetProduct`: evaluates `ring.Poly x gadget.Ciphertext -> RLWE`, where `gadget.Ciphertext` is a matrix of RLWE samples encrypting scaled plaintext by the special RNS basis and a modulus P.
    - `GadgetProductNoModDown`: evaluates `ring.Poly x gadget.Ciphertext -> RLWE` but without the division by P (the result is given mod QP).
    - `GadgetProductSinglePAndBitDecompNoModDown`: evaluates `ring.Poly x gadget.Ciphertext -> RLWE`, where `gadget.Ciphertext` is a matrix of RLWE samples encrypting scaled plaintext by the special RNS basis along with a base-2 basis and an optional prime P.
    - `Relinearize`: reduces the degree of a `rlwe.Ciphertext` to one by homomorphically evaluating the decryption of the higher-degree terms.
    - `KeySwitch`: homomorphically re-encrypts a `rlwe.Ciphertext` under a new secret.
    - `KeyswitchHoisted`: homomorphically re-encrypts a `rlwe.Ciphertext` under a series of new secrets, returning a new ciphertext for each secret.
    - `KeyswitchHoistedNoModDown`: homomorphically re-encrypts a `rlwe.Ciphertext` under a series of new secrets, returning a new ciphertext for each secret, but without the division by P (the result is given mod QP).
    - `Automorphism`: homomorphically evaluates the map `X -> X^k`.
    - `AutomorphismHoisted`: homomorphically evaluates multiple maps of the type `X -> X^k`, returning a new ciphertext for each map.
    - `AutomorphismHoistedNoModDown`: homomorphically evaluates multiple maps of the type `X -> X^k`, returning a new ciphertext for each map, but without the division by P (result is given mod QP).
    - `Trace`: homomorphically evaluates the map `X -> sum((-1)^i * X^{i*n+1}) for n <= i < N`.
    - `ExternalProduct`: evaluates `rlwe.Ciphertext x rgsw.Ciphertext -> rlwe.Ciphertext`.
- RLWE: re-enabled bit-decomposition, on top of RNS decomposition, for the inner-product between `rlwe.Ciphertext` and `gadget.Ciphertext`.
    - This functionality can be enabled by setting `Pow2Base` to the desired power of two basis.
    - This functionality can be used in conjunction with the RNS hybrid decomposition (with a modulus `P`) only when `P` is composed of a single prime.
    - This functionality is disabled if `Pow2Base` is set to zero (default value).
- RLWE: enabled instantiation of `rlwe.Parameters` without the modulus `P`.
- RLWE: revamped the `rlwe.Encryptor` interface and implementing structs:
    - Added the `.EncryptZero` method to generate encryptions of zeros.
    - The `.Encrypt` and `.EncryptZero` now accept `ct interface{}` as their ciphertext argument and determine the type of encryption to be performed according to the runtime type of `ct`.
- RLWE: added the `PRNGEncryptor` type, which supports secret-key encryption from a user-specified PRNG.
- RLWE: `rlwe.KeyGenerator` now uses an `rlwe.Encryptor` internally, to generate secret keys, encryption keys and evaluation keys.
- RLWE: extracted the `rlwe/ringqp` sub-package which provides the `ringqp.Ring` and `ringqp.Poly` types to respectively replace the former types `rlwe.RingQP` and `rlwe.PolyQP`.
- DRLWE: added the `Thresholdizer` and `Combiner` types for t-out-of-N-threshold schemes through Shamir secret-sharing.
- DRLWE: added a `README.md` providing package overview and usage instructions.
- DRLWE: removed the obsolete `CollectivePublicKeyGenerator`, `RelinearizationKeyGenerator`, `RotationKeyGenerator`, `PublicKeySwitchingProtocol` and `KeySwitchingProtocol` interfaces.
- DRLWE: renamed `AggregateShare` methods to `AggregateShares`.
- RGSW: added package `rgsw`, which provides a partial implementation of the RLWE-based RGSW encryption scheme. This includes:
    -  `rgsw.Encryptor` and the `rgsw.Ciphertext` types.
    -  `rgsw.Evaluator` to support the external product `RLWE x RGSW -> RLWE`.
    -  `rgsw/lut` sub-package that provides evaluation of Look-Up-Tables (LUT) on `rlwe.Ciphertext` types.
- BFV: renamed `Encoder.DecodeRingT` to `Encoder.SwitchToRingT` to better reflect the purpose of the method.
- CKKS: fixed `MulAndAdd` correctness for non-identical inputs.
- CKKS: added `advanced.EncodingMatrixLiteral.RepackImag2Real` optional field to repack the imaginary part into the right n real slots.
- CKKS: `Trace` now only takes as input the `logSlots` of the encrypted plaintext.
- CKKS: replaced the public variable `.Scale` with `.scale`, it can now be accessed with `.Scale()` and set to a new value with `.SetScale()`.
- CKKS: renamed the methods `ScalingFactor` and `SetScalingFactor` of the interface `Operand` to `Scale` and `SetScale` respectively.
- CKKS/bootstrapping: renamed method `Bootstrapp` to `Bootstrap`.
- BFV/CKKS: key-switching functionalities (such as rotations, relinearization and key-switching) are now all based on the `rlwe.Evaluator`.
- BFV/CKKS: the parameters now are based on the sub-type `rlwe.Parameters`.
- BFV/CKKS: removed deprecated methods `EncryptFromCRP` and `EncryptFromCRPNew`, users should now use the `PRNGEncryptor` interface.
- BFV/CKKS: fixed a panic happening during the benchmark testing.
- DBFV/DCKKS: removed the `dbfv/dckks.CKGProtocol`, `dbfv/dckks.RKGProtocol` and `dbfv/dckks.RTGProtocol` types. Users should use the corresponding `drlwe` types instead.
- DBFV/DCKKS: `MaskedTransformFunc` is now a struct and takes as additional input to the linear transform two Boolean flags to parameterize if the decoding/encoding process must be done before/after the linear transform.
- DBFV/DCKKS: `refresh` and `maskedTransform` protocols now allow the user to specify the output parameters, enabling parameter switching.
- DCKKS: fixed `dckks.RefreshProtocol` correctness when the output scale is different from the input scale.
- Examples: added `examples/ckks/advanced/lut`, which is an example that performs homomorphic decoding -> LUT -> homomorphic encoding on a `ckks.Ciphertext`.
- Examples: removed `examples/ckks/advanced/rlwe_lwe_bridge_LHHMQ20`, which is replaced by `examples/ckks/advanced/lut`.
- Examples: removed `examples/rlwe/lwe_bridge` since the code of this example is now part of `rlwe.Evaluator` and showcased in `examples/ckks/advanced/lut`.
- CI: revamped Makefile to no longer require github.com/dedis/coding and integrated linting/vet checks.

## [3.0.5]

- CKKS: Baby-Step Giant-Step Polynomial Evaluation Algorithm (BSGSPEA)
    - Added `PolynomialBasis`, a struct to generate powers of monomials. This struct can be marshalled.
    - Renamed former `PolynomialBasis` enumerated type to `BasisType`.
    - `EvaluatePoly` and `EvaluatePolyVector` now both accept pre-computed `PolynomialBasis` as input in addition to `Ciphertext`.
    - Fixed correctness error and panic when a non-relinearized ciphertext and a plaintext were given to `Mul` and `MulAndAdd`.
    - Fixed automatic-scale matching in BSGS that wasn't reliably ensuring that scales between two ciphertext to be added was the same.
    - Improved BSGSPEA with lazy relinearization and lazy rescaling.
    - Overall the precision of the BSGSPEA is greatly improved and its complexity is reduced. This also improves the precision of the bootstrapping.

## [3.0.4] - 2022-04-26

- CKKS: updated the bootstrapping circuit to use the key-encapsulation mechanism of `Bootstrapping for Approximate Homomorphic Encryption with Negligible Failure-Probability by Using Sparse-Secret Encapsulation`. The previous bootstrapping circuit can be run by setting `EphemeralSecretWeight=0`.
- BFV: added the `Evaluator.Rescale` and `Evaluator.RescaleTo` methods to switch BFV ciphertexts to lower levels.
- BFV: all `Evaluator` methods on ciphertext support all arithmetic operations at lower levels, but require that operands are at the same level.
- BFV: the plaintext modulus `T` can now equal to the level-zero modulus Q[0] (i.e., be a factor of the ciphertext modulus `Q`).
- BFV: added the methods `NewCiphertextLvl`, `NewPlaintextLvl`, `NewPlaintextMulLvl`, `Evaluator.AddScalar` and `Evaluator.MulScalarAndAdd`. 
- BFV: merged `[]uint64` and `[]int64` plaintext encoding methods (e.g. `EncodeUint` and `EncodeInt` are replaced by `Encode`) and added the respective `[...]New` methods.
- BFV: added the methods `EvaluatePoly` and `EvaluatePolyVector` for homomorphic polynomial evaluation.
- BFV/RING: moved `RNSScaler` from `ring` to `bfv`.
- RING: removed deprecated `SimpleScaler`.

## [3.0.2] - 2022-02-21

- RING: fixed sparse ternary sampler to properly sample on non-zero poly.

## [3.0.1] - 2022-02-21

- RLWE/CKKS/BFV: added the `H` field and `HammingWeight` method in parameters-related structs, to specify distribution of all secrets in the schemes.
- RLWE/DRLWE: all secrets in the ternary distribution are now sampled with a fixed hamming weight, according to the parameters.
- CKKS: encoder is now about 3.5x faster (without taking the NTT into account).

## [3.0.0] - 2022-02-21
- ALL: renamed the module to `github.com/tuneinsight/v3`.
- RING: renamed `FastBasisExtender` to `BasisExtender`.
- RING: `.PolyToBigint[...](*)` now take as input `gap` which defines the multiples of `X^{i*gap}` to reconstruct.
- RLWE: removed `FastEncryptor`. Encryption without rescaling by `P` is now automatically used by `Encryptor` if no `P` is specified in the parameters.
- RLWE: `NewAdditiveShareBigint` now takes as input the size of the share.
- RLWE/CKKS/BFV: added `.ShallowCopy()`, `.WithKey()` (shallow copy with new key) to `Encryptor` and `Decryptor`.
- BFV/CKKS: added `.ShallowCopy()` to `Encoder` and `EncoderBigComplex` (only CKKS).
- DRLWE/DCKKS/DBFV: added `.ShallowCopy()` to all protocols.
- DLRWE/DCKKS/DBFV: protocols `drlwe.CKSProtocol` and `drlwe.PCKSProtocol` and sub-protocols based on these two protocols now only take a polynomial as input for the share generation instead of the full ciphertext.
- DRLWE/DCKKS/DBFV: uniformized API of share generation and aggregation to `.GenShare(*)` and `.AggregateShare(*)`.

## [2.4.0] - 2022-01-10

- RING: added support for ring operations over the conjugate invariant ring.
- RING: added support for custom NTT via the `NumberTheoreticTransformer` interface.
- RLWE: added support for RLWE primitives over the conjugate invariant ring.
- RLWE: added `encoding.BinaryMarshaler` implementation for `rlwe.Ciphertext` types.
- RLWE: added an example implementation of homomorphic RLWE slot shuffling based on RLWE<->LWE conversion.
- RLWE: increased the maximum supported polynomial degree to 2^17.
- CKKS: Trace does not multiply the output by (N/n)^-1 anymore.
- CKKS: added support for the CKKS scheme over the conjugate invariant ring.
- CKKS: renamed `Scale` to `DefaultScale` in `Parameters` and `ParametersLiteral`.
- CKKS: added the `Evaluator.Average` method.
- CKKS: added `DomainSwitcher` type for conversion between Standard and Conjugate Invariant variants of CKKS.
- CKKS: added support for both  `[]complex128` and `[]float64` as input to `Encoder.Encode*` methods.
- CKKS: added support for `[]float64` as input to `GetPrecisionStats`.
- CKKS: added support for `func(float64)float64` and `func(complex128)complex128` as input to `Approximate`.
- CKKS: uniformized the arguments' position for all methods of the `Encoder` interface.
- CKKS: renamed `Encoder.EncodeNTT/New` to `Encoder.Encode/New`and added  `Encoder.EncodeSlots`, `Encoder.DecodeSlots` and `Encoder.DecodeSlotsPublic`.
- CKKS: added `EncodeSlotsQP` to encode on `rlwe.PolyQP` to support the new `LinearTransform` interface.
- CKKS: improved `Encoder` implementation; it is now much faster when encoding sparse plaintexts.
- CKKS: changed the approximation intervals from `complex128` to `float64`.
- CKKS: renamed `PtDiagMatrix` to `LinearTransform`.
- CKKS: added `LinearTransform.Rotations()` to get the required rotation for the receiver plaintext linear tranform.
- CKKS: added `Parameters.RotationsForLinearTransform` to get the required rotation for the given plaintext linear tranform.
- CKKS: added `NewLinearTransform`, `EncodeNewLinearTransform`, `GenLinearTransform` and `GenLinearTransformBSGS` to allocate and initialize plaintext linear transforms.
- CKKS: removed plaintext linear transforms (old `PtDiagMatrix`) constructors and initializers from `Encoder`.
- CKKS: added `Evaluator.EvaluatePolyVector` to enable efficient evaluation of multiple different polynomials on the same ciphertext.
- CKKS: fixed a bug in the BSGS approach for linear transform where the selection of the ratio bettween giant step and baby step could lead to a ratio of N.
- CKKS: the EvalMod step of the bootstrapping now works for moduli of any size, regardless of `Q[0]` or `MessageRatio`.
- DCKKS: added support for multiparty CKKS over the conjugate invariant ring.
- DCKKS: fixed `MaskedTransformProtocol` correctness for sparse plaintexts.
- Examples: updated the `ckks/sigmoid` example to `ckks/polyeval` example, that now showcases the use of `PolynomialVector`.

## [2.3.0] - 2021-10-12

- RING: added `MapSmallDimensionToLargerDimensionNTT` method which maps from  Y = X^{N/n} to X in the NTT domain.
- RING: `FastBasisExtender` type can now extend the basis of polynomials of any level in base Q to polynomials of any level in base P.
- RING: changed RNS division `Div[floor/round]BylastModulus[NTT]` to `Div[floor/round]BylastModulus[NTT]Lvl` (the level of the last modulus must always be provided).
- RING: RNS division no longer modifies the output polynomial's level, this is to facilitate the usage of memory pools.
- RING: added the method `MFormVector`, which switches a slice of `uint64` into the Montgomery domain.
- RING: RNS scaler (used in BFV) does not modify the input anymore.
- RLWE: `GenSwitchingKey` now accepts secret-keys of different dimensions and level as input to enable re-encryption between different ciphertext degrees.
- RLWE: added `SwitchCiphertextRingDegreeNTT` and `SwitchCiphertextRingDegree` to switch ciphertext ring degrees.
- RLWE: added the `rlwe.RingQP` type to represent the extended ring R_qp.
- RLWE: added the `rlwe.PolyQP` type to represent polynomials in the extended ring R_qp.
- DRLWE: added the `CKGCRP`, `RKGCRP`, `RTGCRP` and `CKSCRP` types to represent the common reference polynomials in these protocols.
- DRLWE: added the `CRS` interface for PRNGs that implement a common reference string among the parties.
- DRLWE: added the `SampleCRP(crs CRS)` method to each protocol types to sample their respective CRP type.
- BFV: changed the plaintext scaling from `floor(Q/T)*m` to `round((Q*m)/T)` to reduce the initial ciphertext noise. 
- CKKS: added the `ckks/advanced` sub-package and moved the homomorphic encoding, decoding and modular reduction into it.
- CKKS: added the `ckks/bootstrapping` sub-package and moved the CKKS bootstrapping into it. This package now mostly relies on the `ckks/advanced` package.
- CKKS: renamed the `ChebyshevInterpolation` type to `Polynomial`.
- CKKS: removed the `EvaluateCheby` method that was redundant with the `EvaluatePoly` one.
- CKKS: optimized the `EvaluatePoly` to account for odd/even polynomials and fixed some small imprecisions in scale management occurring for some specific polynomial degrees.
- CKKS: some advanced methods related to automorphisms are now public to facilitate their external use.
- CKKS: improved the consistency of the API for in-place and `[..]New` methods.
- CKKS: added the method `NewCiphertextAtLevelFromPoly`, which creates a ciphertext at a specific level from two polynomials.
- CKKS: updated precision stats struct, added L2 norm in the statistics and improved the command line prints. 
- CKKS: improved the algorithmic complexity of `MultiplyByDiagMatrixBSGS` and updated the bootstrapping parameters accordingly.
- CKKS: `PermuteNTTHoistedNoModDown` now returns `[phi(P*c0 + c0'), phi(c1')]` instead of `[phi(c0'), phi(c1')]`.
- CKKS: Changed `RotateHoistedNoModDown` to `RotateHoistedNoModDownNew` for consistency.
- DBFV/DCKKS: both now use their respective CRP type for each protocol.
- EXAMPLE: added showcase of the `ckks/advanced` sub-package: a bridge between CKKS and FHEW ciphertexts using homomorphic decoding, ring dimension switching, homomorphic matrix multiplication and homomorphic modular reduction.

## [2.2.0] - 2021-07-15

- Added SECURITY.md
- ALL: when possible, public functions now use `int` instead of `uint64` as parameters and return values.
- ALL: `ring.Ring` are not instantiated once in the parameters and read only. They are then accessed by other structs, like the encryptor or evaluator.
- RING: removed `MulPoly` and its related tests.
- RING: `ring.Ring` is now read-only and thread-safe.
- RING: RNS rescaling API is now in place and can take a different poly as output.
- RING: added `ReadFromDistLvl` and `ReadAndAddFromDistLvl` to Gaussian sampler API.
- RING: added `IsNTT` and `IsMForm` flags in the `ring.Poly` type. For now, these flags are never checked or changed by the `ring` package.
- RLWE: added a new `rlwe` package as common implementation base package for the Lattigo RLWE schemes.
- RLWE: extracted the `rlwe.Parameters` type as common base struct for BFV and CKKS parameters.
- RLWE: extracted the `rlwe.KeyGenerator` type as common key-generator for BFV and CKKS.
- RLWE: extracted the `rlwe.Ciphertext` type as common base struct for BFV and CKKS ciphertexts.
- RLWE: extracted the `rlwe.Plaintext` type as common base struct for BFV and CKKS plaintext.
- RLWE: extracted the `rlwe.Encryptor`  type as common base interface for BFV and CKKS encryptors.
- RLWE: extracted the `rlwe.Decryptor`  type as common base interface for BFV and CKKS decryptors.
- RLWE: extracted the `rlwe.KeySwitcher` type as a common key-switching implementation for BFV and CKKS evaluators.
- RLWE: renamed the `Parameters.Copy()` method to `Parameters.CopyNew()` for consistency.
- RLWE: added `Parameter` struct, that stores the relevant `ring.Ring` instances and has getter methods to access them.
- RLWE: added equality and inclusion check methods for the `rlwe.RotatationKeySet` type.
- RLWE: added tests for encryption, decryption, key-generation and key-switching.
- RLWE: moved keys related marshalling tests of `bfv` and `ckks` packages the `rlwe` package.
- DRLWE: added a new `drlwe` package as a common implementation base for the lattigo multiparty RLWE schemes.
- DRLWE: added tests for the protocols.
- DRLWE: moved keys-related marshalling tests of `dbfv` and `dckks` packages to the `drlwe` package.
- BFV/CKKS: the schemes now use a common implementation for their keys.
- BFV/CKKS: the rotation-keys are now indexed by their corresponding Galois automorphism.
- BFV/CKKS: the `Evaluator` interface now has a single method for all column rotations and one method for the row-rotation/conjugate.
- BFV/CKKS: the relinearization and rotation keys are now passed to the `Evaluator` constructor methods (and no longer to the operations methods).
- BFV/CKKS: added the ParameterLiteral type for literally specifying scheme parameters in Go programs.
- BFV/CKKS: removed the now obsolete `Moduli` and `LogModuli` types and their associated `Parameters` constructors.
- BFV/CKKS: `Parameters` types are now passed by value in most situations.
- BFV/CKKS: added `encoding/json`-compatible JSON serializers and deserializers for the `Parameters` types.
- BFV/CKKS: removed the scheme-specific key types.
- BFV/CKKS: added a `-params=[params json]` flag for all test and bench suites for specifying parameters from the command line.
- DBFV/DCKKS: added a common interface and implementation for each multiparty protocol.
- DBFV/DCKKS: added standalone Encryption-To-Shares (`E2SProtocol`) and Shares-To-Encryption (`S2EProtocol`) protocols for domain switching between encryptions and secret-shares.
- DBFV/DCKKS: generalized the Refresh-and-permute protocol into generic `MaskedTransformProtocol` that accepts an arbitrary linear function.
- DCKKS: public-refresh now takes a target desired output scale, which enables refreshing the ciphertext to the default scale.
- BFV: the moduli of `ringQMul` are now generated based on `N` and`Q`.
- CKKS: added `Parameter` methods that compute the required rotations for relevant `Evaluator` operations.
- CKKS: added methods for performing linear-transformations and improved several aspects listed below.
- CKKS: improved the tests for `CoeffsToSlots` and `SlotsToCoeffs`.

#### CKKS Bootstrapping
- The procedure now allows for a more granular parameterization.
- Added flag in bootstrapping parameters for bit-reversed inputs (with bit-reversed output) CoeffsToSlots and SlotsToCoeffs.
- Added optional Arcsine.
- The procedure now uses the new linear-transformation API.
- `CoeffsToSlots` and `SlotsToCoeffs` are now standalone public functions.

#### New CKKS Evaluator methods 
- `RotateHoisted`: evaluates several rotations on a single ciphertext.
- `LinearTransform`: evaluates one or more `PtDiagMatrix` on a ciphertext using `MultiplyByDiagMatrix` or `MultiplyByDiagMatrixBSGS` according to the encoding of `PtDiagMatrix`.
- `MultiplyByDiagMatrix`: multiplies a ciphertext with a `PtDiagMatrix` using n rotations with single hoisting.
- `MultiplyByDiagMatrixBSGS`: multiplies a ciphertext with a `PtDiagMatrix` using 2sqrt(n) rotations with double-hoisting.
- `InnerSumLog`: optimal log approach that works for any value (not only powers of two) and can be parameterized to inner sum batches of values (sub-vectors).
- `InnerSum`: naive approach that is faster for small values but needs more keys.
- `ReplicateLog`: optimal log approach that works for any value (not only powers of two) and can be parameterized to replicate batches of values (sub-vectors).
- `Replicate`: naive approach that is faster for small values but needs more keys.

#### New CKKS Encoder methods
- `PtDiagMatrix`: struct that represents a linear transformation.
- `EncodeDiagMatrixBSGSAtLvl`: encodes a `PtDiagMatrix` at a given level, with a given scale for the BSGS algorithm.
- `EncodeDiagMatrixAtLvl`: encodes a `PtDiagMatrix` at a given level, with a given scale for a naive evaluation.
- `DecodePublic`: adds Gaussian noise of variance floor(sigma * sqrt(2*pi)) before the decoding step (see SECURITY.md).
- `DecodeCoeffsPublic`: adds Gaussian noise of variance floor(sigma * sqrt(2*pi)) before the decoding step (see SECURITY.md).
- `GetErrSTDFreqDom` : get the error standard deviation in the frequency domain (slots).
- `GetErrSTDTimeDom`: get the error standard deviation in the time domain (coefficients).

#### CKKS Fixes
- `MultByi` now correctly sets the output ciphertext scale.
- `Relinearize` now correctly sets the output ciphertext level.
- matrix-vector multiplication now correctly manages ciphertexts of higher level than the plaintext matrix.
- matrix-vector encoding now properly works for negative diagonal indexes.

#### Others
- PrecisionStats now includes the standard deviation of the error in the slots and coefficients domains.

## [2.1.1] - 2020-12-23

### Added
- BFV/CKKS: added a check for minimum polynomial degree when creating parameters.
- BFV: added the `bfv.Element.Level` method.
- RING: test for sparse ternary sampler.

### Changed
- BFV/CKKS: pk is now (-as + e, a) instead of (-(as + e), a).
- BFV: harmonized the EvaluationKey setter from `SetRelinKeys` to `Set`.
- CKKS: renamed `BootstrappParams` into `BootstrappingParameters`.
- CKKS: the `Evaluator.DropLevel`, `Parameters.SetLogSlots` and `Element.Copy` methods no longer return errors.
- RING: minimum poly degree modulus is 16 to ensure the NTT correctness.
- RING: isPrime has been replaced by big.ProbablyPrime, which is deterministic for integers < 2^64.

### Fixed
- ALL: reduced cyclomatic complexity of several functions.
- ALL: fixed all instances reported by staticcheck and gosec excluding G103 (audit the use of unsafe).
- ALL: test vectors are now generated using the crypto/rand instead of math/rand package.
- ALL: fixed some unhandled errors.
- BFV/CKKS:  improved the documentation: documented several hard-coded values and fixed typos.
- RING: fixed bias in sparse ternary sampling for some parameters.
- RING: tests for the modular reduction algorithms are now deterministic.

## [2.1.0] - 2020-12-11

### Added
- BFV: special-purpose plaintext types (`PlaintextRingT` or `PlaintextMul`) for optimized ct-pt operations. See bfv/encoder.go and bfv/plaintext.go.
- BFV: allocation-free `Encoder` methods.
- RING: `GenNTTPrimes` now takes the value `Nth` (for Nth primitive root) as input rather than `logN`.

### Changed
- BFV: the `Encoder.DecodeUint64` and `Encoder.DecodeInt64` methods now take the output slice as argument.
- CKKS: API of `Evaluator.RotateColumns` becomes `Evaluator.Rotate`.
- CKKS: the change of variable in `Evaluator.EvaluateCheby` isn't done automatically anymore and the user must do it before calling the function to ensure correctness.
- CKKS: when encoding, the number of slots must now be given in log2 basis. This is to prevent errors that would induced by zero values or non power of two values.
- CKKS: new encoder API : `EncodeAtLvlNew` and `EncodeNTTAtLvlNew`, which allow a user to encode a plaintext at a specific level.

### Removed
- CKKS: removed method `Evaluator.EvaluateChebySpecial`.
- BFV: removed `QiMul` field from `bfv.Parameters`. It is now automatically generated.

## [2.0.0] - 2020-10-07

### Performance
- Global 1.5x speed-up across all arithmetic (this does not include sampling).

### Added
- BFV/CKKS: Added fast encryption (directly in Q without the rescaling by P).
- CKKS: Added full-RNS scale-invariant bootstrapping (<https://eprint.iacr.org/2020/1203>).
- CKKS: Added parameterized tests for a range of experiments.
- CKKS: Added arbitrary precision encoding/decoding.
- CKKS: Added scale invariant polynomial evaluation.
- CKKS: Added encoding/decoding for coefficient packing.
- CKKS The user can now choose to encode a plaintext in or out of the NTT domain (the latter option leads to slightly faster encryptions).
- CKKS: Added secret-key gen with error distribution.
- DBFV: Added collective refresh with arbitrary permutation/linear transformation.
- DCKKS: Added collective refresh with arbitrary permutation/linear transformation.
- RING: Added arbitrary precision complex arithmetic, including cos and sin functions.
- RING: Added polynomial interpolation.
- RING: Added polynomial inversion.
- RING: Extracted interface type Scaler for polynomial coefficient scaling.
- RING: Added type RNSScaler as an efficient, cross-platform implementation of the Scaler interface.

### Changed 
- ALL: all tests now use "require".
- BFV/CKKS: Now parameters without P can be used, but the key-switching is disabled.
- BFV/CKKS: Now parameters do not use methods to access internal values.
- BFV/CKKS: New rotations keys optimized for hoisting rotations of the form (-phi^{-1}(s1)a + phi(s0) + e, a).
- BFV: The Decoder uses the RNSScaler implementation of the Scaler interface to perform the t/Q rescaling.
- CKKS: Simplified the code of the hybrid key-switching (does not affect user experience).
- CKKS: The encoding/decoding operations at level 0 are now 500% faster.
- CKKS: The encoder now accepts slices of complex values with length equal to or smaller than the specified number of slots.
- RING: Improved primes finding.
- RING: All Gaussian sampling now uses Ziggurat sampling.
- RING: Revamped polynomial samplers to make them more memory efficient, consistent user friendly, and to enable parallel sampling.
- RING: The SimpleScaler type now use slightly slower but cross-platform big.Int/Float.
- UTILS: Complete revamp of the PRNG (Blake2b XOF), to make it more user friendly and consistent.

### Removed
- BFV/CKKS: Parameters API generation GenFromLogModuli() and GenFromModuli() have been removed and replaced by Gen().
- CKKS: EvaluatePolyFast(.) and EvaluatePolyEco(.) are replaced by EvaluatePoly(.).
- CKKS: EvaluateChebyFast(.) and EvaluateChebyEco(.) are replaced by EvaluatePolyCheby(.).
- CKKS: EvaluateChebyEcoSpecial(.) and EvaluateChebyFastSpecial(.) are replaced by EvaluatePolyChebySpecial(.).
- RING: The Float128 type was removed due to cross-platform incompatibility.

### Fixes
- BFV: Fixed multiplication that was failing when #Qi != #QMul.
- BFV: Fixed a mempool corruption when encrypting from SK.
- CKKS: The function mulrelin now always returns a fully reduced polynomial.
- CKKS: The encoder now correctly checks that the number of slots is a power of two.
- RING: Prevented a rare case of uint64 overflow during prime sampling.
- RING: Prevented a rare case where two identical primes could be returned when sampling primes.

## [1.3.1] - 2020-02-26
### Added
- BFV/CKKS: Added API for encrypting using a CRP (common reference polynomial).
- BFV/CKKS: Added API for encrypting faster (encrypts zero directly in Q instead of QP and does not need to divide by P).
- BFV/CKKS: Parameters can now be created without the modulus P. Doing so disables all key-switching operations.
- CKKS: Added tests for hoisted rotations.
- RING: Added benchmarks for a NTT using purely Barrett reduction for comparison purposes.
### Changed 
- BFV/CKKS: Changed the switching keys from (-as1 + (s0-s1) + e, a) to (-as1 + s0 + e, a). This does not affect the user experience as it only changes the internal behavior, which is kept consistent. However, Rotation and KeySwitching keys generated with older releases will induce wrong results and will need to be re-generated.
### Fixes
- BFV: Fixed EncryptFromSK that was not correctly wiping the memory pool before using it, which lead to back encryptions.
- BFV: Fixed an index out of bound error that would happen during the multiplication if #QMul > #Qi.
- CKKS: Removed some redundant operations in the hoisted rotations.
- CKKS: MulRelin now always returns a fully reduced ciphertext.
- DCKKS: PCKS and CKS now correctly set the scale of the output ciphertext to the scale of the input ciphertext.
- RING: Fixed GenerateNTTPrimes that could return twice the same prime if the initial value was prime.
- RING: The function context.UniformPoly now samples based on the number of moduli of the context rather than based on the input polynomial.

## [1.3.0] - 2019-12-20
### Added
- All schemes: New switching-keys and key-switching algorithm based on the concept presented in https://eprint.iacr.org/2019/688.pdf.
- All schemes: New marshaling interface for all structures.
- BFV/CKKS: New Parameters structs and API that enable a better customization and fine tuning for specific applications.
- CKKS: New API for hoisted rotations, faster than sequential rotations.
- DBFV/DCKKS: Added collective refresh of a ciphertext (decentralized bootstrapping).
- RING: Added Ziggurat sampling, available from the context.
- RING: Enabled dense and sparse ternary polynomials sampling directly from the context.
- RING: New API enabling "level"-wise polynomial arithmetic.
- RING: New API for modulus switching with flooring and rounding.
- UTILS: The package utils now regroups all the utility methods which were previously duplicated among packages.
### Removed
- BFV/CKKS/DBFV/DCKKS: Removed their respective context. Ring context remains public.
- All schemes: Removed key-switching with bit decomposition. This option will however be re-introduced at a later stage since applications using small parameters can be impacted by this change.
- BFV/CKKS/RING: Removed redundant/irrelevant tests and benchmarks.
- BFV: Removed context QP as it is now not used in the multiplication.
- BFV: Removed int encoder, now only batch encoding is supported.
- CKKS: Modulus switching is moved to the Ring package.
- RING: Removed the algorithms that needed Float128 during the BFV multiplication.
- RING: Removed most wrapping methods for bigInt, which are now replaced by the native math/big package.
- RING: Removed ternary sampler, which is now part of the context.
### Changed
- All schemes: Encryptor, Decryptor, Encoder, Evaluator, KeyGenerator are now interface types.
- All schemes: Improved Godoc and error strings.
- ALl schemes: Greatly reduced the number of methods that could return an error.
- All schemes: New tests and benchmarks with fully supported regex.
- All schemes: Coefficient-wise arithmetic using double slices is now substantially faster.
- BFV/CKKS: Changed the name of the underlying ring contexts. Q now represents the ciphertext modulus (with QMul being the extended ciphertext modulus for BFV) and QP represents modulus of the keys (where P is the special primes used during the new key-switching).
- BFV/CKKS/DBFV/DCKKS: The structures are now created using the parameters instead of the context.
- BFV: Quantization during multiplication does not use Float128 any more, resulting in a substantial speed improvement.
- BFV: BatchEncoder has been renamed Encoder.
- CKKS: The scale is now stored as a float64 instead of a power of 2.
- CKKS: Rounding is applied instead of flooring when a real value is converted to an integer value. This change affects the rescaling and the encoding.
- CKKS: Use of one context for all levels, instead of requiring one ring context per level.
- CKKS: New baby-step giant-step algorithm for evaluating polynomials in standard and Chebyshev basis.
- CKKS: Reduced the number of NTT needed during the encryption.
- CKKS: API for MultConst is now MultByConst.
- BFV/CKKS: New API for the rotation-keys generation.
- DBFV/DCKKS: Complete revamp of the API and interfaces enabling a much easier integration into larger systems.
- DBFV/DCKKS: Improved PCKS and CKS using the concept of the new key-switching technique which enables to reduces the added noise.
- DCKKS: All protocols work for ciphertexts at any levels.
- RING: Faster MulScalarBigint (now similar to MulScalar).
- UTILS: PRNG must be keyed to be forward secure.
### Fixes
- All packages: Corrected typos, godoc and golint.
- CKKS: ciphertext rotation now correctly sets the scale of the output ciphertext.
- DBFV/DCKKS: Correctness is now ensured when the same protocol instance is used to generate multiples shares.

## [1.2.0] - 2019-12-01
Internal version, merged with 1.3.0.

## [1.1.0] - 2019-10-01
### Added
- CHANGELOG.md file.
- BFV: New methods on bfvcontext to access information/stored variables.
- BFV: When creating a new batch encoder, an error will now be returned if the plaintext modulus does not allow NTT.
- CKKS: New methods on ckkscontext to access information/stored variables.
- CKKS: Added marshalling and tests for marshalling.
- CKKS: Added default parameters and parameters marshalling.
- CKKS (API change): encryption can now also be done with the secret-key.
- CKKS (API change): new separate struct for the encoder, that will store a small memory pool for temporary elements.
- BFV/CKKS: Operand interface.
- RING: New method MultByVector.
- RING (API change): New Ternary Sampler, which enables to specify the key distribution {-1, 0, 1} -> [(1-p)/2, p, (1-p)/2]; it is faster than the previous implementation.
- GoDoc for BFV, CKKS, Ring, DBFV, and DCKKS.
- README for BFV and CKKS.
- Code cleaning for all packages.
- Minor optimizations for all packages.

### Changed
- Updated README.md.
- BFV (API change): bfvcontext now only accepts as input a struct of the type Parameters, similar to the one used for DefaultParams.
- BFV (API change): Removed bfvcontext from ciphertexts and plaintexts.
- BFV (API change): Encryption can now also be done with the secret-key.
- BFV: The value logQ in the bfvcontext now stores the bit size of the product of the ciphertext's moduli.
- BFV: The information printed by the tests now better conveys the parameters.
- BFV: Updated and optimized default parameters.
- BFV: Updated example with secure parameters.
- CKKS (API change): ckkscontext now only accepts as input a struct of the type Parameters, similar to the one used for DefaultParams.
- CKKS (API change): ckkscontext now requires as input the moduli chain (in bit size), instead of a generic logQ and levels. This allows more fine-grained control on the rescaling and levels.
- CKKS (API change): removed ckkscontext from ciphertexts and plaintexts.
- CKKS: Updated the value logQ in the ckkscontext, which now stores the bit size of the product of the ciphertext's moduli.
- CKKS: Updated the information printed by the tests to more reflect the parameters.
- BFV/CKKS: Greatly simplified the code related to the rotations.
- BFV/CKKS: Reduced the number of instances where an error could be returned and updated information of the returning errors.
- RING (API change): Changed the API of the validation method of the context to better reflect its purpose.
- BFV/CKKS/RING: The copy function will now copy the input poly on the target poly (previously the target was copied on the input).
- Updates on all packages and tests to comply with the API changes in BFV and CKKS.

### Removed
- The evaluator of both BFV and CKKS cannot operate on two plaintexts anymore. They now always return an element of type Ciphertext.
- The contexts of BFV and CKKS will not store their checksum anymore, nor will the evaluator check for context consistency of the input and output elements.

### Fixed
- Fixed overflow occurring in the basis extension when small and large moduli are used together.

## [1.0.0] - 2019-08-17
### Added
- First public release.
