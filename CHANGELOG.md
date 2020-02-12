# Changelog
All notable changes to this project will be documented in this file. 

## [Unreleased]
### Added
- Bootstrapping for CKKS.
- Network layer implementation of protocols supporting Secure Multiparty Computation (SMC).

## [1.3.1] - 2020-02-12
### Added
- BFV/CKKS : added API for encrypting using a CRP (common reference polynomial).
- BFV/CKKS : added API for encrypting faster (encrypts zero directly in Q instead of QP and does not need to divide by P).
- BFV : added hoisted rotations.
- BFV/CKKS : added tests for hoisted rotations.
- RING : added benchmarks for a NTT using purely Barrett reduction for comparison purposes.
### Changed :
- BFV/CKKS : changed the switching keys from (-as1 + (s0-s1) + e, a) to (-as1 + s0 + e, a).
### Fixes
- BFV/CKKS : Fixed EncryptFromSK that was not correctly wiping the memory pool before using it, which lead to back encryptions.
- BFV : Fixed an index out of bound error that would happen during the multiplication if #QMul > #Qi.
- CKKS : removed some redundant operations in the hoisted rotations.

## [1.3.0] - 2019-12-20
### Added
- All schemes : new switching-keys and key-switching algorithm based on the concept presented in https://eprint.iacr.org/2019/688.pdf.
- All schemes : new marshaling interface for all structures.
- BFV/CKKS : new Parameters structs and API enabling a better customization and fine tuning for specific applications.
- CKKS : new API for hoisted rotations, which is faster than sequential rotations.
- DBFV/DCKKS : added collective refresh of a ciphertext (decentralized bootstrapping).
- RING : added Ziggurat sampling, available from the context.
- RING : enabled dense and sparse ternary polynomials sampling directly from the context.
- RING : new API enabling "level" wise polynomial arithmetic.
- RING : new API for modulus switching with flooring and rounding.
- UTILS : utils now regroups all the utility methods which were previously duplicated among packages.
### Removed
- BFV/CKKS/DBFV/DCKKS : removed their respective context. Ring context remains public.
- All schemes : removed key-switching with bit decomposition. This option will however be re-introduced at a later stage since applications using small parameters can suffer from this change.
- BFV/CKKS/RING : removed redudant/irrelevant tests and benchmarks.
- BFV : removed context QP as it is not any more used in the multiplication.
- BFV : removed int encoder, now only batch encoding is supported.
- CKKS : modulus switching is now located in Ring.
- RING : removed the algorithms that needed Float128 during the BFV multiplication.
- RING : removed most wrapping methods for bigInt, which are now replaced by the native math/big package.
- RING : removed ternary sampler, which is now part of the context.
### Changed
- All schemes : Encryptor, Decryptor, Encoder, Evaluator, KeyGenerator are now interface types.
- All schemes : Improved Godoc and error strings.
- ALl schemes : greatly reduced the number of methods that could return an error.
- All schemes : new tests and benchmarks with fully supported regex.
- All schemes : coefficient wise arithmetic using double slices is now substentially faster.
- BFV/CKKS : changed the name of the underlying ring contexts. Q now represents the ciphertext modulus (with QMul being the extended ciphertext modulus for BFV) and QP represents modulus of the keys (P being the special primes used during the new key-switching).
- BFV/CKKS/DBFV/DCKKS : structures are now created using the parameters instead of the context.
- BFV : quantization during multiplication doesn't use Float128 any more, resulting in a substential speed improvement.
- BFV : BatchEncoder has been renamed Encoder.
- CKKS : the scale is now stored as a float64 instead of a power of 2.
- CKKS : rounding is applied instead of flooring when a real value is converted to an integer value. This change affects the rescaling and the encoding.
- CKKS : previously needed one ring context per level, now only uses one context for all levels.
- CKKS : new baby-step giant-step algorithm for evaluating polynomials in standard and Chebyshev basis.
- CKKS : reduced the number of NTT needed during the encryption.
- CKKS : API for MultConst is now MultByConst.
- BFV/CKKS : new API for the rotation-keys generation.
- DBFV/DCKKS : complete revamp of the API and interfaces enabling a much easier integration into larger systems.
- DBFV/DCKKS : improved PCKS and CKS using the concept of the new key-switching technique which enables to reduces the added noise.
- DCKKS : all protocols work for ciphertexts at any levels.
- RING : faster MulScalarBigint (now similar to MulScalar).
- UTILS : PRNG must be keyed to be forward secure.
### Fixes
- All packages : typos, godoc and golint.
- CKKS : ciphertext rotation now correctly sets the scale of the output ciphertext.
- DBFV/DCKKS : correctness is now ensured when the same protocol instance is used to generate multiples shares.

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
- BFV and CKKS: Operand interface.
- Ring: New method MultByVector.
- Ring (API change): new Ternary Sampler, which allows to specify the key distribution {-1, 0, 1} -> [(1-p)/2, p, (1-p)/2]; faster than the previous implementation.
- GoDoc for BFV, CKKS, Ring, DBFV, and DCKKS.
- README for BFV and CKKS.
- Code cleaning for all packages.
- Minor optimizations for all packages.

### Changed
- Updated README.md.
- BFV (API change): bfvcontext now only accepts as input a struct of the type Parameters, similar to the one used for DefaultParams.
- BFV (API change): removed bfvcontext from ciphertexts and plaintexts.
- BFV (API change): encryption can now also be done with the secret-key.
- BFV: The value logQ in the bfvcontext now stores the bit size of the product of the ciphertext's moduli.
- BFV: The information printed by the tests now better conveys the parameters.
- BFV: Updated and optimized default parameters.
- BFV: Updated example with secure parameters.
- CKKS (API change): ckkscontext now only accepts as input a struct of the type Parameters, similar to the one used for DefaultParams.
- CKKS (API change): ckkscontext now requires as input the moduli chain (in bit size), instead of a generic logQ and levels. This allows more fine-grained control on the rescaling and levels.
- CKKS (API change): removed ckkscontext from ciphertexts and plaintexts.
- CKKS: Updated the value logQ in the ckkscontext, which now stores the bit size of the product of the ciphertext's moduli.
- CKKS: Updated the information printed by the tests to more reflect the parameters.
- BFV and CKKS: greatly simplified the code related to the rotations.
- BFV and CKKS: reduced the number of instances where an error could be returned and updated informations of the returning errors.
- Ring (API change): changed the API of the validation methode of the context to better reflect its purpose.
- BFV/CKKS/Ring : the copy function will now copy the input poly on the target poly (previously the target was copied on the input).
- Updates on all packages and tests to comply with the API changes in BFV and CKKS.

### Removed
- The evaluator of both BFV and CKKS cannot anymore operate on two plaintexts. They now always returns an element of type Ciphertext.
- The contexts of BFV and CKKS will not anymore store their checksum, nor will the evaluator check for context consistency of the input and output elements.

### Fixed
- Fixed overflow occurring in the basis extension when small and large moduli are used together.

## [1.0.0] - 2019-08-17
### Added
- First public release.
