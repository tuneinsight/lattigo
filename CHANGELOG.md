# Changelog
All notable changes to this project will be documented in this file. 

## [Unreleased]
### Added
- Bootstrapping for CKKS.
- Modulable CRT decomposition for the key-switching keys.
- Examples for the distributed schemes.
- Network layer implementation of protocols supporting Secure Multiparty Computation (SMC).

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