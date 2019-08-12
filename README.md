# LattiGo: Lattice Cryptography Library in Golang
[![Build Status](https://travis-ci.com/lca1/lattigo.svg?token=kz1BaknyyJcURGZurf6m&branch=dev)](https://travis-ci.com/lca1/lattigo)

[![Build Status](https://travis-ci.com/lca1/lattigo.svg?token=kz1BaknyyJcURGZurf6m&branch=dev)](https://travis-ci.com/lca1/lattigo)

This package provides a toolbox of lattice-based cryptographic primitives for Go. The library is still at an experimental stage and should be used for research purposes only.

The __LattiGo__ subpackages from the lowest to the highest abstraction level and their provided functionalities are as follows:



- `Ring`: CRT accelerated modular arithmetic operations for polynomials, including CRT basis extension, CRT rescaling,  Number Theoretic Transformation (NTT), uniform, Gaussian and Ternary sampling.
	- Post-quantum key exchange - a new hope (<https://eprint.iacr.org/2015/1092>)
	- Faster arithmetic for number-theoretic transforms (<https://arxiv.org/abs/1205.2926>)
	- Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography (<https://eprint.iacr.org/2016/504>)
	- Gaussian sampling in lattice-based cryptography (<https://tel.archives-ouvertes.fr/tel-01245066v2>)

- `BFV`: CRT accellerated Fan-Vercauteren variant of Brakerski's scale invariant homomorphic encryption scheme. Provides modular arithmetic over the integers.
	- Somewhat Practical Fully Homomorphic Encryption (<https://eprint.iacr.org/2012/144>).
	- A Full RNS Variant of FV Like Somewhat Homomorphic Encryption Schemes (<https://eprint.iacr.org/2016/510>)
	- An Improved RNS Variant of the BFV Homomorphic Encryption Scheme (<https://eprint.iacr.org/2018/117>)

- `CKKS`: CRT accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers scheme. Provides approximate arithmetic over the complex numbers.
	- Homomorphic Encryption for Arithmetic of Approximate Numbers (<https://eprint.iacr.org/2016/421>)
	- A Full RNS Variant of Approximate Homomorphic Encryption (<https://eprint.iacr.org/2018/931>)
	- Improved Bootstrapping for Approximate Homomorphic Encryption (<https://eprint.iacr.org/2018/1043>)

- `Distributed BFV/CKKS`: Distributed version of the BFV and CKKS. Provides protocols for the generation of shared private and public keys, and switching keys (evaluation keys, switching keys and rotation keys), as well as a protocol to operate the key-switching on a shared ciphertext.  


## Examples

In the subfolder "examples" you can find runnable go files giving examples on the usage of the library.
In each subpackage you can find additional test files documenting further usage approaches.

## License

__LattiGo__ is licenced under the Apache License version 2.0.
