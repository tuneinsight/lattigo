# Lattigo: lattice-based cryptographic library in Go

<p align="center">
	<img src="https://i.ibb.co/5vXt5mj/logo.png" />
</p>

_The Lattigo library unleashes the potential of lattice-based cryptography in secure multiparty computation for modern software stacks._

[![Build Status](https://travis-ci.com/ldsec/lattigo.svg?token=kz1BaknyyJcURGZurf6m&branch=master)](https://travis-ci.com/ldsec/lattigo)

Lattigo is a Go package implementing lattice-based cryptographic primitives.
The library features:
- A pure Go implementation bringing code-simplicity and easy builds.
- A public interface for an efficient multiprecision polynomial arithmetic layer.
- Comparable performance to state-of-the-art C++ libraries.

Lattigo aims at enabling fast prototyping of secure-multiparty computation solutions based on distributed homomorphic cryptosystems, by harnessing Go's natural concurrency model.

## Library overview

The library comprises the following sub-packages:

- `lattigo/ring`: RNS-accelerated modular arithmetic operations for polynomials, including: RNS basis extension; RNS rescaling;  number theoretic transform (NTT); uniform, Gaussian and ternary sampling.

- `lattigo/bfv`: RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
	
- `lattigo/ckks`: RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers (HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.

- `lattigo/dbfv` and `lattigo/dckks`: Distributed (or threshold) versions of the BFV and CKKS schemes that enable secure multiparty computation solutions with secret-shared secret keys.

- `lattigo/examples`: Executable Go programs demonstrating the usage of the Lattigo library.
                      Note that each subpackage includes test files that further demonstrate the usage of Lattigo primitives.

- `lattigo/utils`: Supporting structures and functions.

## Versions and Roadmap

The first version of Lattigo, (v1.x.x) is meant to be a working prototype.
Hence, there will be backward-incompatible changes within this version. 

See CHANGELOG.md for the current and past versions.

### Upcoming features

- Bootstrapping for CKKS.
- Modulable CRT decomposition for the key-switching keys.
- Examples for the distributed schemes.
- Network layer implementation of protocols supporting Secure Multiparty Computation (SMC).


## Disclaimer

This library is still at an experimental stage and should be used for research purposes only.

## License

Lattigo is licensed under the Apache 2.0 License.

## Contact

If you want to contribute to Lattigo or you have any suggestion, do not hesitate to contact us at [lattigo@listes.epfl.ch](mailto:lattigo@listes.epfl.ch).

## Citing

Please use the following BibTex entry for citing Lattigo:

    @misc{lattigo,
	    title = {Lattigo 1.1.0},
	    howpublished = {Online: \url{http://github.com/ldsec/lattigo}},
	    month = oct,
	    year = 2019,
	    note = {EPFL-LDS}
    }
    


## References

1. Somewhat Practical Fully Homomorphic Encryption (<https://eprint.iacr.org/2012/144>)
1. Computing accross Trust Boundaries using Distributed Homomorphic Cryptography (<https://eprint.iacr.org/2019/961>)
1. A Full RNS Variant of FV Like Somewhat Homomorphic Encryption Schemes (<https://eprint.iacr.org/2016/510>)
1. An Improved RNS Variant of the BFV Homomorphic Encryption Scheme (<https://eprint.iacr.org/2018/117>)
1. Homomorphic Encryption for Arithmetic of Approximate Numbers (<https://eprint.iacr.org/2016/421>)
1. A Full RNS Variant of Approximate Homomorphic Encryption (<https://eprint.iacr.org/2018/931>)
1. Improved Bootstrapping for Approximate Homomorphic Encryption (<https://eprint.iacr.org/2018/1043>)
1. Post-quantum key exchange - a new hope (<https://eprint.iacr.org/2015/1092>)
1. Faster arithmetic for number-theoretic transforms (<https://arxiv.org/abs/1205.2926>)
1. Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography (<https://eprint.iacr.org/2016/504>)
1. Gaussian sampling in lattice-based cryptography (<https://tel.archives-ouvertes.fr/tel-01245066v2>)
