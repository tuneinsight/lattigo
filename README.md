# Lattigo: lattice-based multiparty homomorphic encryption library in Go

<p align="center">
	<img src="logo.png" />
</p>

![Go tests](https://github.com/tuneinsight/lattigo/actions/workflows/ci.yml/badge.svg)

Lattigo is a Go module that implements Ring-Learning-With-Errors-based homomorphic-encryption
primitives and Multiparty-Homomorphic-Encryption-based secure protocols. The library features:
- An implementation of the full-RNS BFV and CKKS schemes and their respective multiparty versions.
- Comparable performance to state-of-the-art C++ libraries.
- Dense-key and sparse-key efficient and high-precision bootstrapping procedures for full-RNS CKKS.
- A pure Go implementation that enables cross-platform builds, including WASM compilation for
  browser clients.

Lattigo is meant to support HE in distributed systems and microservices architectures, for which Go
is a common choice thanks to its natural concurrency model and portability.

## Library overview

The library exposes the following packages:

- `lattigo/ring`: Modular arithmetic operations for polynomials in the RNS basis, including: RNS
  basis extension; RNS rescaling; number theoretic transform (NTT); uniform, Gaussian and ternary
  sampling.

- `lattigo/bfv`: The Full-RNS variant of the Brakerski-Fan-Vercauteren scale-invariant homomorphic
  encryption scheme. It provides modular arithmetic over the integers.
	
- `lattigo/ckks`: The Full-RNS Homomorphic Encryption for Arithmetic for Approximate Numbers (HEAAN,
  a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers (in its classic
  variant) and over the real numbers (in its conjugate-invariant variant).

- `lattigo/dbfv` and `lattigo/dckks`: Multiparty (a.k.a. distributed or threshold) versions of the
  BFV and CKKS schemes that enable secure multiparty computation solutions with secret-shared secret
  keys.

- `lattigo/rlwe` and `lattigo/drlwe`: common base for generic RLWE-based multiparty homomorphic
  encryption. It is imported by the `lattigo/bfv` and `lattigo/ckks` packages.

- `lattigo/examples`: Executable Go programs that demonstrate the use of the Lattigo library. Each
                      subpackage includes test files that further demonstrate the use of Lattigo
                      primitives.

- `lattigo/utils`: Supporting structures and functions.

## Versions and Roadmap

The Lattigo library was originally exclusively developed by the EPFL Laboratory for Data Security
until its version 2.4.0.

Starting with the release of version 3.0.0, Lattigo is maintained and supported by [Tune Insight
SA](https://tuneinsight.com).

Also starting with from version 3.0.0, the module name has changed to
github.com/tuneinsight/lattigo/v3, and the official repository has been moved to
https://github.com/tuneinsight/lattigo. This has the following implications for modules that depend
on Lattigo:
- Modules that require github.com/ldsec/lattigo/v2 will still build correctly.
- To upgrade to a version >= 3.0.0, modules have to require `github.com/tuneinsight/lattigo/v3/`,
  for example by changing the imports to `github.com/tuneinsight/lattigo/v3/[package]` and by
  running `go mod tidy`.


The current version of Lattigo, (v3.x.x) is fast-evolving and in constant development. Consequently,
there will still be backward-incompatible changes within this major version, in addition to many bug
fixes and new features. Hence, we encourage all Lattigo users to update to the latest Lattigo
version.
 

See CHANGELOG.md for the current and past versions.

## License

Lattigo is licensed under the Apache 2.0 License. See LICENSE.

## Contact

If you want to contribute to Lattigo or you have any suggestion, do not hesitate to contact us at
[lattigo@tuneinsight.com](mailto:lattigo@tuneinsight.com).

## Citing

Please use the following BibTex entry for citing Lattigo:

    @misc{lattigo,
	    title = {Lattigo v3.0.0},
	    howpublished = {Online: \url{https://github.com/tuneinsight/lattigo}},
	    month = Feb,
	    year = 2022,
	    note = {EPFL-LDS, Tune Insight SA}
    }
    

## References

1. Efficient Bootstrapping for ApproximateHomomorphic Encryption with Non-Sparse Keys
   (<https://eprint.iacr.org/2020/1203>)
1. Somewhat Practical Fully Homomorphic Encryption (<https://eprint.iacr.org/2012/144>)
1. Multiparty Homomorphic Encryption: From Theory to Practice (<https://eprint.iacr.org/2020/304>)
1. A Full RNS Variant of FV Like Somewhat Homomorphic Encryption Schemes
   (<https://eprint.iacr.org/2016/510>)
1. An Improved RNS Variant of the BFV Homomorphic Encryption Scheme
   (<https://eprint.iacr.org/2018/117>)
1. Homomorphic Encryption for Arithmetic of Approximate Numbers (<https://eprint.iacr.org/2016/421>)
1. A Full RNS Variant of Approximate Homomorphic Encryption (<https://eprint.iacr.org/2018/931>)
1. Improved Bootstrapping for Approximate Homomorphic Encryption
   (<https://eprint.iacr.org/2018/1043>)
1. Better Bootstrapping for Approximate Homomorphic Encryption (<https://epring.iacr.org/2019/688>)
1. Post-quantum key exchange - a new hope (<https://eprint.iacr.org/2015/1092>)
1. Faster arithmetic for number-theoretic transforms (<https://arxiv.org/abs/1205.2926>)
1. Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography
   (<https://eprint.iacr.org/2016/504>)
1. Gaussian sampling in lattice-based cryptography
   (<https://tel.archives-ouvertes.fr/tel-01245066v2>)

The Lattigo logo is a lattice-based version of the original Golang mascot by [Renee
French](http://reneefrench.blogspot.com/).
