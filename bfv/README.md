# BFV

## Overview

The BFV package provides an RNS-accelerated implementation of the Fan-Vercauteren version of Brakerski's (BFV) scale-invariant homomorphic encryption scheme. It enables SIMD modular arithmetic over encrypted vectors or integers.

## Implementation Notes

The proposed implementation is not standard and is built as a wrapper over the `bgv` package, which implements a unified variant of the BFV and BGV schemes. The only practical difference with the textbook BFV is that the plaintext modulus must be coprime with the ciphertext modulus. This is both required for correctness (T^{-1} mod Q must be defined) and for security reasons (if T divides Q then the BGV scheme is not IND-CPA secure anymore).

For additional information, see the `README.md` in the `bgv` package.
