# BFV

## Overview

The BFV package is a placeholder for a RNS-accelerated implementation of the Fan-Vercauteren version of Brakerski's (BFV) scale-invariant homomorphic encryption schemes. It enables SIMD modular arithmetic over encrypted vectors or integers.

## Implementation Notes


The implementation of BFV is part of the `bgv` package which implements a unified variant of the BFV and BGV schemes. The only practical difference with the standard BFV is that the plaintext modulus must be coprime with the ciphertext modulus. This is both required for correctness ($T^{-1}\mod Q$ must be defined) and for security reasons (if $T|Q$ then the BGV scheme is not IND-CPA secure anymore). To instantiate the BFV cryptosystem, generate a new BGV evaluator with the optional scale-invariant parameter set to `true`.
Under the hood, setting the scale-invariant flag replaces the following functions of the BGV evaluator object with their scale-invariant alternatives:

  - `Evaluator.Mul`
  - `Evaluator.MulNew`
  - `Evaluator.MulRelin`
  - `Evaluator.MulRelinNew`
  - `Evaluator.Rescale`

For additional information, see the [`README.md`](../bgv/README.md) in the `bgv` package.

## Noise Growth

The only modification proposed in the implementation that could affect the noise is the multiplication, but in theory the noise should behave the same between the two implementations. 

The experiment that follows empirically verifies the above statement.

We instantiated both version of the schemes `BFV_OLD` (textbook BFV) and `BFV_NEW` (wrapper of the generalized BGV) with the following parameters:

```go
ParametersLiteral{
	LogN: 14,
	Q:    []uint64{0x3fffffa8001, 0x1000090001, 0x10000c8001, 0x10000f0001, 0xffff00001},
	P:    []uint64{0x7fffffd8001},
	T: 0xf60001,
}
```

and recorded the average log2 of the standard deviation, minimum and maximum residual noise after 1024 multiplications between two random ciphertexts (without relinearization) encrypted using a public key:

```
 scheme     std       min       max
BFV_OLD | 41.3617 | 26.7891 | 43.4034
BFV_NEW | 40.7618 | 26.2434 | 42.8023
```

We observe that `BFV_NEW` has on average `0.5` bit less noise, but this is due to a fix in the `ring` package where the `ModDown` operation (RNS division by `P`) changing the division from floored to rounded.
