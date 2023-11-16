# CKKS

The package CKKS is an RNS-accelerated version of the Homomorphic Encryption for Arithmetic of
Approximate Numbers (HEAAN, a.k.a. CKKS) scheme originally proposed by Cheon, Kim, Kim and Song. The
package supports two variants of the scheme: the standard one that encrypts vectors of complex
numbers, and the conjugate-invariant one that encrypts vectors of real numbers, as [proposed by Kim
and Song](https://eprint.iacr.org/2018/952).The `RingType` field of the `Parameter` struct controls
which variant is instantiated:

For `RingType: ring.Standard`, the standard variant of CKKS is used. This requires that all moduli
in the chain are congruent to 1 modulo 2N for N the ring degree. This variant supports packing of up
to N/2 plaintext complex values into a single ciphertext.

For `RingType: ring.ConjugateInvariant`, the conjugate-invariant variant of CKKS is used. This
requires that all moduli in the chain are congruent to 1 modulo 4N for N the ring degree. This
variant supports packing of up to N plaintext real values into a single ciphertext.

## Brief description of the Standard variant

This scheme can be used to do arithmetic over
![equation](https://latex.codecogs.com/gif.latex?%5Cmathbb%7BC%7D%5E%7BN/2%7D). The plaintext space
and the ciphertext space share the same domain

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cmathbb%7BZ%7D_Q%5BX%5D/%28X%5EN%20&plus;%201%29">,
</p>
with <img src="https://latex.codecogs.com/gif.latex?N"> a power of 2.

The batch encoding of this scheme

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cmathbb%7BC%7D%5E%7BN/2%7D%20%5Cleftrightarrow%20Z_Q%5BX%5D/%28X%5EN%20&plus;%201%29">.
</p>

maps an array of complex numbers to a polynomial with the property:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?decode%28encode%28m_1%29%20%5Cotimes%20encode%28m_2%29%29%20%5Capprox%20m_1%20%5Codot%20m_2">,
</p>
where <img src="https://latex.codecogs.com/gif.latex?%5Cotimes"> represents a component-wise product, and <img src="https://latex.codecogs.com/gif.latex?%5Codot"> represents a nega-cyclic convolution.

## Security parameters

![equation](https://latex.codecogs.com/gif.latex?N%20%3D%202%5E%7BlogN%7D): the ring dimension,
which defines the degree of the cyclotomic polynomial, and the number of coefficients of the
plaintext/ciphertext polynomials; it should always be a power of two. This parameter has an impact
on both security and performance (security increases with N and performance decreases with N). It
should be chosen carefully to suit the intended use of the scheme.

![equation](https://latex.codecogs.com/gif.latex?Q): the ciphertext modulus. In Lattigo, it is
chosen to be the product of a chain of small coprime moduli
![equation](https://latex.codecogs.com/gif.latex?q_i) that verify
![equation](https://latex.codecogs.com/gif.latex?q_i%20%5Cequiv%201%20%5Cmod%202N) in order to
enable both the RNS and NTT representation. The used moduli
![equation](https://latex.codecogs.com/gif.latex?q_i) are chosen to be of size 30 to 60 bits for the
best performance. This parameter has an impact on both security and performance (for a fixed
![equation](https://latex.codecogs.com/gif.latex?N), a larger
![equation](https://latex.codecogs.com/gif.latex?Q) implies both lower security and lower
performance). It is closely related to ![equation](https://latex.codecogs.com/gif.latex?N) and
should be carefully chosen to suit the intended use of the scheme.

![equation](https://latex.codecogs.com/gif.latex?%5Csigma): the variance used for the error
polynomials. This parameter is closely tied to the security of the scheme (a larger
![equation](https://latex.codecogs.com/gif.latex?%5Csigma) implies higher security).

## Other parameters

![equation](https://latex.codecogs.com/gif.latex?scale): the plaintext scale. Since complex numbers
are encoded on polynomials with integer coefficients, the original values must be scaled during the
encoding, before being rounded to the nearest integer. The
![equation](https://latex.codecogs.com/gif.latex?scale) parameter is the power of two by which the
values are multiplied during the encoding. It has an impact on the precision of the output and on
the amount of operations a fresh encryption can undergo before overflowing.

## Choosing the right parameters for a given application

There are 3 application-dependent parameters: 
- **LogN**: it determines (a) how many values can be encoded (batched) at once (maximum N/2) in one
  plaintext, and (b) the maximum total modulus bit size (the product of all the moduli) for a given
  security parameter.
- **Modulichain**: it determines how many consecutive scalar and non-scalar multiplications (the
  depth of the arithmetic circuit) can be evaluated before requiring decryption. Since Lattigo
  features an RNS implementation, this parameter requires careful fine-tuning depending on the
  application; i.e., the rescaling procedure can only rescale by one of the RNS modulus at a time,
  whose size has to be chosen when creating the CKKS context. Additionally, the individual size of
  each of the moduli also has an effect on the error introduced during the rescaling, since they
  cannot be powers of 2, so they should be chosen as NTT primes as close as possible to a power of 2
  instead.
- **LogDefaultScale**: it determines the scale of the plaintext, affecting both the precision and the
  maximum allowed depth for a given security parameter.

Configuring parameters for CKKS is very application dependent, requiring a prior analysis of the
circuit to be executed under encryption. The following example illustrates how this parametrization
can be done, showing that it is possible to come up with different parameter sets for a given
circuit, each set having pros and cons.

Let us define the evaluation of an arbitrary smooth function f(x) on an array of ~4000 complex
elements contained in a square of side 2 centered at the complex origin (with values ranging between
-1-1i and 1+1i). We first need to find a good polynomial approximation for the given range. Lattigo
provides an automatic Chebyshev approximation for any given polynomial degree, which can be used for
this purpose (it is also possible to define a different polynomial approximation of lower degree
with an acceptable error).

Let us assume that we find an approximation of degree 5, i.e., *a + bx + cx^3 + dx^5*. This function
can be evaluated with 3 scalar multiplications, 3 additions and 3 non-scalar multiplications,
consuming a total of 4 levels (one for the scalar multiplications and 3 for the non-scalar
multiplications). 

We then need to chose a scale for the plaintext, that will influence both the bit consumption for
the rescaling, and the precision of the computation. If we choose a scale of 2^40, we need to
consume at least 160 bits (4 levels) during the evaluation, and we still need some bits left to
store the final result with an acceptable precision. Let us assume that the output of the
approximation lies always in the square between -20-20i and 20+20i; then, the final modulus must be
at least 5 bits larger than the final scale (to preserve the integer precision).

The following parameters will work for the posed example:

- **LogN** = 13
- **Modulichain** = [45, 40, 40, 40, 40], for a logQ <= 205
- **LogDefaultScale** = 40

But it is also possible to use less levels to have ciphertexts of smaller size and, therefore, a
faster evaluation, at the expense of less precision. This can be achieved by using a scale of 30
bits and squeezing two multiplications in a single level, while pre-computing the last scalar
multiplication already in the plaintext. Instead of evaluating *a + bx + cx^3 + dx^5*, we
pre-multiply the plaintext by d^(1/5) and evaluate *a + b/(d^(1/5))x + c/(d^(3/5)) + x^5*.

The following parameters are enough to evaluate this modified function:

- **LogN** = 13
- **Modulichain** = [35, 60, 60], for a logQ <= 155
- **LogDefaultScale** = 30

To summarize, several parameter sets can be used to evaluate a given function, achieving different
trade-offs for space and time versus precision.

## Choosing secure parameters

The CKKS scheme supports the standard recommended parameters chosen to offer a security of 128 bits
for a secret key with uniform ternary distribution
![equation](https://latex.codecogs.com/gif.latex?s%20%5Cin_u%20%5C%7B-1%2C%200%2C%201%5C%7D%5EN),
according to the Homomorphic Encryption Standards group
(https://homomorphicencryption.org/standard/).

Each set of security parameters is defined by the tuple
![equation](https://latex.codecogs.com/gif.latex?%5C%7Blog_2%28N%29%2C%20log_2%28Q%29%2C%20%5Csigma%5C%7D)
:

- **{12, 109, 3.2}**
- **{13, 218, 3.2}**
- **{14, 438, 3.2}**
- **{15, 881, 3.2}**

As mentioned, setting parameters for CKKS involves not only choosing this tuple, but also defining
the actual moduli chain depending on the application at hand, which is why the provided default
parameter sets have to be fine-tuned, preserving the values of the aforementioned tuples, in order
to maintain the required security level of 128 bits. That is, Lattigo provides a set of default
parameters for CKKS, including example moduli chains, ensuring 128 bit security. The user might want
to choose different values in the moduli chain optimized for a specific application. As long as the
total modulus is equal or below the above values for a given logN, the scheme will still provide a
security of at least 128 bits against the current best known attacks.

Finally, it is worth noting that these security parameters are computed for fully entropic ternary
keys (with probability distribution {1/3,1/3,1/3} for values {-1,0,1}). Lattigo uses this
fully-entropic key configuration by default. It is possible, though, to generate keys with lower
entropy, by modifying their distribution to {(1-p)/2, p, (1-p)/2}, for any p between 0 and 1, which
for p>>1/3 can result in low Hamming weight keys (*sparse* keys). *We recall that it has been shown
that the security of sparse keys can be considerably lower than that of fully entropic keys, and the
CKKS security parameters should be re-evaluated if sparse keys are used*.

## References

1. Homomorphic Encryption for Arithmetic of Approximate Numbers (<https://eprint.iacr.org/2016/421>)
2. A Full RNS Variant of Approximate Homomorphic Encryption (<https://eprint.iacr.org/2018/931>)
3. Approximate Homomorphic Encryption over the Conjugate-invariant Ring (<https://eprint.iacr.org/2018/952>)
4. Approximate Homomorphic Encryption with Reduced Approximation Error (<https://eprint.iacr.org/2020/1118>)
5. On the precision loss in approximate homomorphic encryption (<https://eprint.iacr.org/2022/162>)