# CKKS

The package CKKS is a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers (HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.

## Brief description

This scheme can be used to do arithmetic over ![equation](https://latex.codecogs.com/gif.latex?%5Cmathbb%7BC%7D%5E%7BN/2%7D). The plaintext space and the ciphertetext space share the same domain :

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cmathbb%7BZ%7D_Q%5BX%5D/%28X%5EN%20&plus;%201%29">
</p>

The encoding of this scheme :

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cmathbb%7BC%7D%5E%7BN/2%7D%20%5Cleftrightarrow%20Z_Q%5BX%5D/%28X%5EN%20&plus;%201%29">
</p>

maps an array of complex number to a polynomial with the property :

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?decode%28encode%28m_1%29%20%5Cotimes%20encode%28m_2%29%29%20%5Capprox%20m_1%20%5Codot%20m_2">
</p>


## Security parameters

![equation](https://latex.codecogs.com/gif.latex?N%20%3D%202%5E%7BlogN%7D) : the ring dimension, it should always be a power of two. This parameters has an impact on both security and performances. It should be chosen carefuly to suit the intended use of the scheme.

![equation](https://latex.codecogs.com/gif.latex?Q) : the ciphertext modulus. In this implementation it is chosen to be the a chain of small coprime moduli ![equation](https://latex.codecogs.com/gif.latex?q_i) of 40 to 60 bits verifying ![equation](https://latex.codecogs.com/gif.latex?q_i%20%5Cequiv%201%20%5Cmod%202N) to enable the RNS and NTT representation. This parameter has an impact on both security and performances. It is closely related to ![equation](https://latex.codecogs.com/gif.latex?N) and should be chosen carefuly to suit the intended use of the scheme.

![equation](https://latex.codecogs.com/gif.latex?%5Csigma) : the variance used for the error polynomials. This parameter is closely tied to the security of the scheme.

## Other parameters

![equation](https://latex.codecogs.com/gif.latex?scale) : the plaintext scale. Since complex number are encoded on polynomial with integer coefficients, the values must be scaled during the encoding before being rounded to the nearest integer. This parameter is the power of two by which the values are multiplied during the encoding. It has an impact on the precision of the output, as well as the performances, since a bigger scale allows for less computations.

## Chosing security parameters

The CKKS scheme supports official parameters chosen to offer a security of 128 bits for secret with uniform ternary distribution ![equation](https://latex.codecogs.com/gif.latex?s%20%5Cin_u%20%5C%7B-1%2C%200%2C%201%5C%7D%5EN), according to the Homomorphic Encryption Standardization (https://homomorphicencryption.org/standard/).  

The set of security parameters are defined by the three parameters ![equation](https://latex.codecogs.com/gif.latex?%5C%7Blog_2%28N%29%2C%20log_2%28Q%29%2C%20%5Csigma%5C%7D) :

- {11, 56, 3.2}
- {12, 110, 3.2}
- {13, 219, 3.2}
- {14, 441, 3.2}
- {15, 885, 3.2}

Contrary to BFV we have chosen to not hardcode those parameters. Instead the user must chose logN, the size of each moduli and the number of moduli and make sure that the product of all the moduli is compliant with the security parameters (which are the same as for BFV). For example, if logN = 14 is chosen, the user must make sure that total bisize of its modulus (the product of all the moduli) does not exceed 441. This can be estimated by computing number of moduli times the bitsize of each modulus.
