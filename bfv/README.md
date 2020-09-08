# BFV

The BFV package is an RNS-accelerated implementation of the Fan-Vercauteren version of Brakerski's scale-invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.

## Brief description

This scheme can be used to do arithmetic over &nbsp; ![equation](https://latex.codecogs.com/gif.latex?%5Cmathbb%7BZ%7D_t%5EN).

The plaintext space and the ciphertext space share the same domain

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cmathbb%7BZ%7D_Q%5BX%5D/%28X%5EN%20&plus;%201%29">,
</p>
with <img src="https://latex.codecogs.com/gif.latex?N"> a power of 2.

The batch encoding of this scheme

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cmathbb%7BZ%7D_t%5EN%20%5Cleftrightarrow%20%5Cmathbb%7BZ%7D_Q%5BX%5D/%28X%5EN%20&plus;%201%29">
</p>

maps an array of integers to a polynomial with the property:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?decode%28encode%28m_1%29%20%5Cotimes%20encode%28m_2%29%29%20%3D%20m_1%20%5Codot%20m_2">,
</p>
where  represents &nbsp; ![equation](https://latex.codecogs.com/gif.latex?%24%5Codot%24) &nbsp; a component-wise product,and &nbsp; ![equation](https://latex.codecogs.com/gif.latex?%24%5Cotimes%24) &nbsp; represents a nega-cyclic convolution.

## Security parameters

![equation](https://latex.codecogs.com/gif.latex?N%20%3D%202%5E%7BlogN%7D): the ring dimension, which defines the degree of the cyclotomic polynomial, and the number of coefficients of the plaintext/ciphertext polynomials; it should always be a power of two. This parameter has an impact on both security and performance (security increases with N and performance decreases with N). It should be carefully chosen to suit the intended use of the scheme.

![equation](https://latex.codecogs.com/gif.latex?Q): the ciphertext modulus. In Lattigo, it is chosen to be the product of small coprime moduli ![equation](https://latex.codecogs.com/gif.latex?q_i) that verify ![equation](https://latex.codecogs.com/gif.latex?q_i%20%5Cequiv%201%20%5Cmod%202N) in order to enable both the RNS and NTT representation. The used moduli ![equation](https://latex.codecogs.com/gif.latex?q_i) are chosen to be of size 50 to 60 bits for the best performance. This parameter has an impact on both security and performance (for a fixed ![equation](https://latex.codecogs.com/gif.latex?N), a larger ![equation](https://latex.codecogs.com/gif.latex?Q) implies both lower security and lower performance). It is closely related to ![equation](https://latex.codecogs.com/gif.latex?N) and should be chosen carefully to suit the intended use of the scheme.

![equation](https://latex.codecogs.com/gif.latex?%5Csigma): the variance used for the error polynomials. This parameter is closely tied to the security of the scheme (a larger ![equation](https://latex.codecogs.com/gif.latex?%5Csigma) implies higher security).

## Other parameters

![equation](https://latex.codecogs.com/gif.latex?P): the extended ciphertext modulus. This modulus is used during the multiplication, and it has no impact on the security. It is also defined as the product of small coprime moduli ![equation](https://latex.codecogs.com/gif.latex?p_j) and should be chosen such that ![equation](https://latex.codecogs.com/gif.latex?Q%5Ccdot%20P%20%3E%20Q%5E2) by a small margin (~20 bits). This can be done by using one more small coprime modulus than ![equation](https://latex.codecogs.com/gif.latex?Q).

![equation](https://latex.codecogs.com/gif.latex?t): the plaintext modulus. This parameter defines the maximum value that a plaintext coefficient can take. If a computation leads to a higher value, this value will be reduced modulo the plaintext modulus. It can be initialized with any value, but in order to enable batching, it must be prime and verify ![equation](https://latex.codecogs.com/gif.latex?t%20%5Cequiv%201%20%5Cmod%202N). It has no impact on the security.

## Choosing security parameters

The BFV scheme supports the standard recommended parameters chosen to offer a security of 128 bits for a secret key with uniform ternary distribution ![equation](https://latex.codecogs.com/gif.latex?s%20%5Cin_u%20%5C%7B-1%2C%200%2C%201%5C%7D%5EN), according to the Homomorphic Encryption Standards group (https://homomorphicencryption.org/standard/).  

Each set of parameters is defined by the tuple ![equation](https://latex.codecogs.com/gif.latex?%5C%7Blog_2%28N%29%2C%20log_2%28Q%29%2C%20%5Csigma%5C%7D):

- **{12, 109, 3.2}**
- **{13, 218, 3.2}**
- **{14, 438, 3.2}**
- **{15, 881, 3.2}**

These parameter sets are hard-coded in the file [params.go](https://github.com/ldsec/lattigo/blob/master/bfv/params.go). By default the variance should always be set to 3.2 unless the user is perfectly aware of the security implications of changing this parameter.

Finally, it is worth noting that these security parameters are computed for fully entropic ternary keys (with probability distribution {1/3,1/3,1/3} for values {-1,0,1}). Lattigo uses this fully-entropic key configuration by default. It is possible, though, to generate keys with lower entropy, by modifying their distribution to {(1-p)/2, p, (1-p)/2}, for any p between 0 and 1, which for p>>1/3 can result in low Hamming weight keys (*sparse* keys). *We recall that it has been shown that the security of sparse keys can be considerably lower than that of fully entropic keys, and the BFV security parameters should be re-evaluated if sparse keys are used*.
