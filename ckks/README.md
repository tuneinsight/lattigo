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

## Chosing the right parameters for a given application

There are 3 application dependent parameters : 
- **LogN** : it will dictate how many values can be encoded at once (maximum N/2) and the total modulus bit size (the product of all the moduli) for a given security parameter
- **Modulichain** : it will dictate how many scalar and non-scalar multiplications can be evaluated before having to decrypt (the depth of the circuit). This parameter is the one needing the most attention : since we work with an RNS implementation, the rescaling procedure can only rescale by one of the RNS modulus, whose size has to be chosen when creating the ckkscontext. The individual size of the moduli also has an effect on the error introduced during the rescaling since we do not divide by a power of 2, but by an NTT prime close to one
- **Logscale** : it will dictate the scale of the plaintext, affecting both the precision and the maximum allowed depth for a given security parameter

Configuring parameters for CKKS is very application dependent. Using an concrete example is probably the best way to explain how it should be done. The following example will illustrate why it is needed to first analyze the whole circuit to run before choosing, and that it is possible to come up with different parameters for a same circuit, each having pros and cons.

Let's say we want to evaluate an arbitrary smooth function f(x) on an array of ~4000 values whose value can range between -1-1i and 1+1i. We first need to find a good polynomial approximation for the given range. The scheme provides an automatic Chebyshev approximation, but it might not always be the best choice, one might often come up with a different approximation of lower degree with an acceptable error. 

Assume we find an approximation of degree 5, for example *a + bx + cx^3 + dx^5*. This function can be evaluated with 3 scalar multiplications, 3 additions and 3 non-scalar multiplications, for a total of 4 levels consumed (one for the scalar multiplications, 3 for the non-scalar multiplications). 

We then need to chose a scale for the plaintext, that will both influence our bit consumption for the rescaling, but also the precision of the computations. If we chose a scale of 40, we will consume 160 bits of modulus during the evaluation, and we still need some bits left to store the result. If we assume that the output of the approximation will be between -20-20i and 20+20i, then we need to make sure that our final modulus is at least 5 bits larger than the final scale (to retain integer precision).

The following parameters will work for the example : 

- **LogN** = 13
- **Modulichain** = [45, 40, 40, 40, 40], for a logQ <= 205
- **LogScale** = 40

But it is also possible to use less levels to have ciphertext of smaller size and a faster evaluation, at the expense of less precision, by using a scale of 30 and squeezing two multiplications in a single level and pre-computing the last scalar multiplication already on the plaintext. Instead of evaluating *a + bx + cx^3 + dx^5*, we pre-multply the plaintext by d^(1/5) and evaluate *a + b/(d^(1/5))x + c/(d^(3/5)) + x^5*. 

The following parameters are enough to evaluate this modified function :

- **LogN** = 13
- **Modulichain** = [35, 60, 60], for a logQ <= 155
- **LogScale** = 30

To summary, for a same function we can find different set of parameters, for example one using more space and time for a greater precision, one using much less space and time but with less precision. 

## Chosing secure parameters

The CKKS scheme supports official parameters chosen to offer a security of 128 bits for secret with uniform ternary distribution ![equation](https://latex.codecogs.com/gif.latex?s%20%5Cin_u%20%5C%7B-1%2C%200%2C%201%5C%7D%5EN), according to the Homomorphic Encryption Standardization (https://homomorphicencryption.org/standard/).  

The set of security parameters are defined by the three parameters ![equation](https://latex.codecogs.com/gif.latex?%5C%7Blog_2%28N%29%2C%20log_2%28Q%29%2C%20%5Csigma%5C%7D) :

- **{11, 54, 3.2}**
- **{12, 109, 3.2}**
- **{13, 218, 3.2}**
- **{14, 438, 3.2}**
- **{15, 881, 3.2}**

Setting parameters for CKKS is much more application dependant than for BFV. This is why default parameters might not be of much use. We however provide a set of default params for CKKS ensuring 128 bit security. The user might want to use different values in the moduli chain optimized for a specific application. As long as the total modulus is equal or below the above values for a given logN, the scheme will provite a security of at least 128 bits against the current best known attacks.


