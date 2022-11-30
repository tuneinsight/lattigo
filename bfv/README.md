# BFV

The BFV package is an RNS-accelerated implementation of the Fan-Vercauteren version of Brakerski's
scale-invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.

## Brief description

This scheme can be used to do arithmetic over $\mathbb{Z}_t^N$.

The plaintext space and the ciphertext space share the same domain

$$
\mathbb{Z}_Q[X]/(X^N + 1)
$$

with $N$ a power of 2.

The batch encoding of this scheme

$$
\mathbb{Z}_t^N \leftrightarrow \mathbb{Z}_Q[X]/(X^N + 1)
$$

maps an array of integers to a polynomial with the property:

$$
decode(encode(m_1) \otimes encode(m_2)) = m_1 \odot m_2
$$

where $\odot$ represents a component-wise product, and $\otimes$ represents a nega-cyclic convolution.

## Security parameters

$N = 2^{logN}$: the ring dimension,
which defines the degree of the cyclotomic polynomial, and the number of coefficients of the plaintext/ciphertext polynomials; it should always be a power of two. This parameter has an impact on both security and performance (security increases with N and performance decreases with N). It should be carefully chosen to suit the intended use of the scheme.

$Q$: the ciphertext modulus. In Lattigo, it is chosen to be the product of small coprime moduli $q_i$ that verify $q_i \equiv 1 \mod 2N$ in order to
enable both the RNS and NTT representation. The used moduli $q_i$ are chosen to be of size 50 to 60 bits for the best performance. This parameter has an impact on both security and performance (for a fixed $N$, a larger $Q$ implies both lower security and lower performance). It is closely related to $N$ and should be chosen carefully to suit the intended use of the scheme.

$\sigma$: the variance used for the error polynomials. This parameter is closely tied to the security of the scheme (a larger $\sigma$) implies higher security).

## Other parameters

$P$: the extended ciphertext modulus. This modulus
is used during the multiplication, and it has no impact on the security. It is also defined as the
product of small coprime moduli $p_j$ and should be
chosen such that $Q\cdot P > Q^2$ by a small margin (~20 bits). This can be done by using one more small coprime modulus than $Q$.

$t$: the plaintext modulus. This parameter defines
the maximum value that a plaintext coefficient can take. If a computation leads to a higher value,
this value will be reduced modulo the plaintext modulus. It can be initialized with any value, but
in order to enable batching, it must be prime and verify $t \equiv 1 \mod 2N$. It has no impact
on the security.

## Choosing security parameters

The BFV scheme supports the standard recommended parameters chosen to offer a security of 128 bits
for a secret key with uniform ternary distribution
$s \in_u \\{-1, 0, 1\\}^N $, according to the [Homomorphic Encryption Standards group](https://homomorphicencryption.org/standard/).  

Each set of parameters is defined by the tuple $\\{ log_2(N), log_2(Q), \sigma \\}$:

- **{12, 109, 3.2}**
- **{13, 218, 3.2}**
- **{14, 438, 3.2}**
- **{15, 881, 3.2}**

These parameter sets are hard-coded in the file
[params.go](https://github.com/tuneinsight/lattigo/blob/master/bfv/params.go). By default the
variance should always be set to 3.2 unless the user is perfectly aware of the security implications
of changing this parameter.

Finally, it is worth noting that these security parameters are computed for fully entropic ternary keys (with probability distribution $\\{1/3,1/3,1/3\\}$ for values $\\{-1,0,1\\}$). Lattigo uses this fully-entropic key configuration by default. It is possible, though, to generate keys with lower entropy, by modifying their distribution to $\\{(1-p)/2, p, (1-p)/2 \\}$, for any $p$ between 0 and 1, which for $p\gg 1/3$ can result in low Hamming weight keys (*sparse* keys). *We recall that it has been shown that the security of sparse keys can be considerably lower than that of fully entropic keys, and the BFV security parameters should be re-evaluated if sparse keys are used*.
