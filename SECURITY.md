# Report a Vulnerability
To report a vulnerability please contact us directly using the following email: lattigo@tuneinsight.com.

# Code Review
Lattigo 2.0.0 has been code-reviewed by ELCA in November 2020 and, within the allocated time for the code review, no critical or high-risk issues were found.

# Security of Approximate Homomorphic Encryption
Homomorphic encryption schemes are by definition malleable, and are therefore not secure against chosen ciphertext attacks (CCA security). They can be though secure against chosen plaintext attacks (CPA security).

Classified as an _approximate decryption_ scheme, the CKKS scheme is secure as long as the plaintext result of a decryption is only revealed to entities with knowledge of the secret-key. This is because, given a ciphertext $(-as + m + e, a)$, the decryption outputs a plaintext $m+e$. [Li and Micciancio](https://eprint.iacr.org/2020/1533) show that using this plaintext, it is possible to recover the secret-key with $((-as + m + e) - (m + e)) \cdot a^{-1} = asa^{-1} = s$ (the probability of $a$ being invertible is overwhelming, and if $a$ is not invertible, only a few more samples are required).

This attack demonstrates that, when using an approximate homomorphic encryption scheme, the usual CPA security may not sufficient depending on the application setting. Many applications do not require to share the result with external parties and are not affected by this attack, but the ones that do must take the appropriate steps to ensure that no key-dependent information is leaked. A homomorphic encryption scheme that provides such functionality and that can be secure when releasing decrypted plaintext to external parties is defined to be CPA<sup>D</sup> secure. The corresponding indistinguishability notion (IND-CPA<sup>D</sup>) is defined as "indistinguishability under chosen plaintext attacks with decryption oracles."

# CPA<sup>D</sup> Security for Approximate Homomorphic Encryption
Lattigo implements tools to mitigate _Li and Micciancio_'s attack. In particular, the decoding step of CKKS (and its real-number variant R-CKKS) allows the user to specify the desired fixed-point bit-precision.

Let $\epsilon$ be the scheme error after the decoding step. We compute the bit precision of the output as $\log_{2}(1/\epsilon)$.

If at any point of an application, decrypted values have to be shared with external parties, then the user must ensure that each shared plaintext is first _sanitized_ before being shared. To do so, the user must use the $\textsf{DecodePublic}$ method instead of the usual $\textsf{Decode}$. $\textsf{DecodePublic}$ takes as additional input the desired $\log_{2}(1/\epsilon)$-bit precision and rounds the value by evaluating $y = \lfloor x / \epsilon \rceil \cdot \epsilon$.

Estimating $\text{Pr}[\epsilon < x] \leq 2^{-s}$ of the circuit must be done carefully and we suggest the following process to do so:
 1. Given a security parameter $\lambda$ and a circuit $C$ that takes as inputs length-$n$ vectors $\omega$ following a distribution $\chi$, select the appropriate parameters enabling the homomorphic evaluation of $C(\omega)$, denoted by $H(C(\omega))$, which includes the encoding, encryption, evaluation, decryption and decoding.
 2. Sample input vectors $\omega$ from the distribution $\chi$ and record $\epsilon = C(\omega) - H(C(\omega))$ for each slots. The user should make sure that the underlying circuit computed by $H(C(\cdot))$ is identical to $C(\cdot)$; i.e., if the homomorphic implementation $H(C(\cdot))$ uses polynomial approximations, then $C(\cdot)$ should use them too, instead of using the original exact function. Repeat until enough data points are collected to construct a CDF of $\textsf{Pr}[\epsilon > x]$.
 3. Use the CDF to select the value $\text{E}[\epsilon]$ such that any given slot will fail with probability $2^{-\varepsilon}$ (where $\varepsilon$ is a user-defined security parameter) to reach $\log_{2}(1/\epsilon)$ bits of precision. 
 4. Use the encoder method $\textsf{DecodePublic}$ with the parameter $\log_{2}(1/\epsilon)$ to decode plaintexts that will be published.

Note that, for composability with differential privacy, the variance of the error introduced by the rounding is $\text{Var}[x - \lfloor x \cdot \epsilon \rceil / \epsilon] = \tfrac{\epsilon^2}{12}$ and therefore $\text{Var}[x -  \lfloor x/(\sigma\sqrt{12})\rceil\cdot(\sigma\sqrt{12})] = \sigma^2$.
