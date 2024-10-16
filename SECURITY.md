# Report a Vulnerability
To report a vulnerability please contact us directly using the following email: lattigo@tuneinsight.com.

# Code Review
Lattigo 2.0.0 has been code-reviewed by ELCA in November 2020 and, within the allocated time for the code review, no critical or high-risk issues were found.

# Security of Approximate Homomorphic Encryption
Homomorphic encryption schemes are by definition malleable, and are therefore not secure against chosen ciphertext attacks (CCA security). They can be though secure against chosen plaintext attacks (CPA security).

Classified as an _approximate decryption_ scheme, the CKKS scheme is secure as long as the plaintext result of a decryption is only revealed to entities with knowledge of the secret-key. This is because, given a ciphertext $(-as + m + e, a)$, the decryption outputs a plaintext $m+e$. [Li and Micciancio](https://eprint.iacr.org/2020/1533) show that using this plaintext, it is possible to recover the secret-key with $((-as + m + e) - (m + e)) \cdot a^{-1} = asa^{-1} = s$ (the probability of $a$ being invertible is overwhelming, and if $a$ is not invertible, only a few more samples are required).

This attack demonstrates that, when using an approximate homomorphic encryption scheme, the usual CPA security may not sufficient depending on the application setting. Many applications do not require to share the result with external parties and are not affected by this attack, but the ones that do must take the appropriate steps to ensure that no key-dependent information is leaked. A homomorphic encryption scheme that provides such functionality and that can be secure when releasing decrypted plaintext to external parties is defined to be CPA<sup>D</sup> secure. The corresponding indistinguishability notion (IND-CPA<sup>D</sup>) is defined as "indistinguishability under chosen plaintext attacks with decryption oracles."

## IND-CPA<sup>D</sup> Security for Approximate Homomorphic Encryption
Lattigo implements tools to mitigate _Li and Micciancio_'s attack. In particular, the decoding step of CKKS (and its real-number variant R-CKKS) allows the user to specify the desired fixed-point bit-precision.

Let $\epsilon$ be the scheme error after the decoding step. We compute the bit precision of the output as $\log_{2}(1/\epsilon)$.

If at any point of an application, decrypted values have to be shared with external parties, then the user must ensure that each shared plaintext is first _sanitized_ before being shared. To do so, the user must use the $\textsf{DecodePublic}$ method instead of the usual $\textsf{Decode}$. $\textsf{DecodePublic}$ takes as additional input the desired $\log_{2}(1/\epsilon)$-bit precision and rounds the value by evaluating $y = \lfloor x / \epsilon \rceil \cdot \epsilon$.

Estimating $\text{Pr}[\epsilon < x] \leq 2^{-s}$ of the circuit must be done carefully and we suggest the following process to do so:
 1. Given a security parameter $\lambda$ and a circuit $C$ that takes as inputs length-$n$ vectors $\omega$ following a distribution $\chi$, select the appropriate parameters enabling the homomorphic evaluation of $C(\omega)$, denoted by $H(C(\omega))$, which includes the encoding, encryption, evaluation, decryption and decoding.
 2. Sample input vectors $\omega$ from the distribution $\chi$ and record $\epsilon = C(\omega) - H(C(\omega))$ for each slots. The user should make sure that the underlying circuit computed by $H(C(\cdot))$ is identical to $C(\cdot)$; i.e., if the homomorphic implementation $H(C(\cdot))$ uses polynomial approximations, then $C(\cdot)$ should use them too, instead of using the original exact function. Repeat until enough data points are collected to construct a CDF of $\textsf{Pr}[\epsilon > x]$.
 3. Use the CDF to select the value $\text{E}[\epsilon]$ such that any given slot will fail with probability $2^{-\varepsilon}$ (where $\varepsilon$ is a user-defined security parameter) to reach $\log_{2}(1/\epsilon)$ bits of precision. 
 4. Use the encoder method $\textsf{DecodePublic}$ with the parameter $\log_{2}(1/\epsilon)$ to decode plaintexts that will be published.

Note that, for composability with differential privacy, the variance of the error introduced by the rounding is $\text{Var}[x - \lfloor x \cdot \epsilon \rceil / \epsilon] = \tfrac{\epsilon^2}{12}$ and therefore $\text{Var}[x -  \lfloor x/(\sigma\sqrt{12})\rceil\cdot(\sigma\sqrt{12})] = \sigma^2$.

[Bossuat et al.](https://eprint.iacr.org/2024/853) recent research paper provides tight bounds on noise to optimize the rounding process, minimizing loss in both precision and efficiency.
In Lattigo, we are planning to implement a detailed noise analysis for all basic operations, including bootstrapping, based on this work. To support this, we will provide a noise estimator tool that combines the noise bounds for individual operations, allowing for accurate estimates even for complex circuits. 

# Security of Exact Homomorphic Encryption
In recent papers [Checri et al.](https://eprint.iacr.org/2024/116) and [Cheon et al.](https://eprint.iacr.org/2024/127), the authors revealed new passive key-recovery attacks targeting also the exact FHE cryptosystems, including BFV, BGV, and TFHE. They exploit imperfect correctness and show that BFV, BGV and TFHE are not protected against IND-CPA<sup>D</sup> attackers.

## IND-CPA<sup>D</sup> Security for Exact Homomorphic Encryption
Achieving IND-CPA<sup>D</sup> security for the exact homomorphic encryption schemes requires near-perfect correctness, meaning decryption failures must be exceptionally rare, with a probability lower than $2^{−\lambda}$, where $\lambda$ is a user-defined security parameter. Such failures should be so unlikely that finding one is computationally infeasible.
For exact schemes like BFV and BGV, implemented in Lattigo, near-perfect correctness can be maintained by adjusting scheme parameters to bound decryption noise, though this comes at the cost of performance. The scheme must also control noise growth by limiting the number and type of operations performed at each computation level. To support this, as before we are planning to provide a noise estimator tool.

# Security of Multiparty/Threshold Homomorphic Encryption
Multiparty or Threshold Fully Homomorphic Encryption involve secret-sharing the encryption key among multiple users, requiring their collaboration to decrypt data. The scheme maintains security even if a small number of users are compromised by an attacker. However, since all users receive the decrypted value, and some may be corrupted, the attacker can gain knowledge of the computation's result. This scenario highlights the necessity of IND-CPA<sup>D</sup> security within threshold-FHE. However, it is important to note that IND-CPA<sup>D</sup> security alone is insufficient, as the attacker not only learns the decrypted value but also gains access to shares of the secret key.
All existing solutions for achieving Multiparty/Threshold Fully Homomorphic Encryption necessitate exponential flooding noise at the end of the evaluation process.

It is important to clarify that, for many applications, IND-CPA is sufficient. IND-CPA<sup>D</sup>  model does not always align well with practical use cases. This model fails to address a broad range of real-world adversarial scenarios—for example, an attacker capable of producing or computing ciphertexts is unlikely to be limited to the constraints defined by the IND-CPA<sup>D</sup> model. Furthermore, it introduces unnecessary performance drawbacks for FHE. However, depending on the specific application and threat model, additional security guarantees may be necessary. In such cases, these enhanced guarantees can often be achieved through supplementary application-specific countermeasures or protocol modifications.

# Recommendation for applicative countermeasures
1. FHE ciphertexts are inherently malleable, and this malleability, combined with vulnerabilities such as circular security and decision-to-search attacks, can lead to key-recovery attacks. As a foundational principle, it’s crucial that FHE ciphertexts are transmitted only through private and authenticated channels, encapsulated within traditional cryptographic methods.
2. Most IND-CPA<sup>D</sup> attacks require hundreds of thousands of queries to the evaluation and decryption oracles. By employing ephemeral keys or key rotations, we can limit the number of available plaintext-ciphertext pairs, effectively reducing the attack surface and mitigating these types of attacks.
3. A zero-knowledge proof can be employed to verify both the correctness of the ciphertext and the prover's knowledge of the corresponding plaintext, without revealing any information about the plaintext itself.
4. To ensure the accuracy of the public results, it is recommended that independent parties independently replicate the computation.
5. Circuit Privacy: Ensuring that the output of an FHE computation does not leak any secret information from the evaluator.
6. other physical limitations: firewall, rate control, enclaves


