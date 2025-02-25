# Report a Vulnerability
To report a vulnerability please contact us directly using the following email: lattigo@tuneinsight.com.

# Code Review
Lattigo 2.0.0 was code-reviewed by ELCA in November 2020 and, within the allocated time for the code review, no critical or high-risk issues were found.

# IND-CPA Security

In IND-CPA security, the adversary has access only to an encryption oracle. This means that to attempt breaking the cryptographic scheme, the adversary can select multiple plaintexts and observe the corresponding well-formed ciphertexts produced by the encryption function. IND-CPA serves as the baseline security standard for cryptosystems with provable security.  
Notably, all widely used FHE schemes —such as BFV, BGV, CKKS (implemented in the current version of Lattigo), and TFHE— provide only CPA security, which refers to protection against passive adversaries. In Lattigo, the default parameters provide at least 128-bits of security (according to the [Lattice estimator](https://github.com/malb/lattice-estimator)). 


# IND-CPA-D Security 
IND-CPA-D security strengthens IND-CPA by granting the attacker access to a decryption oracle for ciphertexts where the plaintext is known. This includes ciphertexts the attacker encrypted legitimately, as well as those derived by evaluating circuits of their choosing.
## Approximate Homomorphic Encryption Scheme (CKKS)
Classified as an *approximate decryption* scheme, the CKKS scheme is secure as long as the plaintext result of a decryption is only revealed to entities with knowledge of the secret-key. This is because, given a ciphertext 
$(-as +m + e, a)$
, the decryption outputs a plaintext $m + e$. [Li and Micciancio](https://eprint.iacr.org/2020/1533) show that using this plaintext, it is possible to recover the secret-key with $((- a s + m + e ) - ( m + e ) ) \cdot a^{-1} = a s a^{-1} = s$ (the probability of $a$ being invertible is overwhelming, and if $a$ is not invertible, only a few more samples are required).

This attack demonstrates that, when using an approximate homomorphic encryption scheme, the usual CPA security may not be sufficient depending on the application setting. Many applications do not require sharing the result with external parties and are not affected by this attack, but the ones that do must take the appropriate steps to ensure that no key-dependent information is leaked. A homomorphic encryption scheme that provides such functionality and that can be secure when releasing decrypted plaintext to external parties is defined to be CPA-D secure. The corresponding indistinguishability notion (IND-CPA-D) is defined as "indistinguishability under chosen plaintext attacks with decryption oracles".

The Lattigo API does _not_ provide automatic IND-CPA-D security. However, Lattigo implements the necessary API to mitigate _Li and Micciancio_'s attack. In particular, the decoding step of CKKS (and its real-number variant R-CKKS) allows the user to specify the desired fixed-point bit-precision.

Let $\epsilon$ be the scheme error after the decoding step. We compute the bit precision of the output as $\log_{2}(1/\epsilon)$.

If at any point of an application, decrypted values have to be shared with external parties, then the user must ensure that each shared plaintext is first _sanitized_ before being shared. To do so, the user must use the $\textsf{DecodePublic}$ method instead of the usual $\textsf{Decode}$. $\textsf{DecodePublic}$ takes as additional input the desired $\log_{2}(1/\epsilon)$-bit precision and rounds the value by evaluating $y = \lfloor x / \epsilon \rceil \cdot \epsilon$.

Estimating $\text{Pr}[\epsilon < x] \leq 2^{-s}$ of the circuit must be done carefully and we suggest the following process to do so:
 1. Given an error bound $\varepsilon$ (a user defined security parameter) and a circuit $C$ that takes as inputs length-$n$ vectors $\omega$ following a distribution $\chi$, select the appropriate parameters enabling the homomorphic evaluation of $C(\omega)$, denoted by $H(C(\omega))$, which includes the encoding, encryption, evaluation, decryption and decoding.
 2. Sample input vectors $\omega$ from the distribution $\chi$ and record $\epsilon = C(\omega) - H(C(\omega))$ for each slots. The user should make sure that the underlying circuit computed by $H(C(\cdot))$ is identical to $C(\cdot)$; i.e., if the homomorphic implementation $H(C(\cdot))$ uses polynomial approximations, then $C(\cdot)$ should use them too, instead of using the original exact function. Repeat until enough data points are collected to construct a CDF of $\textsf{Pr}[\epsilon > x]$.
 3. Use the CDF to select the value $\text{E}[\epsilon]$ such that any given slot will fail with probability $2^{-\varepsilon}$ to reach $\log_{2}(1/\epsilon)$ bits of precision. 
 4. Use the encoder method $\textsf{DecodePublic}$ with the parameter $\log_{2}(1/\epsilon)$ to decode plaintexts that will be published.

Note that, for composability with differential privacy, the variance of the error introduced by the rounding is $\text{Var}[x - \lfloor x \cdot \epsilon \rceil / \epsilon] = \tfrac{\epsilon^2}{12}$ and therefore $\text{Var}[x -  \lfloor x/(\sigma\sqrt{12})\rceil\cdot(\sigma\sqrt{12})] = \sigma^2$.

[Bossuat et al.](https://eprint.iacr.org/2024/853) provide tight bounds on noise to optimize the rounding process, minimizing loss in both precision and efficiency. Integration of the corresponding noise estimator in Lattigo is planned. 

## Exact Homomorphic Encryption Scheme
[Checri et al.](https://eprint.iacr.org/2024/116) and [Cheon et al.](https://eprint.iacr.org/2024/127), revealed new passive key-recovery attacks targeting also the exact FHE cryptosystems, including BFV, BGV, and TFHE. They exploit imperfect correctness and show that BFV, BGV and TFHE are not protected against IND-CPA-D attackers.

Achieving IND-CPA-D security for the exact homomorphic encryption schemes requires near-perfect correctness, meaning decryption failures must be exceptionally rare, with a probability lower than $2^{−\lambda}$, where $\lambda$ is a user-defined security parameter. Such failures should be so unlikely that finding one is computationally infeasible.
For exact schemes like BFV and BGV, implemented in Lattigo, near-perfect correctness can be maintained by adjusting scheme parameters to bound decryption noise, though this comes at the cost of performance. The circuit implementation must also control noise growth by limiting the number and type of operations performed at each computation level.

## Multiparty/Threshold Homomorphic Encryption

Multiparty or Threshold Fully Homomorphic Encryption involves secret-sharing the encryption key among multiple users, requiring their collaboration to decrypt data. The scheme maintains security even if a defined threshold number of users are compromised by an attacker. However, when decrypting, a party has access to the raw decrypted value before rounding, even with exact schemes. Since the receiver(s) (i.e. the decrypting party) may be corrupted, the attacker can mount the same kind of attacks, as described by [Checri et al.](https://eprint.iacr.org/2024/116) Fortunately, proper noise flooding as described by [Mouchet et al.](https://eprint.iacr.org/2020/304) (Section IV.E) thwarts these attacks. However, such methods require exponential noise, which will affect the performance and/or the correctness of the protocol.


## Limitations of the IND-CPA-D model

It is important to clarify that, for many applications, IND-CPA is sufficient.The IND-CPA-D model does not always align well with practical use cases. This model fails to address a broad range of real-world adversarial scenarios-for example, an attacker capable of producing or computing ciphertexts is unlikely to be limited to the constraints defined by the IND-CPA-D model. Furthermore, it introduces unnecessary performance drawbacks for FHE. However, depending on the specific application and threat model, additional security guarantees may be necessary. In such cases, these enhanced guarantees can often be achieved through supplementary application-specific countermeasures or protocol modifications.

## Recommendation for applicative countermeasures
1. As already said, the inherent malleability of FHE ciphertexts, along with weaknesses such as circular security and decision-to-search attacks, creates potential avenues for key-recovery attacks. As a foundational principle, it’s crucial that FHE ciphertexts are transmitted only through private and authenticated channels, encapsulated within traditional cryptographic methods.  
2. Most IND-CPA-D attacks require hundreds of thousands of queries to the evaluation and decryption oracles. By employing ephemeral keys or key rotations, we can limit the number of available plaintext-ciphertext pairs, effectively reducing the attack surface and mitigating these types of attacks.  
3.  Augmenting the FHE scheme with some kind of verifiable computation mechanism (e.g. SNARKs) might help protect against stronger adversaries (see e.g., [Canard et al.](https://eprint.iacr.org/2024/812)).  
4. Other physical limitations: firewall, rate control, enclave

# IND-CCA Security 

Fully Homomorphic Encryption (FHE) ciphertexts are malleable by design. This malleability, when combined with vulnerabilities such as circular security and decision-to-search attacks, can result in trivial key-recovery reaction attacks. This implies that FHE schemes are **not** secure against chosen ciphertext attacks (CCA security). Many intermediate security notions between IND-CPA and IND-CCA exist; an extensive summary can be found in the work of [Canard et al.](https://eprint.iacr.org/2024/812).

# Circuit Privacy
As discussed above, the FHE schemes proposed in Lattigo guarantee some notion of confidentiality (e.g. IND-CPA) but they do not necessarily hide information about the underlying computation. In other words, they do not provide **circuit privacy**. More precisely, schemes like CKKS and BFV/BGV leak information on the computation performed because the homomorphic operations introduce some structure in the ciphertext. Then, an adversary with knowledge of ciphertexts can infer details about the circuit (e.g. its depth, the type of operations, constants…). For instance, given a ciphertext $\text{ct}_1$ and its transformed version $\text{ct}_2$ after a homomorphic multiplication by a constant $C$, an adversary can deduce $C$.  

Circuit privacy usually requires techniques such as noise flooding or rerandomization, which in turn implies a significant loss in performance. 
