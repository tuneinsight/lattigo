# BGV

The BGV package provides a unified RNS-accelerated variant of the Fan-Vercauteren version of the Brakerski's scale invariant homomorphic encryption scheme (BFV) and Brakerski-Gentry-Vaikuntanathan (BGV) homomorphic encryption scheme. It enables SIMD modular arithmetic over encrypted vectors or integers.

## Implementation Notes

The proposed implementation is not standard and provides all the functionalities of the BFV and BGV schemes under a unfied scheme.
This enabled by the equivalency between the LSB and MSB encoding when T is coprime to Q (Appendix A of <https://eprint.iacr.org/2013/372>).

### Intuition

The textbook BGV scheme encodes the plaintext in the LSB and the encryption is done by the error by $T$:

$$\textsf{Encrypt}_{s}(\textsf{Encode}(m)) = [-as + m + Te, a]_{Q_{\ell}}$$ where $$Q_{\ell} = \prod_{i=0}^{L} q_{i}$$


 The decoding process is then carried out by taking the decrypted plaintext $[m + Te]_{Q_{\ell}}$ and taking it modulo $T$ which vanishes the error.

The only non-linear part of the BGV scheme is its modulus switch and that this operation is identical to a CKKS-style rescaling (quantization of the ciphertext by $\frac{1}{q_{\ell}}$) with a pre- and post-processing:

1) Multiply the ciphertext by $T^{-1}\mod Q_{\ell}$ (switch from LSB to MSB encoding)
2) Apply the CKKS-style rescaling (division by $q_{\ell}$)
3) Multiply the ciphertext by $T \mod Q_{\ell-1}$ (switch from MSB to LSB encoding)

Since the modulus switch is the only non-linear part of the BGV scheme, we can move this pre- and post- processing in the encoding step, i.e. instead of scaling the error by T we scale the plaintext by $T^{-1} mod Q_{\ell}$.

### Functionalities

The above change enables an implementation of the BGV scheme with an MSB encoding, which is essentially the BFV scheme. In other words, if $T\not|Q$ then the BFV and BGV schemes are indistinguishable up to a plaintext scaling factor of $T^{-1}\mod Q$. 

It can also be seen as a variant of the BGV scheme with two tensoring operations:
- The BGV-style tensoring with a noise growth proportional to the current noise
- The BFV-style tensoring with a noise growth invariant to the current noise