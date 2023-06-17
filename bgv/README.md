# BGV

The BGV package provides a unified RNS-accelerated variant of the Fan-Vercauteren version of the Brakerski's scale invariant homomorphic encryption scheme (BFV) and Brakerski-Gentry-Vaikuntanathan (BGV) homomorphic encryption scheme. It enables SIMD modular arithmetic over encrypted vectors or integers.

## Implementation Notes

The proposed implementation is not standard and provides all the functionalities of the BFV and BGV schemes under a unfied scheme.
This enabled by the equivalency between the LSB and MSB encoding when T is coprime to Q (Appendix A of <https://eprint.iacr.org/2013/372>).

### Intuition

The textbook BGV scheme encodes the plaintext in the LSB and scales the error by T. The decoding process is then carried out by taking the decrypted plaintext (which is modulo Q, the ciphertext modulus) and taking it modulo T, which vanishes the error.

The only non-linear part of the BGV scheme is its modulus switch and tha this operation is identical to a CKKS-style rescaling (quantization of the ciphertext by 1/qi) with a pre- and post-processing:

1) Multiply the ciphertext by T^{-1} mod Q (switch from LSB to MSB encoding)
2) Apply the CKKS-style rescaling (truncate the lower bits)
3) Multiply the ciphertext by T mod Q (switch from MSB to LSB encoding)

Since the modulus switch is the only non-linear part of the BGV scheme, we can move this pre- and post- processing in the encoding step, i.e. instead of scaling the error by T we scale the plaintext by T^{-1} mod Q.

### Functionalities

The above change enables an implementation of the BGV scheme with an MSB encoding, which is essentially the BFV scheme. In other words, if T is coprime to Q then the BFV and BGV schemes are indistinguishable. 

It can also be seen as a variant of the BGV scheme with two tensoring operations:
- The BGV-style tensoring with a noise growth proportional to the current noise
- The BFV-style tensoring with a noise growth invariant to the current noise