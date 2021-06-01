## Description of the MKCKKS package
This package contains an implementation of the "special modulus" variant of the CKKS/BFV-based MKHE scheme proposed by Chen & al. in their 2019 paper: "Efficient Multi-Key Homomorphic Encryptionwith Packed Ciphertexts with Applicationto Oblivious Neural Network Inference".

This package contain code that is common to the ```mkbfv``` and ```mkckks``` package such as keys, decryption algorithm, key generation and relinearization.

### Relinearization

The relinearization algorithm is described in annexe A of the paper. Both the evaluation key and the public key must be provided to the relinearization algorithm.

### References

1. Efficient Multi-Key Homomorphic Encryption with Packed Ciphertext with Application to Oblivious Neural Network Inference (<https://eprint.iacr.org/2019/524>)