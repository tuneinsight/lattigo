## Description of the MKRLWE package
This package contains an implementation of the "special modulus" variant of the CKKS/BFV-based MKHE scheme proposed by Chen & al. in their 2019 paper: "Efficient Multi-Key Homomorphic Encryptionwith Packed Ciphertexts with Applicationto Oblivious Neural Network Inference".

This package contain code that is common to the ```mkbfv``` and ```mkckks``` package such as keys, decryption algorithm, key generation and relinearization.


### Keys

This package contains an implementation of several keys needed in MKHE schemes. All the keys in the MK setting hold an identifier for the participant holding the key, which can be determined by the programmer. The following keys can be found : 
 - Public Keys : A pair (public key, crs) represented by the type ```mkrlwe.MKPublicKey```
 
 - Secret Keys : the type ```mkrlwe.MKSecretKey``` wraps the type ```rlwe.SecretKey```, along with an identifier for the participant
 
 - Relinearization Keys : the type ```mkrlwe.RelinearizationKey``` wraps the type ```rlwe.SwitchingKey```, encapsulating the vectors d0 and d2 in [1] along with a vector representing the vector d1 in [1]. These keys are used in the context of relinearization of extended ciphertexts after a multiplication.
 
 - Rotation Keys : the type ```mkrlwe.RotationKey``` wraps the type ```rlwe.SwitchingKey```, along with an identifier for the participant. The keys are useful to perform the rotation operation on multi-key ciphertexts based on the evaluation of Galois Automorphism. 



### Key Generation

The ```mkrlwe.KeyGen``` algorithm can be used to generate the keys in packages implementing rlwe encryption schemes. Rotation Keys need to be generated separately with the rotation value specified : for each rotation value, one rotation key is generated. 
The key generation procedures in the implementation use **Modulus Raising** for performance optimization.


### Relinearization

Relinearization implements the algorithm 3 in the appendix of [1], using **Modulus Raising**. 
The algorithm linearizes each entry of an extended ciphertext by multiplying it with the public key and relinearization key of the participants involved in the generation of the ```MKCiphertext```. This method is more efficient in terms of storage and noise growth than a faster implementation with a shared relinearization key between participants [1].


### Decryption

The decryption has two phases. First, all participants compute a share of the decryption using the ```MKDecryptor.PartDec``` function.
Then they send it to all other participants and merge all the shares using the ```MKDecryptor.MergeDec``` function to recover the final result.
Decryption example code snippets can be found in the README files for schemes that build upon the mkrlwe primitives, such as mkckks and mkbfv, under ```lattigo/mkckks``` and ```lattigo/mkbfv``` 




### References

1. Efficient Multi-Key Homomorphic Encryption with Packed Ciphertext with Application to Oblivious Neural Network Inference (<https://eprint.iacr.org/2019/524>)