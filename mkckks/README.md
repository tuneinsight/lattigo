## Description of the MKCKKS package
This package contains an implementation of the "special modulus" variant of the CKKS-based MKHE scheme proposed by Chen & al. in their 2019 paper: "Efficient Multi-Key Homomorphic Encryptionwith Packed Ciphertexts with Application to Oblivious Neural Network Inference".

### What is a multi-key homomorphic encryption scheme

A MKHE scheme is a type of MPC protocol in which each participant encrypts its input to the circuit with its own public key. An evaluator can then compute a function of these inputs and returns an encrypted result. A decryption algorithm is then run to recover the plaintext result.

In this specific scheme, the decryption algorithm is a multi-party protocol in the sense that one participant cannot decrypt the result without exchanging its decryption share with the other participants.

### Setup
In the multi-key setting, each participant generates its own pair of public and private key. Ciphertexts encrypted with the public key can then be used to homomorphically compute circuits. 
For multiplications, the relinearization key and the public key must be provided to the evaluator and for the rotation, the rotation key must be provided.

The ```mkrlwe.KeyGen``` algorithm can be used to generate the keys. The ckks parameters as well as the common reference string must be provided.
The crs has to be shared with all participants:
```go
// the KeyGen algorithm generates the secret, public and relinearization key
keys := mkrlwe.KeyGen(&params.Parameters, crs)
```


### Encryption

The encryption is the same as in the ```ckks``` package. Since the public keys returned by ```mkrlwe.KeyGen``` are decomposed in the RNS basis, only the first component of the ```mkrlwe.MKPublicKey``` must be used to create the encryptor:
```go
//create an encryptor from the first component of the public key
pk := new(rlwe.PublicKey)
pk.Value[0] = keys.PublicKey.Key[0].Poly[0] // b[0]
pk.Value[1] = keys.PublicKey.Key[1].Poly[0] // a[0]

encryptor := ckks.NewEncryptorFromPk(params, pk)
```
It is important to verify that the public key that is sent to the evaluator is the one returned by the ```mkrlwe.KeyGen``` and not the temporary one created from the first compoenents and only used to create the ```ckks.Encryptor```.


Alongside the ```ckks.Encryptor```, a ```ckks.Encoder``` can be used in order to, first, wrap the values in a ```ckks.Plaintext``` and then encrypt them into standard ```ckks.Ciphertexts```:
```go
encryptor := ckks.NewEncryptorFromPk(params, pk)
encoder := ckks.NewEncoder(params)

value := newTestValue(&params, complex(-1, -1), complex(1, 1))
plaintext := encoder.EncodeNTTAtLvlNew(params.MaxLevel(), value, params.LogSlots())

cipher := encryptor2.EncryptFastNew(plaintext)
```

It is also possible to simply use a ```ckks.KeyGenerator``` and a ```ckks.Encryptor``` to create the ciphertext if one wants to reuse an already existing ```rlwe.SecretKey```:
```go
keygen := ckks.NewKeyGenerator(*params)
sk, pk := keygen.GenKeyPair()
encryptorPK := ckks.NewEncryptorFromPk(*params, pk)
ciphertext = encryptorPK.EncryptNew(plaintext)
```

But, since the public key and relinearization key of the multi-key scheme are not compatible with the CKKS scheme, it is necessary to use the ```mkbfv.KeyGenWithSecretKey``` function to create the missing keys from the ```rlwe.SecretKey```:
```go
// a is the crs common to all participants
keys := mkrlwe.KeyGenWithSecretKey(&params.Parameters, a, sk) 
```

### Ciphertexts

The ciphertexts are the same as the one in the ckks package except in the evaluator. The evaluator uses ```MKCiphertexts```, ciphertexts containing data from multiple participants while the ciphertexts that comes out of the encryptor are classical ```ckks.Ciphertext```.

This makes it possible to compute something using the ```dckks``` or ```ckks```package and then switch to the multi-key setting.

The ```MKCiphertexts``` used by the evaluator grow linealry with the number of participants. Each one of them is made up of c0 (containing the encrypted values) and all the ciphertext "shares" of the involved participants.
For example, the result of a binary operation between two freshly encrypted ciphertexts (c01, c1)  and (c02, c2) result in a ```MKCiphertexts``` (c0, c1, c2).

### Evaluator

The evaluator is similar to the one in the ckks in its usage. The only difference is that it converts the ```ckks.Ciphertext``` in ```mkckks.MKCiphertexts``` using a conversion function. Then it must decide on an indexing method for each participant (this can be done using IP addresses, public keys, certificates etc...). 
At the end of the evaluation phase, the ciphertexts must be converted back to ```ckks.Ciphertext``` and then sent to the participants for the collective decryption procedure.

Creation of an evaluator:
```go
// create an evaluator
evaluator := mkckks.NewMKEvaluator(&params)
```

Simple circuit evaluation (a+b)*(a-b):
```go
res1 := evaluator.Sub(ciphers[0], ciphers[1])
res2 := evaluator.Add(ciphers[0], ciphers[1])
res := evaluator.Mul(res1, res2)
evaluator.RelinInPlace(res, evalKeys, pubKeys)
```

It is important to index all material (ciphertexts, keys) received by the participants or else the homomorphic operations won't work. The choice of the indexing method is left to the programmer:
```go
rlk1.PeerID = 1
rlk2.PeerID = 2

pk1.PeerID = 1
pk2.PeerID = 2

relinKeys := []*mkrlwe.MKRelinearizationKey{rlk1, rlk2}
pubKeys := []*mkrlwe.MKPublicKey{pk1, pk2}
```

The conversion from ```ckks.Ciphertexts``` to ```mkckks.MKCiphertext``` and back can be easily performed by the evaluator.

```go
// convert the ckks ciphertexts into multi key ciphertexts
ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, ids)

// evaluate circuit....

// convert the multi key result into ckks ciphertexts for all participants before sending back the results
resCKKS := evaluator.ConvertToCKKSCiphertext(res)
```

### Decryption

The decryption has two phases. First, all participants compute a share of the decryption using the ```MKDecryptor.PartDec``` function.
Then they send it to all other participants and merge all the shares using the ```MKDecryptor.MergeDec``` function to recover the final result.

Partial decryption:
```go
part1 := decryptor1.PartDec(ckksCipher1.Element, ckksCipher1.Level(), keys1.SecretKey)
part2 := decryptor2.PartDec(ckksCipher2.Element, ckksCipher2.Level(), keys2.SecretKey)
```

Final decryption and decoding done after collecting all the decryption sahres:
```go
// Final decryption using the partial shares
decrypted := decryptor.MergeDec(ckksCipher1.Element, ckksCipher1.Level(), []*ring.Poly{part1, part2})

// decode
pt := ckks.NewPlaintext(params, ckksCipher1.Level(), ckksCipher1.Scale())
pt.SetValue(decrypted)

finalValues := encoder.Decode(pt, params.LogSlots())
```

### Examples and use cases 

To see an example and use case of the mkckks package, under ```lattigo/examples/mkckks``` you can find a simple circuit using the mkckks scheme.

### Tests and Benchmarks

To run the tests simply type ```go test -v``` and to run the benchmarks type ```go test -bench MKCKKS -run=^$ -benchmem -timeout 99999s```

### Performances

Here is a short table with the asymptotic memory and time consumption of the different operations with respect to the number of participants.

| Operation  | Asymptotic Time | Asymptotic Memory |
| ------------- | ------------- | ------------- |
| Addition/Subtraction  | O(n)  | O(n) |
| Multiplication  | O(n^2)  | O(n^2)  |
| Relinearization  | O(n^2)  | O(n^2)  |
| Rotation  | O(n)  | O(n)  |
| Partial Decryption  | O(1)  | O(1)  |
| Merge Decryption  | O(n)  | O(n)  |


### References

1. Efficient Multi-Key Homomorphic Encryption with Packed Ciphertext with Application to Oblivious Neural Network Inference (<https://eprint.iacr.org/2019/524>)