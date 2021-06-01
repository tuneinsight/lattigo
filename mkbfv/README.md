## Description of the MKCKKS package
This package contains an implementation of the "special modulus" variant of the BFV-based MKHE scheme proposed by Chen & al. in their 2019 paper: "Efficient Multi-Key Homomorphic Encryptionwith Packed Ciphertexts with Applicationto Oblivious Neural Network Inference".

In this scheme, each participant run the key generation. Then they all encrypt their data using their private key and send the encrypted data and public material (public keys and evaluation keys) to an evaluator. The evaluator computes a circuit homomorphically and sends the result to the participants. The participants then have to collectively compute the decryption.

### Setup
In the multi-key setting, each participant generates its own pair of public and private key. Ciphertexts encrypted with the public key can then be used to homomorphically compute circuits. 
For multiplications, the evaluation key and the public key must be provided to the evaluator and for the rotation, the rotation key must be provided.

#### Setup Example

```go
	keys1 := mkrlwe.KeyGen(&params.Parameters, crs)
	encryptor1 := mkbfv.NewMKEncryptor(keys1.PublicKey, &params)
	encoder1 := bfv.NewEncoder(params)
	decryptor1 := mkrlwe.NewMKDecryptor(&params.Parameters, 0.6)

	value1 := mkrlwe.GetRandomPoly(&params.Parameters, ringT).Coeffs[0]
	plaintext1 := bfv.NewPlaintext(params)
	encoder1.EncodeUint(value1, plaintext1)

	cipher1 := encryptor1.Encrypt(plaintext1)
	evk1 := keys1.EvalKey
	pk1 := keys1.PublicKey

```

### Ciphertexts

The ciphertexts are the same as the one in the bfv package except in the evaluator. The evaluator uses ```MKCiphertexts```, ciphertexts containing data from multiple participants while the ciphertexts that comes out of the encryptor are classical ```bfv.Ciphertext```.

This makes it possible to compute something using the ```dbfv``` or ```bfv```package and then switch to the multi-key setting.

### Evaluator

The evaluator is similar to the one in the bfv in its usage. The only difference is that it converts the ```bfv.Ciphertext``` in ```mkbfv.MKCiphertexts``` using a conversion function. Then it must decide on an indexing method for each participant (this can be done using IP addresses, public keys, certificates etc...). At the end of the evaluation phase, the ciphertexts must be converted back to ```bfv.Ciphertext``` and then sent to the participants for the collective decryption procedure.

#### Evaluation example

```go
	// create an evaluator
	evaluator := mkbfv.NewMKEvaluator(&params)

	// decide on an indexing method for the participants and their public material and ciphertexts
	ids := []uint64{1, 2}
	evk1.PeerID = 1
	evk2.PeerID = 2
	pk1.PeerID = 1
	pk2.PeerID = 2
	evalKeys := []*mkrlwe.MKEvaluationKey{evk1, evk2}
	pubKeys := []*mkrlwe.MKPublicKey{pk1, pk2}

	// convert the bfv ciphertexts into multi key ciphertexts
	ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, ids)

	// evaluate circuit
	res1 := evaluator.Sub(ciphers[0], ciphers[1])
	res2 := evaluator.Add(ciphers[0], ciphers[1])
	res := evaluator.Mul(res1, res2)
	evaluator.RelinInPlace(res, evalKeys, pubKeys)

	// convert the multi key result into bfv ciphertexts for all participants
	resBFV := evaluator.ConvertToBFVCiphertext(res)
```

### Decryption

The decryption has two phases. First, all participants compute a share of the decryption using the ```MKDecryptor.PartDec``` function.
Then they send it to all other participants and merge all the shares using the ```MKDecryptor.MergeDec``` function to recover the final result.

#### Decryption example

```go
	part1 := decryptor1.PartDec(bfvCipher1.Element, bfvCipher1.Level(), keys1.SecretKey)
	part2 := decryptor2.PartDec(bfvCipher2.Element, bfvCipher2.Level(), keys2.SecretKey)

	// Final decryption using the partial shares
	decrypted := decryptor1.MergeDec(bfvCipher1.Element, bfvCipher1.Level(), []*ring.Poly{part1, part2})

	// decode
	pt := bfv.NewPlaintext(params)
	pt.SetValue(decrypted)

	finalValues := encoder1.DecodeUintNew(pt)
```


### Tests and Benchmarks

To run the tests simply type ```go test -v``` and to run the benchmarks type ```go test -bench MKBFV -run=^$ -benchmem -timeout 99999s```

### Performances

Relinearization and multiplications are quadratic in the number of participants in both time and memory.

### References

1. Efficient Multi-Key Homomorphic Encryption with Packed Ciphertext with Application to Oblivious Neural Network Inference (<https://eprint.iacr.org/2019/524>)