## Description of the MKCKKS package
This package contains an implementation of the "special modulus" variant of the CKKS-based MKHE scheme proposed by Chen & al. in their 2019 paper: "Efficient Multi-Key Homomorphic Encryptionwith Packed Ciphertexts with Applicationto Oblivious Neural Network Inference".



### Participant

The Participant interface encapsulate all the local computation like Encryption, decryption and keys generation (including evaluation keys).
A standard error greater than 3.2 must be provided to create a new Participant. This will be used to compute the partial decryption of a ciphertext.

#### Use example

```go
        // get random value
		value := newTestValue(params, complex(-1, -1), complex(1, 1))

        // security parameter for partial decryption
        sigmaSmudging := 6.0 

        // create CRS and participant
        a := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)

		p := NewParticipant(params, sigmaSmudging, a)

		// encrypt
		cipher := p.Encrypt(value)

        //-----------Perform some homomorphic operations-------------

		// decrypt
		partialDec := p.GetPartialDecryption(cipher)
		decrypted := p.Decrypt(cipher, []*ring.Poly{partialDec})
```

It is possible to switch from the classic CKKS setting to a multi key setting by creating a participant from an already existing secret key.
It is also possible to use a ckks.Ciphertext and wrap it in a MKCiphertext.

```go

		ciphertext1 = encryptorPK.EncryptFastNew(plaintext)

		// setup keys and public parameters
		a := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)
		part1 := NewParticipantFromSecretKey(params, 6.0, a, sk)
		part2 := NewParticipant(params, 6.0, a)

		// perform addition with a mkckks ciphertext
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))
		ciphertext2 := part2.Encrypt(values2)
		evaluator := NewMKEvaluator(params)
		res := evaluator.Add(ciphertext2, &MKCiphertext{Ciphertexts: ciphertext1, PeerID: []uint64{part1.GetID()}})
```

### Evaluator

The evaluator is similar to the one in the ckks package. 

#### Use example

```go
        // Create Evaluator
		evaluator := NewMKEvaluator(params)

        // Gather public keys and evaluation keys from all participants involved
		evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey(), participants[2].GetEvaluationKey(), participants[3].GetEvaluationKey()}
		publicKeys := []*mkrlwe.MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey(), participants[2].GetPublicKey(), participants[3].GetPublicKey()}

        // Multiply
		resCipher1 := evaluator.Mul(cipher1, cipher2)
		resCipher2 := evaluator.Mul(cipher3, cipher4)

        // Relinearize
		evaluator.RelinInPlace(resCipher1, evalKeys[:2], publicKeys[:2])
		evaluator.RelinInPlace(resCipher2, evalKeys[2:], publicKeys[2:])

        // Drop the level and rescale to reduce the noise
		evaluator.DropLevel(resCipher1, 1)
		evaluator.Rescale(resCipher1, resCipher1)
		evaluator.DropLevel(resCipher2, 1)
		evaluator.Rescale(resCipher2, resCipher2)

        // Add the ciphhertexts resulting from the multiplication
		resCipher := evaluator.Add(resCipher1, resCipher2)

```
