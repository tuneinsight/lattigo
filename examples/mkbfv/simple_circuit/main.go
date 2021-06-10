package main

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkbfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// example of computation of a simple circuit with 2 parties
func main() {

	paramLit := bfv.DefaultParams[0]
	params, err := bfv.NewParametersFromLiteral(paramLit)

	if err != nil {
		panic("Couldn't retrieve default bfv parameters")
	}

	// generation of crs
	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}

	crs := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)

	// Plaintext ring
	ringT, err := ring.NewRing(params.N(), []uint64{params.T()})
	if err != nil {
		panic("Couldn't instanciate ringT")
	}

	// Participant 1
	keys1 := mkrlwe.KeyGen(&params.Parameters, crs)

	//create an encryptor from the first component of the public key
	pk1 := new(rlwe.PublicKey)
	pk1.Value[0] = keys1.PublicKey.Key[0].Poly[0] // b[0]
	pk1.Value[1] = keys1.PublicKey.Key[1].Poly[0] // a[0]

	encryptor1 := bfv.NewEncryptorFromPk(params, pk1)
	encoder1 := bfv.NewEncoder(params)
	decryptor1 := mkrlwe.NewMKDecryptor(&params.Parameters)

	value1 := mkrlwe.GetRandomPoly(&params.Parameters, ringT).Coeffs[0]
	plaintext1 := bfv.NewPlaintext(params)
	encoder1.EncodeUint(value1, plaintext1)

	// data that will be sent to the evaluator
	cipher1 := encryptor1.EncryptFastNew(plaintext1)
	evk1 := keys1.RelinKey
	pubkey1 := keys1.PublicKey

	// Participant 2
	keys2 := mkrlwe.KeyGen(&params.Parameters, crs)

	pk2 := new(rlwe.PublicKey)
	pk2.Value[0] = keys2.PublicKey.Key[0].Poly[0] // b[0]
	pk2.Value[1] = keys2.PublicKey.Key[1].Poly[0] // a[0]

	encryptor2 := bfv.NewEncryptorFromPk(params, pk2)
	encoder2 := bfv.NewEncoder(params)
	decryptor2 := mkrlwe.NewMKDecryptor(&params.Parameters)

	value2 := mkrlwe.GetRandomPoly(&params.Parameters, ringT).Coeffs[0]
	plaintext2 := bfv.NewPlaintext(params)
	encoder2.EncodeUint(value2, plaintext2)

	cipher2 := encryptor2.EncryptFastNew(plaintext2)
	evk2 := keys2.RelinKey
	pubkey2 := keys2.PublicKey

	// Evaluator: evaluates (c1 - c2) * (c1 + c2)
	evaluator := mkbfv.NewMKEvaluator(&params)

	// decide on an indexing method for the participants and their public material and ciphertexts
	ids := []uint64{1, 2}
	evk1.PeerID = 1
	evk2.PeerID = 2
	pubkey1.PeerID = 1
	pubkey2.PeerID = 2
	evalKeys := []*mkrlwe.MKRelinearizationKey{evk1, evk2}
	pubKeys := []*mkrlwe.MKPublicKey{pubkey1, pubkey2}

	// convert the bfv ciphertexts into multi key ciphertexts
	ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, ids)

	// evaluate circuit
	res1 := evaluator.Sub(ciphers[0], ciphers[1])
	res2 := evaluator.Add(ciphers[0], ciphers[1])
	res := evaluator.Mul(res1, res2)
	evaluator.RelinInPlace(res, evalKeys, pubKeys)

	// convert the multi key result into bfv ciphertexts for all participants
	resBFV := evaluator.ConvertToBFVCiphertext(res)

	// Partial Decryption done by each participants
	bfvCipher1 := resBFV[0]
	bfvCipher2 := resBFV[1]

	part1 := decryptor1.PartDec(bfvCipher1.Element, bfvCipher1.Level(), keys1.SecretKey, 6.0)
	part2 := decryptor2.PartDec(bfvCipher2.Element, bfvCipher2.Level(), keys2.SecretKey, 6.0)

	// Final decryption using the partial shares
	decrypted := decryptor1.MergeDec(bfvCipher1.Element, bfvCipher1.Level(), []*ring.Poly{part1, part2})

	// decode
	pt := bfv.NewPlaintext(params)
	pt.SetValue(decrypted)

	finalValues := encoder1.DecodeUintNew(pt)

	// perform the operation in plaintext space and check correctness

	p1 := ringT.NewPoly()
	p2 := ringT.NewPoly()
	tmp1 := ringT.NewPoly()
	tmp2 := ringT.NewPoly()

	copy(p1.Coeffs[0], value1)
	copy(p2.Coeffs[0], value2)

	ringT.Add(p1, p2, tmp1)
	ringT.Sub(p1, p2, tmp2)
	ringT.MulCoeffs(tmp1, tmp2, tmp1)

	for i, v := range tmp1.Coeffs[0] {
		if v != finalValues[i] {
			panic("Circuit Error")
		}
	}

	println("Computation Successful")

}
