package main

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkbfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
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
	encryptor1 := mkbfv.NewMKEncryptor(keys1.PublicKey, &params)
	encoder1 := bfv.NewEncoder(params)
	decryptor1 := mkrlwe.NewMKDecryptor(&params.Parameters, 0.6)

	value1 := mkrlwe.GetRandomPoly(&params.Parameters, ringT).Coeffs[0]
	plaintext1 := bfv.NewPlaintext(params)
	encoder1.EncodeUint(value1, plaintext1)

	cipher1 := encryptor1.Encrypt(plaintext1)
	evk1 := keys1.EvalKey
	pk1 := keys1.PublicKey

	// Participant 2
	keys2 := mkrlwe.KeyGen(&params.Parameters, crs)
	encryptor2 := mkbfv.NewMKEncryptor(keys2.PublicKey, &params)
	encoder2 := bfv.NewEncoder(params)
	decryptor2 := mkrlwe.NewMKDecryptor(&params.Parameters, 0.6)

	value2 := mkrlwe.GetRandomPoly(&params.Parameters, ringT).Coeffs[0]
	plaintext2 := bfv.NewPlaintext(params)
	encoder2.EncodeUint(value2, plaintext2)

	cipher2 := encryptor2.Encrypt(plaintext2)
	evk2 := keys2.EvalKey
	pk2 := keys2.PublicKey

	// Evaluator: evaluates (c1 - c2) * (c1 + c2)
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

	// Partial Decryption done by each participants
	bfvCipher1 := resBFV[0]
	bfvCipher2 := resBFV[1]

	part1 := decryptor1.PartDec(bfvCipher1.Element, bfvCipher1.Level(), keys1.SecretKey)
	part2 := decryptor2.PartDec(bfvCipher2.Element, bfvCipher2.Level(), keys2.SecretKey)

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
