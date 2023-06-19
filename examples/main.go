package main

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func main() {
	var params bfv.Parameters
	var err error
	if params, err = bfv.NewParametersFromLiteral(bfv.ParametersLiteral{
		LogN: 14,
		Q:    []uint64{0x3fffffa8001, 0x1000090001, 0x10000c8001, 0x10000f0001, 0xffff00001},
		P:    []uint64{0x7fffffd8001},
		T:    0xf60001,
	}); err != nil {
		panic(err)
	}

	kgen := bfv.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	ecd := bfv.NewEncoder(params)
	enc := bfv.NewEncryptor(params, pk)
	dec := bfv.NewDecryptor(params, sk)

	ct2 := bfv.NewCiphertext(params, 2, params.MaxLevel())

	eval := bfv.NewEvaluator(params, nil)

	var variance, min, max float64

	n := 1024
	for i := 0; i < n; i++ {
		coeffs0, ct0 := NewTestVector(params, ecd, enc)
		coeffs1, ct1 := NewTestVector(params, ecd, enc)

		eval.Mul(ct0, ct1, ct2)
		params.RingT().MulCoeffsBarrett(coeffs0, coeffs1, coeffs0)

		v, mi, ma := Noise(ct2, ecd, dec, eval)

		variance += v
		min += mi
		max += ma

		fmt.Println(i, variance, min, max)
	}

}

func NewTestVector(params bfv.Parameters, ecd *bfv.Encoder, enc rlwe.EncryptorInterface) (coeffs *ring.Poly, ct *rlwe.Ciphertext) {

	var prng sampling.PRNG
	var err error
	if prng, err = sampling.NewPRNG(); err != nil {
		panic(err)
	}

	uSampler := ring.NewUniformSampler(prng, params.RingT())

	coeffs = uSampler.ReadNew()
	pt := bfv.NewPlaintext(params, params.MaxLevel())
	ecd.Encode(coeffs.Coeffs[0], pt)
	ct = enc.EncryptNew(pt)
	return
}

func Noise(ct *rlwe.Ciphertext, ecd *bfv.Encoder, dec *rlwe.Decryptor, eval *bfv.Evaluator) (float64, float64, float64) {

	have := make([]uint64, 1<<14)
	pt := dec.DecryptNew(ct)
	ecd.Decode(pt, have)
	ecd.Encode(have, pt)
	return rlwe.Norm(eval.SubNew(ct, pt), dec)
}
