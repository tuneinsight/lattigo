package mkbfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// Add adds the cyphertexts component wise and expend their list of involved peers
func Add(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, ringQ *ring.Ring, params *bfv.Parameters) {

	eval := bfv.NewEvaluator(params, bfv.EvaluationKey{}) // declare evaluator to get access to homomorphic operations

	eval.Add(c1.ciphertexts.Element, c2.ciphertexts.Element, out.ciphertexts)

}

// MultSharedRelinKey will compute the homomorphic multiplication and relinearize the resulting cyphertext using pre computed Relin key
func MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters) {

	tensorAndRescale(c1, c2, out, params)

	// Call Relin on the resulting ciphertext
	RelinearizationWithSharedRelinKey(relinKey, out)

}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func MultRelinDynamic() {
	// TODO: implement multiplication

	// Call Relin alg 2
}

func tensorAndRescale(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, params *bfv.Parameters) {

	ringQ := GetRingQ(params)

	// Compute tensor product of ciphertexts (goes from Rq^(k+1) -> Rq^(k+1)^2 )
	for _, c := range c1.ciphertexts.Value() {
		for _, cIN := range c2.ciphertexts.Value() { // TODO: take tensor product implemented in evaluator.go ?

			current := ringQ.NewPoly()
			ringQ.MulCoeffs(c, cIN, current)

			out.ciphertexts.Element.SetValue(append(out.ciphertexts.Element.Value(), current))
		}
	}

	// rescale by a factor t/q
	quantize(out.ciphertexts.Element, params)

}

// multiply by t/q the given ciphertext
func quantize(ctOut *bfv.Element, params *bfv.Parameters) {

	t := params.T()
	ringQ := GetRingQ(params)
	ringQMul := GenRingQMul(params)

	levelQ := uint64(len(ringQ.Modulus) - 1)
	levelQMul := uint64(len(ringQMul.Modulus) - 1)

	pHalf := new(big.Int).Rsh(ringQMul.ModulusBigint, 1)

	convertor := ring.NewFastBasisExtender(ringQ, ringQMul)

	c2Q1 := make([]*ring.Poly, 4)
	c2Q2 := make([]*ring.Poly, 4)

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.Value() {

		ringQ.InvNTTLazy(c2Q1[i], c2Q1[i])
		ringQMul.InvNTTLazy(c2Q2[i], c2Q2[i])

		// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
		convertor.ModDownSplitQP(levelQ, levelQMul, c2Q1[i], c2Q2[i], c2Q2[i])

		// Centers (ct(x)Q -> P)/Q by (P-1)/2 and extends ((ct(x)Q -> P)/Q) to the basis Q
		ringQMul.AddScalarBigint(c2Q2[i], pHalf, c2Q2[i])
		convertor.ModUpSplitPQ(levelQMul, c2Q2[i], ctOut.Value()[i])
		ringQ.SubScalarBigint(ctOut.Value()[i], pHalf, ctOut.Value()[i])

		// Option (2) (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
		ringQ.MulScalar(ctOut.Value()[i], t, ctOut.Value()[i])
	}
}
