package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKEvaluator is a wrapper for the bfv evaluator
type MKEvaluator interface {
	Add(c1 *MKCiphertext, c2 *MKCiphertext, params *bfv.Parameters) *MKCiphertext
	MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters) *MKCiphertext
	MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, params *bfv.Parameters) *MKCiphertext
}

type mkEvaluator struct {
	bfvEval bfv.Evaluator
	ringQ   *ring.Ring
}

// NewMKEvaluator returns an evaluator for the multi key bfv scheme.
func NewMKEvaluator(params *bfv.Parameters) MKEvaluator {
	r := GetRingQ(params)
	return &mkEvaluator{bfvEval: bfv.NewEvaluator(params, bfv.EvaluationKey{}), ringQ: r}
}

// Add adds the cyphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Add(c1 *MKCiphertext, c2 *MKCiphertext, params *bfv.Parameters) *MKCiphertext {

	out := NewMKCiphertext(c1.peerIDs, eval.ringQ, params)

	val := make([]*ring.Poly, len(c1.peerIDs))

	for i := uint64(0); i < uint64(len(c1.peerIDs)); i++ {
		val[i] = eval.ringQ.NewPoly()
		eval.ringQ.Add(c1.ciphertexts.Value()[i], c2.ciphertexts.Value()[i], val[i])
	}

	out.ciphertexts.SetValue(val)

	return out
}

// MultSharedRelinKey will compute the homomorphic multiplication and relinearize the resulting cyphertext using pre computed Relin key
func (eval *mkEvaluator) MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters) *MKCiphertext {

	out := NewMKCiphertext(c1.peerIDs, eval.ringQ, params)

	eval.bfvEval.Mul(c1.ciphertexts.Element.Ciphertext(), c2.ciphertexts.Element.Ciphertext(), out.ciphertexts)

	// Call Relin on the resulting ciphertext
	RelinearizationWithSharedRelinKey(relinKey, out)

	return out
}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func (eval *mkEvaluator) MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, params *bfv.Parameters) *MKCiphertext {

	out := NewMKCiphertext(c1.peerIDs, eval.ringQ, params)

	// TODO: resize here since out is in Rq^(k+1)^2 ?  Mul seems to be incompatible with multi key ciphertexts....
	out.ciphertexts.Resize(params, c1.ciphertexts.Degree()+c2.ciphertexts.Degree())

	eval.bfvEval.Mul(c1.ciphertexts.Element.Ciphertext(), c2.ciphertexts.Element.Ciphertext(), out.ciphertexts)

	// Call Relin alg 2
	RelinearizationOnTheFly(evalKeys, publicKeys, out, params)

	return out
}

// tensor computes the tensor product between 2 ciphhertexts and returns the result in out
func (eval *mkEvaluator) tensor(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext) {
	// TODO implement tensor product
}

// quantize multiplies the values of an element by t/q
func (eval *mkEvaluator) quantize(ctOut *bfv.Element) { //TODO: adapt fromm evaluator.go

	levelQ := uint64(len(eval.ringQ.Modulus) - 1)
	levelQMul := uint64(len(eval.ringQMul.Modulus) - 1)

	c2Q1 := eval.poolQ[2]
	c2Q2 := eval.poolQmul[2] // what is the state of these variables after tensor ?

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.Value() {
		eval.ringQ.InvNTTLazy(c2Q1[i], c2Q1[i])
		eval.ringQMul.InvNTTLazy(c2Q2[i], c2Q2[i])

		// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
		eval.baseconverterQ1Q2.ModDownSplitQP(levelQ, levelQMul, c2Q1[i], c2Q2[i], c2Q2[i])

		// Centers (ct(x)Q -> P)/Q by (P-1)/2 and extends ((ct(x)Q -> P)/Q) to the basis Q
		eval.ringQMul.AddScalarBigint(c2Q2[i], eval.pHalf, c2Q2[i])
		eval.baseconverterQ1Q2.ModUpSplitPQ(levelQMul, c2Q2[i], ctOut.value[i])
		eval.ringQ.SubScalarBigint(ctOut.value[i], eval.pHalf, ctOut.value[i])

		// Option (2) (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
		eval.ringQ.MulScalar(ctOut.Value()[i], eval.t, ctOut.Value()[i])
	}

}
