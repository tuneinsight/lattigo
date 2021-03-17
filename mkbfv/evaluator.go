package mkbfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKEvaluator is a wrapper for the bfv evaluator
type MKEvaluator interface {
	Add(c1 *MKCiphertext, c2 *MKCiphertext, params *bfv.Parameters) *MKCiphertext
	MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters) *MKCiphertext
	MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, params *bfv.Parameters) *MKCiphertext
}

type mkEvaluator struct {
	bfvEval         bfv.Evaluator
	params          *bfv.Parameters
	ringQ           *ring.Ring
	ringQMul        *ring.Ring
	pHalf           *big.Int
	samplerGaussian *ring.GaussianSampler
	polyPoolQ1      []*ring.Poly
	polyPoolQ2      []*ring.Poly
	convertor       *ring.FastBasisExtender
}

// NewMKEvaluator returns an evaluator for the multi key bfv scheme.
func NewMKEvaluator(params *bfv.Parameters) MKEvaluator {

	ringQ := GetRingQ(params)
	ringQMul := GetRingQMul(params)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	sampler := GetGaussianSampler(params, ringQ, prng)
	convertor := ring.NewFastBasisExtender(ringQ, ringQMul)

	pHalf := new(big.Int).Rsh(ringQMul.ModulusBigint, 1)

	return &mkEvaluator{
		bfvEval:         bfv.NewEvaluator(params, bfv.EvaluationKey{}),
		params:          params,
		ringQ:           ringQ,
		ringQMul:        ringQMul,
		pHalf:           pHalf,
		samplerGaussian: sampler,
		convertor:       convertor}
}

// Add adds the ciphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Add(c1 *MKCiphertext, c2 *MKCiphertext, params *bfv.Parameters) *MKCiphertext {

	out := NewMKCiphertext(c1.peerIDs, eval.ringQ, params)

	val := make([]*ring.Poly, len(c1.peerIDs)+1)

	for i := uint64(0); i < uint64(len(c1.peerIDs)+1); i++ {
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

func (eval *mkEvaluator) modUpAndNTT(ct *bfv.Element, cQ, cQMul []*ring.Poly) {
	levelQ := uint64(len(eval.ringQ.Modulus) - 1)
	for i := range ct.Value() {
		eval.convertor.ModUpSplitQP(levelQ, ct.Value()[i], cQMul[i])
		eval.ringQ.NTTLazy(ct.Value()[i], cQ[i])
		eval.ringQMul.NTTLazy(cQMul[i], cQMul[i])
	}
}

// tensor computes the tensor product between 2 ciphertexts and returns the result in out
// c1 and c2 must have be of dimension k+1, where k = #participants
// out has dimensions (k+1)**2
func (eval *mkEvaluator) tensor(ct0, ct1 *bfv.Element, out *bfv.Element) {

	outputDegree := (ct0.Degree() + 1) * (ct1.Degree() + 1)

	c0Q1 := make([]*ring.Poly, ct0.Degree()+1)
	c0Q2 := make([]*ring.Poly, ct0.Degree()+1)

	for i := uint64(0); i < ct0.Degree()+1; i++ {
		c0Q1[i] = eval.ringQ.NewPoly()
		c0Q2[i] = eval.ringQMul.NewPoly()
	}

	eval.modUpAndNTT(ct0, c0Q1, c0Q2) // split ct0 in ringQ and ringQMul

	c2Q1 := make([]*ring.Poly, outputDegree) // prepare output
	c2Q2 := make([]*ring.Poly, outputDegree)

	for i := uint64(0); i < outputDegree; i++ {

		c2Q1[i] = eval.ringQ.NewPoly()
		c2Q2[i] = eval.ringQMul.NewPoly()
	}

	// Squaring case
	if ct0 == ct1 {

		c00Q1 := make([]*ring.Poly, ct0.Degree()+1)
		c00Q2 := make([]*ring.Poly, ct0.Degree()+1)

		for i := range ct0.Value() {

			c00Q1[i] = eval.ringQ.NewPoly()
			c00Q2[i] = eval.ringQMul.NewPoly()

			eval.ringQ.MForm(c0Q1[i], c00Q1[i])
			eval.ringQMul.MForm(c0Q2[i], c00Q2[i])
		}

		for i := uint64(0); i < ct0.Degree()+1; i++ {
			for j := i + 1; j < ct0.Degree()+1; j++ {
				eval.ringQ.MulCoeffsMontgomery(c00Q1[i], c0Q1[j], c2Q1[i+j])
				eval.ringQMul.MulCoeffsMontgomery(c00Q2[i], c0Q2[j], c2Q2[i+j])

				eval.ringQ.Add(c2Q1[i+j], c2Q1[i+j], c2Q1[i+j])
				eval.ringQMul.Add(c2Q2[i+j], c2Q2[i+j], c2Q2[i+j])
			}
		}

		for i := uint64(0); i < ct0.Degree()+1; i++ {
			eval.ringQ.MulCoeffsMontgomeryAndAdd(c00Q1[i], c0Q1[i], c2Q1[i<<1])
			eval.ringQMul.MulCoeffsMontgomeryAndAdd(c00Q2[i], c0Q2[i], c2Q2[i<<1])
		}

		// Normal case
	} else {

		c1Q1 := make([]*ring.Poly, ct1.Degree()+1)
		c1Q2 := make([]*ring.Poly, ct1.Degree()+1)

		for i := uint64(0); i < ct1.Degree()+1; i++ {
			c1Q1[i] = eval.ringQ.NewPoly()
			c1Q2[i] = eval.ringQMul.NewPoly()
		}

		eval.modUpAndNTT(ct1, c1Q1, c1Q2)

		for i := range ct0.Value() {
			eval.ringQ.MForm(c0Q1[i], c0Q1[i])
			eval.ringQMul.MForm(c0Q2[i], c0Q2[i])
			for j := range ct1.Value() {
				eval.ringQ.MulCoeffsMontgomeryAndAdd(c0Q1[i], c1Q1[j], c2Q1[i+j])
				eval.ringQMul.MulCoeffsMontgomeryAndAdd(c0Q2[i], c1Q2[j], c2Q2[i+j])
			}
		}
	}
}

// quantize multiplies the values of an element by t/q
func (eval *mkEvaluator) quantize(c2Q1, c2Q2 []*ring.Poly, ctOut *bfv.Element) {

	levelQ := uint64(len(eval.ringQ.Modulus) - 1)
	levelQMul := uint64(len(eval.ringQMul.Modulus) - 1)

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.Value() { // will iterate on (k + 1)^2 values.. TODO: check corectness !!!!!
		eval.ringQ.InvNTTLazy(c2Q1[i], c2Q1[i])
		eval.ringQMul.InvNTTLazy(c2Q2[i], c2Q2[i])

		// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
		eval.convertor.ModDownSplitQP(levelQ, levelQMul, c2Q1[i], c2Q2[i], c2Q2[i])

		// Centers (ct(x)Q -> P)/Q by (P-1)/2 and extends ((ct(x)Q -> P)/Q) to the basis Q
		eval.ringQMul.AddScalarBigint(c2Q2[i], eval.pHalf, c2Q2[i])
		eval.convertor.ModUpSplitPQ(levelQMul, c2Q2[i], ctOut.Value()[i])
		eval.ringQ.SubScalarBigint(ctOut.Value()[i], eval.pHalf, ctOut.Value()[i])

		// Option (2) (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
		eval.ringQ.MulScalar(ctOut.Value()[i], eval.params.T(), ctOut.Value()[i])
	}

}
