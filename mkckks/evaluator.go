package mkckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKEvaluator is a wrapper for the ckks evaluator
type MKEvaluator interface {
	Add(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey) *MKCiphertext
	Rescale(c *MKCiphertext, out *MKCiphertext)
}

type mkEvaluator struct {
	ckksEval        ckks.Evaluator
	params          *ckks.Parameters
	ringQ           *ring.Ring
	ringQMul        *ring.Ring
	pHalf           *big.Int
	samplerGaussian *ring.GaussianSampler
	polyPoolQ1      []*ring.Poly
	polyPoolQ2      []*ring.Poly
	convertor       *ring.FastBasisExtender
}

// NewMKEvaluator returns an evaluator for the multi key ckks scheme.
func NewMKEvaluator(params *ckks.Parameters) MKEvaluator {

	if params == nil {
		panic("Cannot create evaluator with uninitilized parameters")
	}

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
		ckksEval:        ckks.NewEvaluator(params, ckks.EvaluationKey{}),
		params:          params,
		ringQ:           ringQ,
		ringQMul:        ringQMul,
		pHalf:           pHalf,
		samplerGaussian: sampler,
		convertor:       convertor}
}

// Add adds the ciphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Add(c0 *MKCiphertext, c1 *MKCiphertext) *MKCiphertext {

	if c0 == nil || c1 == nil || c0.ciphertexts == nil || c1.ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}
	if c0.ciphertexts.Degree() != c1.ciphertexts.Degree() {
		panic("Ciphertexts must be of same degree before addition")
	}

	out := NewMKCiphertext(c0.peerIDs, eval.ringQ, eval.params, c0.ciphertexts.Level())

	var tmp0, tmp1 *ckks.Element

	level := utils.MinUint64(utils.MinUint64(c0.ciphertexts.Level(), c1.ciphertexts.Level()), out.ciphertexts.Level())

	val := make([]*ring.Poly, len(c1.peerIDs)+1)
	scale0 := c0.ciphertexts.Scale()
	scale1 := c1.ciphertexts.Scale()

	if scale1 > scale0 {

		tmp0 = ckks.NewElement()

		if uint64(scale1/scale0) != 0 {
			eval.ckksEval.MultByConst(c0.ciphertexts.Ciphertext(), uint64(scale1/scale0), tmp0.Ciphertext())
		}

		tmp1 = c1.ciphertexts.El()

	} else if scale0 > scale1 {

		tmp1 = ckks.NewElement()

		if uint64(scale0/scale1) != 0 {
			eval.ckksEval.MultByConst(c1.ciphertexts.Ciphertext(), uint64(scale0/scale1), tmp1.Ciphertext())
		}

		tmp0 = c0.ciphertexts.El()

	} else {
		tmp0 = c0.ciphertexts.El()
		tmp1 = c1.ciphertexts.El()
	}

	for i := uint64(0); i < c0.ciphertexts.Degree()+1; i++ {
		val[i] = eval.ringQ.NewPoly()
		eval.ringQ.AddLvl(level, tmp0.Value()[i], tmp1.Value()[i], val[i])
	}

	out.ciphertexts.SetScale(utils.MaxFloat64(scale0, scale1))
	out.ciphertexts.SetValue(val)
	out.peerIDs = c0.peerIDs

	return out
}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func (eval *mkEvaluator) MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey) *MKCiphertext {

	nbrElements := c1.ciphertexts.Degree() + 1 // k+1

	outputDegree := nbrElements * nbrElements // (k+1)**2

	el1 := c1.ciphertexts.Element
	el2 := c2.ciphertexts.Element
	level := utils.MinUint64(el1.Level(), el2.Level())

	out := new(MKCiphertext)
	out.ciphertexts = ckks.NewCiphertext(eval.params, outputDegree-1, level, el1.Scale()*el2.Scale())
	out.peerIDs = c1.peerIDs

	if !el1.IsNTT() {
		panic("cannot MulRelinDynamic: op0 must be in NTT")
	}

	if !el2.IsNTT() {
		panic("cannot MulRelinDynamic: op1 must be in NTT")
	}

	ringQ := eval.ringQ

	tmp1 := ringQ.NewPoly()
	tmp2 := ringQ.NewPoly()

	for i, v1 := range c1.ciphertexts.Value() {

		for j, v2 := range c2.ciphertexts.Value() {

			ringQ.MFormLvl(level, v1, tmp1)
			ringQ.MFormLvl(level, v2, tmp2)

			ringQ.MulCoeffsMontgomeryLvl(level, tmp1, tmp2, out.ciphertexts.Ciphertext().Value()[int(nbrElements)*i+j])
		}
	}

	// Call Relin alg 2
	RelinearizationOnTheFly(evalKeys, publicKeys, out, eval.params)

	return out
}

// Rescale takes a ciphertext at level l reduces it until it reaches its original scale
// this function is the same as in ckks/evaluator.go
func (eval *mkEvaluator) Rescale(c *MKCiphertext, out *MKCiphertext) {

	eval.ckksEval.Rescale(c.ciphertexts, eval.params.Scale(), c.ciphertexts)
}
