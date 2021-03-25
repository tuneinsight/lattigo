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
	MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey) *MKCiphertext
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

// MultSharedRelinKey will compute the homomorphic multiplication and relinearize the resulting cyphertext using pre computed Relin key
func (eval *mkEvaluator) MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey) *MKCiphertext {

	out := eval.tensor(c1.ciphertexts.Element, c2.ciphertexts.Element)

	// Call Relin on the resulting ciphertext
	RelinearizationWithSharedRelinKey(relinKey, out)

	return out
}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func (eval *mkEvaluator) MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey) *MKCiphertext {

	out := eval.tensor(c1.ciphertexts.Element, c2.ciphertexts.Element)

	// Call Relin alg 2
	RelinearizationOnTheFly(evalKeys, publicKeys, out, eval.params)

	return out
}

func (eval *mkEvaluator) modUpAndNTT(ct *ckks.Element, cQ, cQMul []*ring.Poly) {
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
func (eval *mkEvaluator) tensor(ct0, ct1 *ckks.Element) *MKCiphertext {

	nbrElements := ct0.Degree() + 1 // k+1

	outputDegree := nbrElements * nbrElements

	out := new(MKCiphertext)
	out.ciphertexts = ckks.NewCiphertext(eval.params, outputDegree-1, ct0.Level(), ct0.Scale())

	c0Q1 := make([]*ring.Poly, nbrElements)
	c0Q2 := make([]*ring.Poly, nbrElements)

	for i := uint64(0); i < nbrElements; i++ {
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

		c00Q1 := make([]*ring.Poly, nbrElements)
		c00Q2 := make([]*ring.Poly, nbrElements)

		for i := range ct0.Value() {

			c00Q1[i] = eval.ringQ.NewPoly()
			c00Q2[i] = eval.ringQMul.NewPoly()

			eval.ringQ.MForm(c0Q1[i], c00Q1[i])
			eval.ringQMul.MForm(c0Q2[i], c00Q2[i])
		}

		for i := uint64(0); i < nbrElements; i++ {
			for j := i + 1; j < nbrElements; j++ {
				eval.ringQ.MulCoeffsMontgomery(c00Q1[i], c0Q1[j], c2Q1[nbrElements*i+j])
				eval.ringQMul.MulCoeffsMontgomery(c00Q2[i], c0Q2[j], c2Q2[nbrElements*i+j])

				eval.ringQ.Add(c2Q1[i+j], c2Q1[i+j], c2Q1[nbrElements*i+j])
				eval.ringQMul.Add(c2Q2[i+j], c2Q2[i+j], c2Q2[nbrElements*i+j])
			}
		}

		for i := uint64(0); i < ct0.Degree()+1; i++ {
			eval.ringQ.MulCoeffsMontgomeryAndAdd(c00Q1[i], c0Q1[i], c2Q1[i<<1])
			eval.ringQMul.MulCoeffsMontgomeryAndAdd(c00Q2[i], c0Q2[i], c2Q2[i<<1])
		}

		// Normal case
	} else {

		c1Q1 := make([]*ring.Poly, nbrElements)
		c1Q2 := make([]*ring.Poly, nbrElements)

		for i := uint64(0); i < nbrElements; i++ {
			c1Q1[i] = eval.ringQ.NewPoly()
			c1Q2[i] = eval.ringQMul.NewPoly()
		}

		eval.modUpAndNTT(ct1, c1Q1, c1Q2)

		for i := range ct0.Value() {
			eval.ringQ.MForm(c0Q1[i], c0Q1[i])
			eval.ringQMul.MForm(c0Q2[i], c0Q2[i])
			for j := range ct1.Value() {
				eval.ringQ.MulCoeffsMontgomeryAndAdd(c0Q1[i], c1Q1[j], c2Q1[int(nbrElements)+1*i+j])
				eval.ringQMul.MulCoeffsMontgomeryAndAdd(c0Q2[i], c1Q2[j], c2Q2[int(nbrElements)*i+j])
			}
		}
	}

	return out
}

// Rescale takes a ciphertext at level l reduces it until it reaches its original scale
// this function is the same as in ckks/evaluator.go
func (eval *mkEvaluator) Rescale(c *MKCiphertext, out *MKCiphertext) {

	eval.ckksEval.Rescale(c.ciphertexts, eval.params.Scale(), c.ciphertexts)
}
