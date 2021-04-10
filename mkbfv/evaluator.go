package mkbfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKEvaluator is a wrapper for the bfv evaluator
type MKEvaluator interface {
	AddNew(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	Add(c1 *MKCiphertext, c2 *MKCiphertext, cout *MKCiphertext)
	Sub(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	AddPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext
	SubPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext
	Neg(c *MKCiphertext) *MKCiphertext
	MultPlaintext(pt *bfv.PlaintextMul, c *MKCiphertext) *MKCiphertext
	MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey) *MKCiphertext
	MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey) *MKCiphertext
	TensorAndRescale(ct0, ct1 *bfv.Element) *MKCiphertext
	NewPlaintextFromValue([]uint64) *bfv.Plaintext
	NewPlaintextMulFromValue([]uint64) *bfv.PlaintextMul
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
	encoder         bfv.Encoder
}

// NewMKEvaluator returns an evaluator for the multi key bfv scheme.
func NewMKEvaluator(params *bfv.Parameters) MKEvaluator {

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
		bfvEval:         bfv.NewEvaluator(params, bfv.EvaluationKey{}),
		params:          params,
		ringQ:           ringQ,
		ringQMul:        ringQMul,
		pHalf:           pHalf,
		samplerGaussian: sampler,
		convertor:       convertor,
		encoder:         bfv.NewEncoder(params)}
}

// AddNew adds the ciphertexts component wise and expend their list of involved peers. Returns a new ciphertext
func (eval *mkEvaluator) AddNew(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	if c1 == nil || c2 == nil || c1.ciphertexts == nil || c2.ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := NewMKCiphertext(padded1.peerIDs, eval.ringQ, eval.params)

	out.ciphertexts = eval.bfvEval.AddNew(padded1.ciphertexts, padded2.ciphertexts)

	return out
}

// Add adds the ciphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Add(c1 *MKCiphertext, c2 *MKCiphertext, cout *MKCiphertext) {

	if c1 == nil || c2 == nil || cout == nil || c1.ciphertexts == nil || c2.ciphertexts == nil || cout.ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	eval.bfvEval.Add(padded1.ciphertexts, padded2.ciphertexts, cout.ciphertexts)
}

// Sub substracts the ciphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Sub(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	if c1 == nil || c2 == nil || c1.ciphertexts == nil || c2.ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := NewMKCiphertext(padded1.peerIDs, eval.ringQ, eval.params)

	out.ciphertexts = eval.bfvEval.SubNew(padded1.ciphertexts, padded2.ciphertexts)

	return out
}

// AddPlaintext adds the paintext to the ciphertexts component wise
func (eval *mkEvaluator) AddPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext {

	if c == nil || pt == nil || c.ciphertexts == nil || pt.Value() == nil {
		panic("Uninitialized inputs")
	}

	if pt.Degree() != 0 {
		panic("Plaintext must have degree 0")
	}

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params)
	val := make([]*ring.Poly, len(c.peerIDs)+1)

	// copy values
	for i := uint64(1); i < uint64(len(c.peerIDs)+1); i++ {
		val[i] = c.ciphertexts.Value()[i].CopyNew()
	}

	// add the plaintext value in c0
	val[0] = eval.ringQ.NewPoly()
	eval.ringQ.Add(c.ciphertexts.Value()[0], pt.Value()[0], val[0])

	out.ciphertexts.SetValue(val)

	return out
}

// SubPlaintext subtracts the plaintext to the ciphertext component wise
func (eval *mkEvaluator) SubPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext {

	if c == nil || pt == nil || c.ciphertexts == nil || pt.Value() == nil {
		panic("Uninitialized inputs")
	}

	if pt.Degree() != 0 {
		panic("Plaintext must have degree 0")
	}

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params)
	val := make([]*ring.Poly, len(c.peerIDs)+1)

	// copy values
	for i := uint64(1); i < uint64(len(c.peerIDs)+1); i++ {
		val[i] = c.ciphertexts.Value()[i].CopyNew()
	}

	// subtract the plaintext value to c0
	val[0] = eval.ringQ.NewPoly()
	eval.ringQ.Sub(c.ciphertexts.Value()[0], pt.Value()[0], val[0])

	out.ciphertexts.SetValue(val)

	return out
}

// Neg returns the additive inverse of a cyphertext
func (eval *mkEvaluator) Neg(c *MKCiphertext) *MKCiphertext {

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params)

	out.ciphertexts = eval.bfvEval.NegNew(c.ciphertexts)

	return out
}

// MultPlaintext multiplies a plaintext and a ciphertext
func (eval *mkEvaluator) MultPlaintext(pt *bfv.PlaintextMul, c *MKCiphertext) *MKCiphertext {

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params)

	out.ciphertexts = eval.bfvEval.MulNew(c.ciphertexts, pt)

	return out
}

// MultSharedRelinKey will compute the homomorphic multiplication and relinearize the resulting cyphertext using pre computed Relin key
func (eval *mkEvaluator) MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey) *MKCiphertext {

	out := eval.TensorAndRescale(c1.ciphertexts.Element, c2.ciphertexts.Element)

	// Call Relin on the resulting ciphertext
	RelinearizationWithSharedRelinKey(relinKey, out, eval.params)

	return out
}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func (eval *mkEvaluator) MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey) *MKCiphertext {

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := eval.TensorAndRescale(padded1.ciphertexts.Element, padded2.ciphertexts.Element)
	// Call Relin alg 2
	RelinearizationOnTheFly(evalKeys, publicKeys, out, eval.params)
	out.peerIDs = padded1.peerIDs
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
func (eval *mkEvaluator) TensorAndRescale(ct0, ct1 *bfv.Element) *MKCiphertext {

	nbrElements := ct0.Degree() + 1 // k+1

	outputDegree := nbrElements * nbrElements // (k+1)**2

	out := new(MKCiphertext)
	out.ciphertexts = bfv.NewCiphertext(eval.params, outputDegree-1)

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
				eval.ringQ.MulCoeffsMontgomeryAndAdd(c0Q1[i], c1Q1[j], c2Q1[int(nbrElements)*i+j])
				eval.ringQMul.MulCoeffsMontgomeryAndAdd(c0Q2[i], c1Q2[j], c2Q2[int(nbrElements)*i+j])
			}
		}
	}

	eval.quantize(c2Q1, c2Q2, out.ciphertexts.Element)
	return out
}

// quantize multiplies the values of an element by t/q
func (eval *mkEvaluator) quantize(c2Q1, c2Q2 []*ring.Poly, ctOut *bfv.Element) {

	levelQ := uint64(len(eval.ringQ.Modulus) - 1)
	levelQMul := uint64(len(eval.ringQMul.Modulus) - 1)

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.Value() { // will iterate on (k + 1)^2 values.. TODO: check corectness !!!!!
		println("i = ", i)
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

// NewPlaintextFromValue returns a plaintext in ringQ scaled by Q/t
func (eval *mkEvaluator) NewPlaintextFromValue(value []uint64) *bfv.Plaintext {

	plaintext := bfv.NewPlaintext(eval.params)

	// Encode
	eval.encoder.EncodeUint(value, plaintext.Plaintext())

	return plaintext.Plaintext()
}

// NewPlaintextMulFromValue returns a plaintext containing the provided values. This plaintext should only be used for multiplication
func (eval *mkEvaluator) NewPlaintextMulFromValue(value []uint64) *bfv.PlaintextMul {

	plaintext := bfv.NewPlaintextMul(eval.params)

	eval.encoder.EncodeUintMul(value, plaintext)

	return plaintext
}
