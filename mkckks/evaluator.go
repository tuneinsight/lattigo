package mkckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKEvaluator is a wrapper for the ckks evaluator
type MKEvaluator interface {
	AddNew(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	Add(c0 *MKCiphertext, c1 *MKCiphertext, cout *MKCiphertext)
	Sub(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	AddPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext
	SubPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext
	Neg(c *MKCiphertext) *MKCiphertext
	MultPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext
	MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey) *MKCiphertext
	Rescale(c *MKCiphertext, out *MKCiphertext)
	NewPlaintextFromValue([]complex128) *ckks.Plaintext
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
	encoder         ckks.Encoder
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
		convertor:       convertor,
		encoder:         ckks.NewEncoder(params),
	}
}

// AddNew adds the ciphertexts component wise and expend their list of involved peers. A new ciphertext is returned
func (eval *mkEvaluator) AddNew(c0 *MKCiphertext, c1 *MKCiphertext) *MKCiphertext {

	if c0 == nil || c1 == nil || c0.ciphertexts == nil || c1.ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c0, c1, eval.params)

	out := NewMKCiphertext(padded1.peerIDs, eval.ringQ, eval.params, c0.ciphertexts.Level())

	out.ciphertexts = eval.ckksEval.AddNew(padded1.ciphertexts, padded2.ciphertexts)
	return out
}

// Add adds the ciphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Add(c0 *MKCiphertext, c1 *MKCiphertext, cout *MKCiphertext) {

	if c0 == nil || cout == nil || c1 == nil || c0.ciphertexts == nil || c1.ciphertexts == nil || cout.ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c0, c1, eval.params)

	eval.ckksEval.Add(padded1.ciphertexts, padded2.ciphertexts, cout.ciphertexts)
}

// Sub returns the component wise substraction of 2 ciphertexts
func (eval *mkEvaluator) Sub(c0 *MKCiphertext, c1 *MKCiphertext) *MKCiphertext {

	if c0 == nil || c1 == nil || c0.ciphertexts == nil || c1.ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c0, c1, eval.params)

	out := NewMKCiphertext(padded1.peerIDs, eval.ringQ, eval.params, c0.ciphertexts.Level())

	out.ciphertexts = eval.ckksEval.SubNew(padded1.ciphertexts, padded2.ciphertexts)
	return out
}

// AddPlaintext adds the paintext to the ciphertexts component wise
func (eval *mkEvaluator) AddPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext {

	if c == nil || pt == nil || c.ciphertexts == nil || pt.Value() == nil {
		panic("Uninitialized inputs")
	}

	if pt.Degree() != 0 {
		panic("Plaintext must have degree 0")
	}

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params, c.ciphertexts.Level())

	out.ciphertexts = eval.ckksEval.AddNew(c.ciphertexts, pt)

	return out
}

// SubPlaintext subtracts the plaintext from the ciphertext component wise
func (eval *mkEvaluator) SubPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext {

	if c == nil || pt == nil || c.ciphertexts == nil || pt.Value() == nil {
		panic("Uninitialized inputs")
	}

	if pt.Degree() != 0 {
		panic("Plaintext must have degree 0")
	}

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params, c.ciphertexts.Level())

	out.ciphertexts = eval.ckksEval.SubNew(c.ciphertexts, pt)

	return out
}

// Neg returns the additive inverse of a cyphertext
func (eval *mkEvaluator) Neg(c *MKCiphertext) *MKCiphertext {

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params, c.ciphertexts.Level())

	out.ciphertexts = eval.ckksEval.NegNew(c.ciphertexts)

	return out
}

// MultPlaintext multiplies a plaintext and a ciphertext
func (eval *mkEvaluator) MultPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext {

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params, c.ciphertexts.Level())

	out.ciphertexts.SetScale(pt.Scale() * c.ciphertexts.Scale())

	val := make([]*ring.Poly, len(c.peerIDs)+1)

	level := utils.MinUint64(c.ciphertexts.Level(), pt.Level())

	tmp := eval.ringQ.NewPoly()
	eval.ringQ.MFormLvl(level, pt.Value()[0], tmp)

	for i, v := range c.ciphertexts.Value() {
		val[i] = eval.ringQ.NewPoly()
		eval.ringQ.MulCoeffsMontgomeryLvl(level, tmp, v, val[i])
	}

	out.ciphertexts.SetValue(val)

	return out
}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func (eval *mkEvaluator) MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey) *MKCiphertext {

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	nbrElements := padded1.ciphertexts.Degree() + 1 // k+1

	outputDegree := nbrElements * nbrElements // (k+1)**2

	el1 := padded1.ciphertexts.Element
	el2 := padded2.ciphertexts.Element
	level := utils.MinUint64(el1.Level(), el2.Level())

	out := new(MKCiphertext)
	out.ciphertexts = ckks.NewCiphertext(eval.params, outputDegree-1, level, el1.Scale()*el2.Scale())
	out.peerIDs = padded1.peerIDs

	if !el1.IsNTT() {
		panic("cannot MulRelinDynamic: op0 must be in NTT")
	}

	if !el2.IsNTT() {
		panic("cannot MulRelinDynamic: op1 must be in NTT")
	}

	ringQ := eval.ringQ

	tmp1 := ringQ.NewPoly()
	tmp2 := ringQ.NewPoly()

	for i, v1 := range el1.Value() {

		ringQ.MFormLvl(level, v1, tmp1)

		for j, v2 := range el2.Value() {

			ringQ.MFormLvl(level, v2, tmp2)

			ringQ.MulCoeffsMontgomeryLvl(level, tmp1, tmp2, out.ciphertexts.Ciphertext().Value()[int(nbrElements)*i+j])
		}
	}

	// Call Relin alg 2
	RelinearizationOnTheFly(evalKeys, publicKeys, out, eval.params)

	return out
}

// Rescale takes a ciphertext at level l reduces it until it reaches its original 
// this function is the same as in ckks/evaluator.go
func (eval *mkEvaluator) Rescale(c *MKCiphertext, out *MKCiphertext) {

	eval.ckksEval.Rescale(c.ciphertexts, eval.params.Scale(), c.ciphertexts)
}

// NewPlaintextFromValue returns a plaintext from the provided values
func (eval *mkEvaluator) NewPlaintextFromValue(values []complex128) *ckks.Plaintext {

	return eval.encoder.EncodeNTTAtLvlNew(eval.params.MaxLevel(), values, eval.params.LogSlots())
}
