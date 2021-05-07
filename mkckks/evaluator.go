package mkckks

import (
	"sort"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKEvaluator is a wrapper for the ckks evaluator
type MKEvaluator interface {
	Add(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	Sub(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	AddPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext
	SubPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext
	Neg(c *MKCiphertext) *MKCiphertext
	MultPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext
	Mul(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	RelinInPlace(ct *MKCiphertext, evalKeys []*mkrlwe.MKEvaluationKey, publicKeys []*mkrlwe.MKPublicKey)
	Rescale(c *MKCiphertext, out *MKCiphertext)
	Rotate(c *MKCiphertext, n int, keys []*mkrlwe.MKEvalGalKey) *MKCiphertext
	NewPlaintextFromValue([]complex128) *ckks.Plaintext
	DropLevel(ct *MKCiphertext, levels uint64)
}

type mkEvaluator struct {
	ckksEval        ckks.Evaluator
	params          *ckks.Parameters
	ringQ           *ring.Ring
	ringP           *ring.Ring
	samplerGaussian *ring.GaussianSampler
	convertor       *ring.FastBasisExtender
	encoder         ckks.Encoder
}

// NewMKEvaluator returns an evaluator for the multi key ckks scheme.
func NewMKEvaluator(params *ckks.Parameters) MKEvaluator {

	if params == nil {
		panic("Cannot create evaluator with uninitilized parameters")
	}

	ringQ := mkrlwe.GetRingQ(&params.Parameters)
	ringP := mkrlwe.GetRingP(&params.Parameters)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	sampler := mkrlwe.GetGaussianSampler(&params.Parameters, ringQ, prng)
	convertor := ring.NewFastBasisExtender(ringQ, ringP)

	return &mkEvaluator{
		ckksEval:        ckks.NewEvaluator(*params, rlwe.EvaluationKey{}),
		params:          params,
		ringQ:           ringQ,
		ringP:           ringP,
		samplerGaussian: sampler,
		convertor:       convertor,
		encoder:         ckks.NewEncoder(*params),
	}
}

// Add adds the ciphertexts component wise and expend their list of involved peers. A new ciphertext is returned
func (eval *mkEvaluator) Add(c0 *MKCiphertext, c1 *MKCiphertext) *MKCiphertext {

	if c0 == nil || c1 == nil || c0.Ciphertexts == nil || c1.Ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c0, c1, eval.params)

	out := NewMKCiphertext(padded1.PeerID, eval.ringQ, eval.params, padded1.Ciphertexts.Level(), c0.Ciphertexts.Scale())

	out.Ciphertexts = eval.ckksEval.AddNew(padded1.Ciphertexts, padded2.Ciphertexts)
	return out
}

// Sub returns the component wise substraction of 2 ciphertexts
func (eval *mkEvaluator) Sub(c0 *MKCiphertext, c1 *MKCiphertext) *MKCiphertext {

	if c0 == nil || c1 == nil || c0.Ciphertexts == nil || c1.Ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c0, c1, eval.params)

	out := NewMKCiphertext(padded1.PeerID, eval.ringQ, eval.params, padded1.Ciphertexts.Level(), padded1.Ciphertexts.Scale())

	out.Ciphertexts = eval.ckksEval.SubNew(padded1.Ciphertexts, padded2.Ciphertexts)
	return out
}

// AddPlaintext adds the paintext to the ciphertexts component wise
func (eval *mkEvaluator) AddPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext {

	if c == nil || pt == nil || c.Ciphertexts == nil || pt.Value == nil {
		panic("Uninitialized inputs")
	}

	if pt.Degree() != 0 {
		panic("Plaintext must have degree 0")
	}

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params, c.Ciphertexts.Level(), c.Ciphertexts.Scale())

	out.Ciphertexts = eval.ckksEval.AddNew(c.Ciphertexts, pt)

	return out
}

// SubPlaintext subtracts the plaintext from the ciphertext component wise
func (eval *mkEvaluator) SubPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext {

	if c == nil || pt == nil || c.Ciphertexts == nil || pt.Value == nil {
		panic("Uninitialized inputs")
	}

	if pt.Degree() != 0 {
		panic("Plaintext must have degree 0")
	}

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params, c.Ciphertexts.Level(), c.Ciphertexts.Scale())

	out.Ciphertexts = eval.ckksEval.SubNew(c.Ciphertexts, pt)

	return out
}

// Neg returns the additive inverse of a cyphertext
func (eval *mkEvaluator) Neg(c *MKCiphertext) *MKCiphertext {

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params, c.Ciphertexts.Level(), c.Ciphertexts.Scale())

	out.Ciphertexts = eval.ckksEval.NegNew(c.Ciphertexts)

	return out
}

// MultPlaintext multiplies a plaintext and a ciphertext
func (eval *mkEvaluator) MultPlaintext(pt *ckks.Plaintext, c *MKCiphertext) *MKCiphertext {

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params, c.Ciphertexts.Level(), c.Ciphertexts.Scale())

	out.Ciphertexts.SetScale(pt.Scale() * c.Ciphertexts.Scale())

	val := make([]*ring.Poly, len(c.PeerID)+1)

	level := utils.MinUint64(c.Ciphertexts.Level(), pt.Level())

	tmp := eval.ringQ.NewPoly()
	eval.ringQ.MFormLvl(level, pt.Value[0], tmp)

	for i, v := range c.Ciphertexts.Value {
		val[i] = eval.ringQ.NewPoly()
		eval.ringQ.MulCoeffsMontgomeryLvl(level, tmp, v, val[i])
	}

	out.Ciphertexts.SetValue(val)

	return out
}

// Mul will compute the tensor product and output a ciphertext with degree (k+1)**2
func (eval *mkEvaluator) Mul(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	nbrElements := padded1.Ciphertexts.Degree() + 1 // k+1

	outputDegree := nbrElements * nbrElements // (k+1)**2

	el1 := padded1.Ciphertexts
	el2 := padded2.Ciphertexts
	level := utils.MinUint64(el1.Level(), el2.Level())

	out := new(MKCiphertext)
	out.Ciphertexts = ckks.NewCiphertext(*eval.params, outputDegree-1, level, el1.Scale()*el2.Scale())
	out.PeerID = padded1.PeerID

	if !el1.IsNTT {
		panic("cannot MulRelin: op0 must be in NTT")
	}

	if !el2.IsNTT {
		panic("cannot MulRelin: op1 must be in NTT")
	}

	ringQ := eval.ringQ

	tmp1 := ringQ.NewPoly()

	resCipher := out.Ciphertexts.Value

	for i, v1 := range el1.Value {

		ringQ.MFormLvl(level, v1, tmp1)

		for j, v2 := range el2.Value {

			index := int(nbrElements)*i + j
			ringQ.MulCoeffsMontgomeryLvl(level, tmp1, v2, resCipher[index])
		}
	}

	return out
}

// Relinearize a ciphertext after a multiplication
func (eval *mkEvaluator) RelinInPlace(ct *MKCiphertext, evalKeys []*mkrlwe.MKEvaluationKey, publicKeys []*mkrlwe.MKPublicKey) {

	sort.Slice(evalKeys, func(i, j int) bool { return evalKeys[i].PeerID < evalKeys[j].PeerID })
	sort.Slice(publicKeys, func(i, j int) bool { return publicKeys[i].PeerID < publicKeys[j].PeerID })

	checkParticipantsEvalKey(ct.PeerID, evalKeys)
	checkParticipantsPubKey(ct.PeerID, publicKeys)

	mkrlwe.Relinearization(evalKeys, publicKeys, &ct.Ciphertexts.Value, &eval.params.Parameters, ct.Ciphertexts.Level())
}

// Rescale takes a ciphertext at level l reduces it until it reaches its original
// this function is the same as in ckks/evaluator.go
func (eval *mkEvaluator) Rescale(c *MKCiphertext, out *MKCiphertext) {

	eval.ckksEval.Rescale(c.Ciphertexts, eval.params.Scale(), c.Ciphertexts)
}

// DropLevel drops the level of the given ciphertext by levels. No rescaling is applied
func (eval *mkEvaluator) DropLevel(ct *MKCiphertext, levels uint64) {

	eval.ckksEval.DropLevel(ct.Ciphertexts, levels)
}

// Rotate rotate the columns of the ciphertext by n to the left and return the result in a new ciphertext
func (eval *mkEvaluator) Rotate(c *MKCiphertext, n int, keys []*mkrlwe.MKEvalGalKey) *MKCiphertext {

	sort.Slice(keys, func(i, j int) bool { return keys[i].PeerID < keys[j].PeerID })

	checkParticipantsGalKey(c.PeerID, keys)

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params, c.Ciphertexts.Level(), c.Ciphertexts.Scale())

	galEl := eval.params.GaloisElementForColumnRotationBy(n)

	level := c.Ciphertexts.Level()

	ringQP := mkrlwe.GetRingQP(&eval.params.Parameters)

	k := uint64(len(c.PeerID))

	res := make([]*ring.Poly, k+1)

	restmpQ := make([]*ring.Poly, k+1)
	restmpP := make([]*ring.Poly, k+1)
	restmpQ[0] = eval.ringQ.NewPoly()
	restmpP[0] = eval.ringP.NewPoly()

	for i := uint64(0); i < k+1; i++ {
		res[i] = eval.ringQ.NewPoly()
	}

	gal0Q := mkrlwe.NewDecomposedPoly(eval.ringQ, eval.params.Beta())
	gal0P := mkrlwe.NewDecomposedPoly(eval.ringP, eval.params.Beta())
	gal1Q := mkrlwe.NewDecomposedPoly(eval.ringQ, eval.params.Beta())
	gal1P := mkrlwe.NewDecomposedPoly(eval.ringP, eval.params.Beta())

	for i := uint64(1); i <= k; i++ {

		prepareGaloisEvaluationKey(i, level, uint64(len(eval.ringQ.Modulus)), eval.params.Beta(), keys, gal0Q, gal0P, gal1Q, gal1P)

		permutedCipher := eval.ringQ.NewPoly() // apply rotation to the ciphertext
		index := ring.PermuteNTTIndex(galEl, ringQP.N)
		ring.PermuteNTTWithIndexLvl(level, c.Ciphertexts.Value[i], index, permutedCipher)

		decomposedPermutedQ, decomposedPermutedP := mkrlwe.GInverse(permutedCipher, &eval.params.Parameters, level)

		res0P := mkrlwe.Dot(decomposedPermutedP, gal0P, eval.ringP) // dot product and add in c0''
		res0Q := mkrlwe.DotLvl(level, decomposedPermutedQ, gal0Q, eval.ringQ)

		eval.ringP.Add(restmpP[0], res0P, restmpP[0])
		eval.ringQ.AddLvl(level, restmpQ[0], res0Q, restmpQ[0])

		restmpP[i] = mkrlwe.Dot(decomposedPermutedP, gal1P, eval.ringP) // dot product and put in ci''
		restmpQ[i] = mkrlwe.DotLvl(level, decomposedPermutedQ, gal1Q, eval.ringQ)

	}

	// finalize computation of c0'
	index := ring.PermuteNTTIndex(galEl, ringQP.N)
	ring.PermuteNTTWithIndexLvl(level, c.Ciphertexts.Value[0], index, res[0])

	tmpModDown := eval.ringQ.NewPoly()
	eval.convertor.ModDownSplitNTTPQ(level, restmpQ[0], restmpP[0], tmpModDown)
	eval.ringQ.AddLvl(level, res[0], tmpModDown, res[0])

	// finalize computation of ci'
	for i := uint64(1); i <= k; i++ {

		eval.convertor.ModDownSplitNTTPQ(level, restmpQ[i], restmpP[i], tmpModDown)
		eval.ringQ.CopyLvl(level, tmpModDown, res[i])
	}

	out.Ciphertexts.SetValue(res)

	return out
}

// prepare galois evaluation keys for operations in split crt basis
func prepareGaloisEvaluationKey(j, level, modulus, beta uint64, galKeys []*mkrlwe.MKEvalGalKey, gal0Q, gal0P, gal1Q, gal1P *mkrlwe.MKDecomposedPoly) {

	for u := uint64(0); u < beta; u++ {

		gal0Q.Poly[u].Coeffs = galKeys[j-1].Key[0].Poly[u].Coeffs[:level+1]

		gal0P.Poly[u].Coeffs = galKeys[j-1].Key[0].Poly[u].Coeffs[modulus:]

		gal1Q.Poly[u].Coeffs = galKeys[j-1].Key[1].Poly[u].Coeffs[:level+1]

		gal1P.Poly[u].Coeffs = galKeys[j-1].Key[1].Poly[u].Coeffs[modulus:]

	}

}

// NewPlaintextFromValue returns a plaintext from the provided values
func (eval *mkEvaluator) NewPlaintextFromValue(values []complex128) *ckks.Plaintext {

	return eval.encoder.EncodeNTTAtLvlNew(eval.params.MaxLevel(), values, eval.params.LogSlots())
}

func getCiphertextIndex(peerID uint64, ct *MKCiphertext) uint64 {

	for i, id := range ct.PeerID {

		if id == peerID {
			return uint64(i + 1)
		}
	}

	return 0
}

func checkParticipantsEvalKey(peerID []uint64, evalKeys []*mkrlwe.MKEvaluationKey) {

	for i, id := range peerID {
		if id != evalKeys[i].PeerID {
			panic("Incorrect evaluation keys for the given ciphertexts")
		}
	}
}

func checkParticipantsPubKey(peerID []uint64, pubKeys []*mkrlwe.MKPublicKey) {

	for i, id := range peerID {
		if id != pubKeys[i].PeerID {
			panic("Incorrect public keys for the given ciphertexts")
		}
	}
}

func checkParticipantsGalKey(peerID []uint64, galKeys []*mkrlwe.MKEvalGalKey) {

	for i, id := range peerID {
		if id != galKeys[i].PeerID {
			panic("Incorrect galois evaluation keys for the given ciphertexts")
		}
	}
}
