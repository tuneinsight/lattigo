package mkbfv

import (
	"math/big"
	"sort"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKEvaluator is a wrapper for the bfv evaluator
type MKEvaluator interface {
	Add(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	Sub(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	AddPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext
	SubPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext
	Neg(c *MKCiphertext) *MKCiphertext
	MultPlaintext(pt *bfv.PlaintextMul, c *MKCiphertext) *MKCiphertext
	Mul(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext
	RelinInPlace(ct *MKCiphertext, evalKeys []*mkrlwe.MKEvaluationKey, publicKeys []*mkrlwe.MKPublicKey)
	Rotate(c *MKCiphertext, n int, keys []*mkrlwe.MKEvalGalKey) *MKCiphertext
	TensorAndRescale(ct0, ct1 *bfv.Ciphertext) *MKCiphertext
	NewPlaintextFromValue([]uint64) *bfv.Plaintext
	NewPlaintextMulFromValue([]uint64) *bfv.PlaintextMul
}

type mkEvaluator struct {
	bfvEval         bfv.Evaluator
	params          *bfv.Parameters
	ringQ           *ring.Ring
	ringP           *ring.Ring
	ringQMul        *ring.Ring
	pHalf           *big.Int
	samplerGaussian *ring.GaussianSampler
	polyPoolQ1      []*ring.Poly
	polyPoolQ2      []*ring.Poly
	convertorQQMul  *ring.FastBasisExtender
	convertorQP     *ring.FastBasisExtender
	encoder         bfv.Encoder
}

// NewMKEvaluator returns an evaluator for the multi key bfv scheme.
func NewMKEvaluator(params *bfv.Parameters) MKEvaluator {

	if params == nil {
		panic("Cannot create evaluator with uninitilized parameters")
	}

	ringQ := mkrlwe.GetRingQ(&params.Parameters)
	ringP := mkrlwe.GetRingP(&params.Parameters)
	ringQMul := mkrlwe.GetRingQMul(&params.Parameters)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	sampler := mkrlwe.GetGaussianSampler(&params.Parameters, ringQ, prng)
	convertorQQMul := ring.NewFastBasisExtender(ringQ, ringQMul)
	convertorQP := ring.NewFastBasisExtender(ringQ, ringP)

	pHalf := new(big.Int).Rsh(ringQMul.ModulusBigint, 1)

	return &mkEvaluator{
		bfvEval:         bfv.NewEvaluator(*params, rlwe.EvaluationKey{}),
		params:          params,
		ringQ:           ringQ,
		ringP:           ringP,
		ringQMul:        ringQMul,
		pHalf:           pHalf,
		samplerGaussian: sampler,
		convertorQQMul:  convertorQQMul,
		convertorQP:     convertorQP,
		encoder:         bfv.NewEncoder(*params)}
}

// Add adds the ciphertexts component wise and expend their list of involved peers. Returns a new ciphertext
func (eval *mkEvaluator) Add(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := NewMKCiphertext(padded1.PeerID, eval.ringQ, eval.params)

	out.Ciphertexts = eval.AddNew(padded1.Ciphertexts, padded2.Ciphertexts)

	return out
}

// AddNew adds op0 to op1 and creates a new element ctOut to store the result.
func (eval *mkEvaluator) AddNew(op0, op1 bfv.Operand) (ctOut *bfv.Ciphertext) {
	ctOut = bfv.NewCiphertext(*eval.params, utils.MaxUint64(op0.Degree(), op1.Degree()))
	eval.AddWithNil(op0, op1, ctOut)
	return
}

// Add adds op0 to op1 and returns the result in ctOut.
func (eval *mkEvaluator) AddWithNil(op0, op1 bfv.Operand, ctOut *bfv.Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))
	eval.evaluateInPlaceBinary(el0, el1, elOut, false, eval.ringQ.Add)
}

func (eval *mkEvaluator) getElemAndCheckBinary(op0, op1, opOut bfv.Operand, opOutMinDegree uint64) (el0, el1, elOut *rlwe.Element) {

	if op0 == nil && op1 == nil && opOut == nil {
		el0, el1, elOut = nil, nil, nil
		return
	}

	if op0 == nil && op1 == nil && opOut != nil {
		el0, el1, elOut = nil, nil, opOut.El()
		return
	}

	if op0 == nil && op1 != nil && opOut == nil {
		el0, el1, elOut = nil, op1.El(), nil
		return
	}

	if op0 == nil && op1 != nil && opOut != nil {
		el0, el1, elOut = nil, op1.El(), opOut.El()
		return
	}

	if op0 != nil && op1 == nil && opOut == nil {
		el0, el1, elOut = op0.El(), nil, nil
		return
	}

	if op0 != nil && op1 == nil && opOut != nil {
		el0, el1, elOut = op0.El(), nil, opOut.El()
		return
	}

	if op0 != nil && op1 != nil && opOut == nil {
		el0, el1, elOut = op0.El(), op1.El(), nil
		return
	}

	if op0 != nil && op1 != nil && opOut != nil {
		el0, el1, elOut = op0.El(), op1.El(), opOut.El()
		return
	}

	return
}

// Sub substracts the ciphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Sub(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	if c1 == nil || c2 == nil || c1.Ciphertexts == nil || c2.Ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := NewMKCiphertext(padded1.PeerID, eval.ringQ, eval.params)

	out.Ciphertexts = eval.SubNew(padded1.Ciphertexts, padded2.Ciphertexts)

	return out
}

// SubNew subtracts op1 from op0 and creates a new element ctOut to store the result.
func (eval *mkEvaluator) SubNew(op0, op1 bfv.Operand) (ctOut *bfv.Ciphertext) {
	ctOut = bfv.NewCiphertext(*eval.params, utils.MaxUint64(op0.Degree(), op1.Degree()))
	eval.SubWithNil(op0, op1, ctOut)
	return
}

// AddPlaintext adds the paintext to the ciphertexts component wise
func (eval *mkEvaluator) AddPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext {

	if c == nil || pt == nil || c.Ciphertexts == nil || pt.Value == nil {
		panic("Uninitialized inputs")
	}

	if pt.Degree() != 0 {
		panic("Plaintext must have degree 0")
	}

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params)
	val := make([]*ring.Poly, len(c.PeerID)+1)

	// copy values
	for i := uint64(1); i < uint64(len(c.PeerID)+1); i++ {
		val[i] = c.Ciphertexts.Value[i].CopyNew()
	}

	// add the plaintext value in c0
	val[0] = eval.ringQ.NewPoly()
	eval.ringQ.Add(c.Ciphertexts.Value[0], pt.Value[0], val[0])

	out.Ciphertexts.SetValue(val)

	return out
}

// SubPlaintext subtracts the plaintext to the ciphertext component wise
func (eval *mkEvaluator) SubPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext {

	if c == nil || pt == nil || c.Ciphertexts == nil || pt.Value == nil {
		panic("Uninitialized inputs")
	}

	if pt.Degree() != 0 {
		panic("Plaintext must have degree 0")
	}

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params)
	val := make([]*ring.Poly, len(c.PeerID)+1)

	// copy values
	for i := uint64(1); i < uint64(len(c.PeerID)+1); i++ {
		val[i] = c.Ciphertexts.Value[i].CopyNew()
	}

	// subtract the plaintext value to c0
	val[0] = eval.ringQ.NewPoly()
	eval.ringQ.Sub(c.Ciphertexts.Value[0], pt.Value[0], val[0])

	out.Ciphertexts.SetValue(val)

	return out
}

// Neg returns the additive inverse of a cyphertext
func (eval *mkEvaluator) Neg(c *MKCiphertext) *MKCiphertext {

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params)

	out.Ciphertexts = eval.bfvEval.NegNew(c.Ciphertexts)

	return out
}

// MultPlaintext multiplies a plaintext and a ciphertext
func (eval *mkEvaluator) MultPlaintext(pt *bfv.PlaintextMul, c *MKCiphertext) *MKCiphertext {

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params)

	out.Ciphertexts = eval.bfvEval.MulNew(c.Ciphertexts, pt)

	return out
}

// Mul will compute the homomorphic multiplication. No relinearization is done.
func (eval *mkEvaluator) Mul(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	if c1 == c2 { // squaring case
		out := eval.TensorAndRescale(c1.Ciphertexts, c1.Ciphertexts)
		out.PeerID = c1.PeerID
		return out
	}

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := eval.TensorAndRescale(padded1.Ciphertexts, padded2.Ciphertexts)
	out.PeerID = padded1.PeerID
	return out
}

// Relinearize a ciphertext after a multiplication
func (eval *mkEvaluator) RelinInPlace(ct *MKCiphertext, evalKeys []*mkrlwe.MKEvaluationKey, publicKeys []*mkrlwe.MKPublicKey) {

	sort.Slice(evalKeys, func(i, j int) bool { return evalKeys[i].PeerID < evalKeys[j].PeerID })
	sort.Slice(publicKeys, func(i, j int) bool { return publicKeys[i].PeerID < publicKeys[j].PeerID })

	checkParticipantsEvalKey(ct.PeerID, evalKeys)
	checkParticipantsPubKey(ct.PeerID, publicKeys)

	//pass ciphertext in NTT domain
	for _, v := range ct.Ciphertexts.Value {
		eval.ringQ.NTT(v, v)
	}

	mkrlwe.Relinearization(evalKeys, publicKeys, &ct.Ciphertexts.Value, &eval.params.Parameters, uint64(len(eval.ringQ.Modulus)-1))

	//pass ciphertext out of NTT domain
	for _, v := range ct.Ciphertexts.Value {
		eval.ringQ.InvNTT(v, v)
	}
}

func (eval *mkEvaluator) modUpAndNTT(ct *bfv.Ciphertext, cQ, cQMul []*ring.Poly) {
	levelQ := uint64(len(eval.ringQ.Modulus) - 1)
	for i := range ct.Value {
		if ct.Value[i] != nil {
			eval.convertorQQMul.ModUpSplitQP(levelQ, ct.Value[i], cQMul[i])
			eval.ringQ.NTTLazy(ct.Value[i], cQ[i])
		}
		if cQMul[i] != nil {
			eval.ringQMul.NTTLazy(cQMul[i], cQMul[i])
		}

	}
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

// tensor computes the tensor product between 2 ciphertexts and returns the result in out
// c1 and c2 must have be of dimension k+1, where k = #participants
// out has dimensions (k+1)**2
func (eval *mkEvaluator) TensorAndRescale(ct0, ct1 *bfv.Ciphertext) *MKCiphertext {

	nbrElements := ct0.Degree() + 1 // k+1

	outputDegree := nbrElements * nbrElements // (k+1)**2
	out := new(MKCiphertext)
	out.Ciphertexts = bfv.NewCiphertext(*eval.params, outputDegree-1)

	c0Q1 := make([]*ring.Poly, nbrElements)
	c0Q2 := make([]*ring.Poly, nbrElements)

	for i := uint64(0); i < nbrElements; i++ {
		c0Q1[i] = eval.ringQ.NewPoly()
		c0Q2[i] = eval.ringQMul.NewPoly()
	}
	eval.modUpAndNTT(ct0, c0Q1, c0Q2)        // split ct0 in ringQ and ringQMul
	c2Q1 := make([]*ring.Poly, outputDegree) // prepare output
	c2Q2 := make([]*ring.Poly, outputDegree)

	for i := uint64(0); i < outputDegree; i++ {

		c2Q1[i] = eval.ringQ.NewPoly()
		c2Q2[i] = eval.ringQMul.NewPoly()
	}

	c1Q1 := make([]*ring.Poly, nbrElements)
	c1Q2 := make([]*ring.Poly, nbrElements)

	for i := uint64(0); i < nbrElements; i++ {
		c1Q1[i] = eval.ringQ.NewPoly()
		c1Q2[i] = eval.ringQMul.NewPoly()
	}

	eval.modUpAndNTT(ct1, c1Q1, c1Q2)

	for i := range ct0.Value {

		if ct0.Value[i] != nil {
			eval.ringQ.MForm(c0Q1[i], c0Q1[i])
			eval.ringQMul.MForm(c0Q2[i], c0Q2[i])

			for j := range ct1.Value {
				if ct1.Value[j] != nil {
					eval.ringQ.MulCoeffsMontgomeryAndAdd(c0Q1[i], c1Q1[j], c2Q1[int(nbrElements)*i+j])
					eval.ringQMul.MulCoeffsMontgomeryAndAdd(c0Q2[i], c1Q2[j], c2Q2[int(nbrElements)*i+j])
				} else {
					c2Q1[int(nbrElements)*i+j].Zero()
					c2Q2[int(nbrElements)*i+j].Zero()
				}
			}
		} else {
			for j := range ct1.Value {
				c2Q1[int(nbrElements)*i+j].Zero()
				c2Q2[int(nbrElements)*i+j].Zero()
			}
		}

	}

	eval.quantize(c2Q1, c2Q2, out.Ciphertexts.Element)

	return out
}

// quantize multiplies the values of an element by t/q
func (eval *mkEvaluator) quantize(c2Q1, c2Q2 []*ring.Poly, ctOut *rlwe.Element) {

	levelQ := uint64(len(eval.ringQ.Modulus) - 1)
	levelQMul := uint64(len(eval.ringQMul.Modulus) - 1)

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.Value {

		eval.ringQ.InvNTTLazy(c2Q1[i], c2Q1[i])
		eval.ringQMul.InvNTTLazy(c2Q2[i], c2Q2[i])

		// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
		eval.convertorQQMul.ModDownSplitQP(levelQ, levelQMul, c2Q1[i], c2Q2[i], c2Q2[i])

		// Centers (ct(x)Q -> P)/Q by (P-1)/2 and extends ((ct(x)Q -> P)/Q) to the basis Q
		eval.ringQMul.AddScalarBigint(c2Q2[i], eval.pHalf, c2Q2[i])
		eval.convertorQQMul.ModUpSplitPQ(levelQMul, c2Q2[i], ctOut.Value[i])
		eval.ringQ.SubScalarBigint(ctOut.Value[i], eval.pHalf, ctOut.Value[i])

		// Option (2) (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
		eval.ringQ.MulScalar(ctOut.Value[i], eval.params.T(), ctOut.Value[i])

	}
}

// Rotate rotate the columns of the ciphertext by n to the left and return the result in a new ciphertext
func (eval *mkEvaluator) Rotate(c *MKCiphertext, n int, keys []*mkrlwe.MKEvalGalKey) *MKCiphertext {

	sort.Slice(keys, func(i, j int) bool { return keys[i].PeerID < keys[j].PeerID })

	checkParticipantsGalKey(c.PeerID, keys)

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params)

	galEl := eval.params.GaloisElementForColumnRotationBy(n)

	level := uint64(len(eval.ringQ.Modulus)) - 1

	ringQP := mkrlwe.GetRingQP(&eval.params.Parameters)
	ringQ := mkrlwe.GetRingQ(&eval.params.Parameters)

	k := uint64(len(c.PeerID))

	res := make([]*ring.Poly, k+1)

	restmpQ := make([]*ring.Poly, k+1)
	restmpP := make([]*ring.Poly, k+1)

	restmpQ[0] = eval.ringQ.NewPoly()
	restmpP[0] = eval.ringP.NewPoly()

	for i := uint64(0); i < k+1; i++ {
		res[i] = eval.ringQ.NewPoly()
	}

	// pass ciphertext in NTT
	for _, v := range c.Ciphertexts.Value {
		ringQ.NTT(v, v)
	}

	gal0Q := mkrlwe.NewDecomposedPoly(eval.ringQ, eval.params.Beta())
	gal0P := mkrlwe.NewDecomposedPoly(eval.ringP, eval.params.Beta())
	gal1Q := mkrlwe.NewDecomposedPoly(eval.ringQ, eval.params.Beta())
	gal1P := mkrlwe.NewDecomposedPoly(eval.ringP, eval.params.Beta())

	for i := uint64(1); i <= k; i++ {

		prepareGaloisEvaluationKey(i, level, eval.params.Beta(), keys, gal0Q, gal0P, gal1Q, gal1P)

		permutedCipher := eval.ringQ.NewPoly() // apply rotation to the ciphertext
		index := ring.PermuteNTTIndex(galEl, ringQP.N)
		ring.PermuteNTTWithIndexLvl(level, c.Ciphertexts.Value[i], index, permutedCipher)

		decomposedPermutedQ, decomposedPermutedP := mkrlwe.GInverse(permutedCipher, &eval.params.Parameters, uint64(len(eval.ringQ.Modulus)-1))

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
	eval.convertorQP.ModDownSplitNTTPQ(level, restmpQ[0], restmpP[0], tmpModDown)
	eval.ringQ.AddLvl(level, res[0], tmpModDown, res[0])

	// finalize computation of ci'
	for i := uint64(1); i <= k; i++ {

		eval.convertorQP.ModDownSplitNTTPQ(level, restmpQ[i], restmpP[i], tmpModDown)
		eval.ringQ.CopyLvl(level, tmpModDown, res[i])
	}

	// pass input ciphertext and result out of NTT domain
	for _, v := range res {
		ringQ.InvNTT(v, v)
	}
	for _, v := range c.Ciphertexts.Value {
		ringQ.InvNTT(v, v)
	}

	out.Ciphertexts.SetValue(res)

	return out
}

// prepare galois evaluation keys for operations in split crt basis
func prepareGaloisEvaluationKey(j, level, beta uint64, galKeys []*mkrlwe.MKEvalGalKey, gal0Q, gal0P, gal1Q, gal1P *mkrlwe.MKDecomposedPoly) {

	for u := uint64(0); u < beta; u++ {

		gal0Q.Poly[u].Coeffs = galKeys[j-1].Key[0].Poly[u].Coeffs[:level+1]

		gal0P.Poly[u].Coeffs = galKeys[j-1].Key[0].Poly[u].Coeffs[level+1:]

		gal1Q.Poly[u].Coeffs = galKeys[j-1].Key[1].Poly[u].Coeffs[:level+1]

		gal1P.Poly[u].Coeffs = galKeys[j-1].Key[1].Poly[u].Coeffs[level+1:]

	}
}

// NewPlaintextFromValue returns a plaintext in ringQ scaled by Q/t
func (eval *mkEvaluator) NewPlaintextFromValue(value []uint64) *bfv.Plaintext {

	plaintext := bfv.NewPlaintext(*eval.params)

	// Encode
	eval.encoder.EncodeUint(value, plaintext)

	return plaintext
}

// NewPlaintextMulFromValue returns a plaintext containing the provided values. This plaintext should only be used for multiplication
func (eval *mkEvaluator) NewPlaintextMulFromValue(value []uint64) *bfv.PlaintextMul {

	plaintext := bfv.NewPlaintextMul(*eval.params)

	eval.encoder.EncodeUintMul(value, plaintext)

	return plaintext
}

// evaluateInPlaceBinary applies the provided function in place on el0 and el1 and returns the result in elOut.
func (eval *mkEvaluator) evaluateInPlaceBinary(el0, el1, elOut *rlwe.Element, isSub bool, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {

	for i := uint64(0); i < el0.Degree()+1; i++ {

		if el0.Value[i] == nil && el1.Value[i] == nil {
			elOut.Value[i] = nil
		}
		if el0.Value[i] == nil && el1.Value[i] != nil {

			if !isSub {
				elOut.Value[i] = el1.Value[i]
			} else {
				eval.ringQ.Neg(el1.Value[i], elOut.Value[i])
			}
		}
		if el0.Value[i] != nil && el1.Value[i] == nil {
			elOut.Value[i] = el0.Value[i]
		}
		if el0.Value[i] != nil && el1.Value[i] != nil {
			evaluate(el0.Value[i], el1.Value[i], elOut.Value[i])
		}

	}

}

// SubWithNil subtracts op1 from op0 and returns the result in cOut.
func (eval *mkEvaluator) SubWithNil(op0, op1 bfv.Operand, ctOut *bfv.Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))
	eval.evaluateInPlaceBinary(el0, el1, elOut, true, eval.ringQ.Sub)

	if el0 != nil && el1 != nil && elOut != nil {
		if el0.Degree() < el1.Degree() {
			for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
				eval.ringQ.Neg(ctOut.Value[i], ctOut.Value[i])
			}
		}
	}
}
