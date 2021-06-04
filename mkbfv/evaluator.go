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
	RelinInPlace(ct *MKCiphertext, relinKeys []*mkrlwe.MKRelinearizationKey, publicKeys []*mkrlwe.MKPublicKey)
	Rotate(c *MKCiphertext, n int, keys []*mkrlwe.MKRotationKey) *MKCiphertext
	NewPlaintextFromValue([]uint64) *bfv.Plaintext
	NewPlaintextMulFromValue([]uint64) *bfv.PlaintextMul
	ConvertToMKCiphertext(ct []*bfv.Ciphertext, ids []uint64) []*MKCiphertext
	ConvertToBFVCiphertext(mkCT *MKCiphertext) []*bfv.Ciphertext
}

type mkEvaluator struct {
	bfvEval         bfv.Evaluator
	params          *bfv.Parameters
	ringQ           *ring.Ring
	ringP           *ring.Ring
	ringQMul        *ring.Ring
	ringQP          *ring.Ring
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
	ringQP := mkrlwe.GetRingQP(&params.Parameters)

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
		ringQP:          ringQP,
		pHalf:           pHalf,
		samplerGaussian: sampler,
		convertorQQMul:  convertorQQMul,
		convertorQP:     convertorQP,
		encoder:         bfv.NewEncoder(*params)}
}

// ConvertToMKCiphertext takes a slice of bfv ciphertexts and their ID and convert them to multi key ciphertexts
// ciphers will be ordered with respect to their IDs (ascending order)
func (eval *mkEvaluator) ConvertToMKCiphertext(ct []*bfv.Ciphertext, ids []uint64) []*MKCiphertext {

	if len(ids) != len(ct) {
		panic("ids and ciphertexts must be of ssame length")
	}

	res := make([]*MKCiphertext, len(ct))

	for i, c := range ct {

		if c.Degree() != 1 {
			panic("Cannot convert ciphertexts of degree different than 1")
		}

		newCipher := new(MKCiphertext)
		newCipher.Ciphertexts = c
		newCipher.PeerID = []uint64{ids[i]}
		res[i] = newCipher
	}

	return res
}

// ConvertToBFVCiphertext transforms and MKCiphertext into bfv Ciphertexts containing c0 and a ciphertext part.
// ciphers are outputed in ascending id order
func (eval *mkEvaluator) ConvertToBFVCiphertext(mkCT *MKCiphertext) []*bfv.Ciphertext {

	res := make([]*bfv.Ciphertext, len(mkCT.PeerID))

	c0 := mkCT.Ciphertexts.Value[0]

	for i, v := range mkCT.Ciphertexts.Value {

		if i != 0 {
			newCipher := new(bfv.Ciphertext)
			newCipher.Element = new(rlwe.Element)
			newCipher.Value = []*ring.Poly{c0, v}
			res[i-1] = newCipher
		}
	}

	return res
}

// Add adds the ciphertexts component wise and expend their list of involved peers. Returns a new ciphertext
func (eval *mkEvaluator) Add(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	if c1 == nil || c2 == nil || c1.Ciphertexts == nil || c2.Ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := NewMKCiphertext(padded1.PeerID, eval.ringQ, eval.params)
	eval.evaluateInPlaceBinary(padded1.Ciphertexts, padded2.Ciphertexts, out.Ciphertexts, false, eval.ringQ.Add)

	return out
}

// Sub substracts the ciphertexts component wise and expend their list of involved peers. Returns a new ciphertext
func (eval *mkEvaluator) Sub(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	if c1 == nil || c2 == nil || c1.Ciphertexts == nil || c2.Ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := NewMKCiphertext(padded1.PeerID, eval.ringQ, eval.params)
	eval.evaluateInPlaceBinary(padded1.Ciphertexts, padded2.Ciphertexts, out.Ciphertexts, true, eval.ringQ.Sub)

	return out
}

// AddPlaintext adds the paintext to the ciphertexts component wise
func (eval *mkEvaluator) AddPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext {

	out := new(MKCiphertext)
	out.PeerID = c.PeerID

	out.Ciphertexts = eval.bfvEval.AddNew(c.Ciphertexts, pt)

	return out
}

// SubPlaintext subtracts the plaintext to the ciphertext component wise
func (eval *mkEvaluator) SubPlaintext(pt *bfv.Plaintext, c *MKCiphertext) *MKCiphertext {

	out := new(MKCiphertext)
	out.PeerID = c.PeerID

	out.Ciphertexts = eval.bfvEval.SubNew(c.Ciphertexts, pt)

	return out
}

// Neg returns the additive inverse of a cyphertext
func (eval *mkEvaluator) Neg(c *MKCiphertext) *MKCiphertext {

	out := new(MKCiphertext)
	out.PeerID = c.PeerID

	out.Ciphertexts = eval.bfvEval.NegNew(c.Ciphertexts)

	return out
}

// MultPlaintext multiplies a plaintext and a ciphertext
func (eval *mkEvaluator) MultPlaintext(pt *bfv.PlaintextMul, c *MKCiphertext) *MKCiphertext {

	out := new(MKCiphertext)
	out.PeerID = c.PeerID

	out.Ciphertexts = eval.bfvEval.MulNew(c.Ciphertexts, pt)

	return out
}

// Mul will compute the homomorphic multiplication. No relinearization is done.
func (eval *mkEvaluator) Mul(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	if c1 == c2 { // squaring case
		out := eval.tensorAndRescale(c1.Ciphertexts, c1.Ciphertexts)
		out.PeerID = c1.PeerID
		return out
	}

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	out := eval.tensorAndRescale(padded1.Ciphertexts, padded2.Ciphertexts)
	out.PeerID = padded1.PeerID

	return out
}

// Relinearize a ciphertext after a multiplication
func (eval *mkEvaluator) RelinInPlace(ct *MKCiphertext, relinKeys []*mkrlwe.MKRelinearizationKey, publicKeys []*mkrlwe.MKPublicKey) {

	sort.Slice(relinKeys, func(i, j int) bool { return relinKeys[i].PeerID < relinKeys[j].PeerID })
	sort.Slice(publicKeys, func(i, j int) bool { return publicKeys[i].PeerID < publicKeys[j].PeerID })

	checkParticipantsEvalKey(ct.PeerID, relinKeys)
	checkParticipantsPubKey(ct.PeerID, publicKeys)

	//pass ciphertext in NTT domain
	for _, v := range ct.Ciphertexts.Value {
		if v != nil {
			eval.ringQ.NTT(v, v)
		}
	}

	mkrlwe.Relinearization(relinKeys, publicKeys, &ct.Ciphertexts.Value, &eval.params.Parameters, uint64(len(eval.ringQ.Modulus)-1))

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

func checkParticipantsEvalKey(peerID []uint64, relinKeys []*mkrlwe.MKRelinearizationKey) {

	for i, id := range peerID {
		if id != relinKeys[i].PeerID {
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

func checkParticipantsGalKey(peerID []uint64, rotKeys []*mkrlwe.MKRotationKey) {

	for i, id := range peerID {
		if id != rotKeys[i].PeerID {
			panic("Incorrect galois evaluation keys for the given ciphertexts")
		}
	}
}

// tensorAndRescale computes the tensor product between 2 ciphertexts, scale it by t/q and returns the result in out
// c1 and c2 must have be of dimension k+1, where k = #participants
// out has dimensions (k+1)**2
func (eval *mkEvaluator) tensorAndRescale(ct0, ct1 *bfv.Ciphertext) *MKCiphertext {

	nbrElements := ct0.Degree() + 1 // k+1

	outputDegree := nbrElements * nbrElements // (k+1)**2
	out := new(MKCiphertext)
	out.Ciphertexts = new(bfv.Ciphertext)
	out.Ciphertexts.Element = new(rlwe.Element)
	out.Ciphertexts.Value = make([]*ring.Poly, outputDegree)

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
					c2Q1[int(nbrElements)*i+j] = nil
					c2Q2[int(nbrElements)*i+j] = nil
				}
			}
		} else {
			for j := range ct1.Value {
				c2Q1[int(nbrElements)*i+j] = nil
				c2Q2[int(nbrElements)*i+j] = nil
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

		if c2Q1[i] != nil {

			ctOut.Value[i] = eval.ringQ.NewPoly()
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
}

// Rotate rotate the columns of the ciphertext by n to the left and return the result in a new ciphertext
func (eval *mkEvaluator) Rotate(c *MKCiphertext, n int, keys []*mkrlwe.MKRotationKey) *MKCiphertext {

	sort.Slice(keys, func(i, j int) bool { return keys[i].PeerID < keys[j].PeerID })

	checkParticipantsGalKey(c.PeerID, keys)

	out := NewMKCiphertext(c.PeerID, eval.ringQ, eval.params)

	galEl := eval.params.GaloisElementForColumnRotationBy(n)

	level := uint64(len(eval.ringQ.Modulus)) - 1

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
		eval.ringQ.NTT(v, v)
	}

	gal0QP := rlwe.NewSwitchingKey(eval.params.Parameters)
	gal1QP := rlwe.NewSwitchingKey(eval.params.Parameters)

	for i := uint64(1); i <= k; i++ {

		prepareRotationKey(i, level, eval.params.Beta(), keys, gal0QP, gal1QP)

		permutedCipher := eval.ringQ.NewPoly() // apply rotation to the ciphertext
		index := ring.PermuteNTTIndex(galEl, eval.ringQP.N)
		ring.PermuteNTTWithIndexLvl(level, c.Ciphertexts.Value[i], index, permutedCipher)

		decomposedPermutedQP := mkrlwe.GInverseKeySwitch(permutedCipher, &eval.params.Parameters, uint64(len(eval.ringQ.Modulus)-1))
		res0Q, res0P := mkrlwe.DotSwk(level, decomposedPermutedQP, gal0QP, eval.ringQ, eval.ringP, eval.params.Beta())

		eval.ringP.Add(restmpP[0], res0P, restmpP[0])
		eval.ringQ.AddLvl(level, restmpQ[0], res0Q, restmpQ[0])

		restmpQ[i], restmpP[i] = mkrlwe.DotSwk(level, decomposedPermutedQP, gal1QP, eval.ringQ, eval.ringP, eval.params.Beta())
	}

	// finalize computation of c0'
	index := ring.PermuteNTTIndex(galEl, eval.ringQP.N)
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
		eval.ringQ.InvNTT(v, v)
	}
	for _, v := range c.Ciphertexts.Value {
		eval.ringQ.InvNTT(v, v)
	}

	out.Ciphertexts.SetValue(res)

	return out
}

// prepare rotation keys for operations in split crt basis
func prepareRotationKey(j, level, beta uint64, rotKeys []*mkrlwe.MKRotationKey, gal0QP, gal1QP *rlwe.SwitchingKey) {

	for u := uint64(0); u < beta; u++ {
		gal0QP.Value[u][0].Coeffs = rotKeys[j-1].Key.Value[u][0].Coeffs[:level+1]
		gal0QP.Value[u][1].Coeffs = rotKeys[j-1].Key.Value[u][0].Coeffs[level+1:]
		gal1QP.Value[u][0].Coeffs = rotKeys[j-1].Key.Value[u][1].Coeffs[:level+1]
		gal1QP.Value[u][1].Coeffs = rotKeys[j-1].Key.Value[u][1].Coeffs[level+1:]

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
func (eval *mkEvaluator) evaluateInPlaceBinary(el0, el1, elOut *bfv.Ciphertext, isSub bool, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {

	for i := uint64(0); i < el0.Degree()+1; i++ {

		if el0.Value[i] == nil && el1.Value[i] == nil {
			elOut.Value[i] = nil
		} else if el0.Value[i] == nil && el1.Value[i] != nil {

			if !isSub {
				elOut.Value[i] = el1.Value[i]
			} else {
				eval.ringQ.Neg(el1.Value[i], elOut.Value[i])
			}
		} else if el0.Value[i] != nil && el1.Value[i] == nil {
			elOut.Value[i] = el0.Value[i]
		} else {
			evaluate(el0.Value[i], el1.Value[i], elOut.Value[i])
		}

	}

}
