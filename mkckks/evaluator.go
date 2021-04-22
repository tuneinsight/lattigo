package mkckks

import (
	"math/big"
	"sort"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
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
	RotateNew(c *MKCiphertext, n int, keys []*MKEvalGalKey) *MKCiphertext
	SwitchKeysNew(ct *MKCiphertext, switchingKey *rlwe.SwitchingKey) (ctOut *MKCiphertext)
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

	out := NewMKCiphertext(padded1.peerIDs, eval.ringQ, eval.params, padded1.ciphertexts.Level())

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

	out := NewMKCiphertext(padded1.peerIDs, eval.ringQ, eval.params, padded1.ciphertexts.Level())

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

	sort.Slice(evalKeys, func(i, j int) bool { return evalKeys[i].peerID < evalKeys[j].peerID })
	sort.Slice(publicKeys, func(i, j int) bool { return publicKeys[i].peerID < publicKeys[j].peerID })

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

// RotateNew rotate the columns of the ciphertext by n to the left and return the result in a new ciphertext
func (eval *mkEvaluator) RotateNew(c *MKCiphertext, n int, keys []*MKEvalGalKey) *MKCiphertext {

	sort.Slice(keys, func(i, j int) bool { return keys[i].peerID < keys[j].peerID })

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params, c.ciphertexts.Level())

	galEl := eval.params.GaloisElementForColumnRotationBy(n)

	level := c.ciphertexts.Level()

	ringP := GetRingP(eval.params)
	ringQP := GetRingQP(eval.params)

	k := uint64(len(c.peerIDs))

	res := make([]*ring.Poly, k+1)

	restmpQ := make([]*ring.Poly, k+1)
	restmpP := make([]*ring.Poly, k+1)

	for i := uint64(0); i < k+1; i++ {
		restmpQ[i] = eval.ringQ.NewPoly()
		restmpP[i] = ringP.NewPoly()
		res[i] = eval.ringQ.NewPoly()
	}

	for i := uint64(1); i <= k; i++ {

		gal0Q, gal0P, gal1Q, gal1P := prepareGaloisEvaluationKey(i, level, uint64(len(eval.ringQ.Modulus)), eval.params.Beta(), keys)

		permutedCipher := eval.ringQ.NewPoly()

		index := ring.PermuteNTTIndex(galEl, ringQP.N)
		ring.PermuteNTTWithIndexLvl(level, c.ciphertexts.Value()[i], index, permutedCipher)

		decomposedPermutedQ, decomposedPermutedP := GInverse(permutedCipher, eval.params, level)

		res0P := Dot(decomposedPermutedP, gal0P, ringP)
		res0Q := DotLvl(level, decomposedPermutedQ, gal0Q, eval.ringQ)

		ringP.Add(restmpP[0], res0P, restmpP[0])
		eval.ringQ.AddLvl(level, restmpQ[0], res0Q, restmpQ[0])

		restmpP[i] = Dot(decomposedPermutedP, gal1P, ringP)
		restmpQ[i] = DotLvl(level, decomposedPermutedQ, gal1Q, eval.ringQ)
	}

	// finalize computation of c0'
	index := ring.PermuteNTTIndex(galEl, ringQP.N)
	ring.PermuteNTTWithIndexLvl(level, c.ciphertexts.Value()[0], index, res[0])

	tmpModDown := eval.ringQ.NewPoly()
	eval.convertor.ModDownSplitNTTPQ(level, restmpQ[0], restmpP[0], tmpModDown)
	eval.ringQ.AddLvl(level, res[0], tmpModDown, res[0])

	// finalize computation of ci'
	for i := uint64(1); i <= k; i++ {

		eval.convertor.ModDownSplitNTTPQ(level, restmpQ[i], restmpP[i], tmpModDown)
		eval.ringQ.AddLvl(level, res[i], tmpModDown, res[i])
	}

	out.ciphertexts.SetValue(res)

	return out
}

// prepare galois evaluation keys for operations in split crt basis
func prepareGaloisEvaluationKey(j, level, modulus, beta uint64, galKeys []*MKEvalGalKey) (gal0Q, gal0P, gal1Q, gal1P *MKDecomposedPoly) {

	gal0Q = new(MKDecomposedPoly)
	gal0Q.poly = make([]*ring.Poly, beta)
	gal0P = new(MKDecomposedPoly)
	gal0P.poly = make([]*ring.Poly, beta)

	gal1Q = new(MKDecomposedPoly)
	gal1Q.poly = make([]*ring.Poly, beta)
	gal1P = new(MKDecomposedPoly)
	gal1P.poly = make([]*ring.Poly, beta)

	for u := uint64(0); u < beta; u++ {
		gal0Q.poly[u] = galKeys[j-1].key[0].poly[u].CopyNew()
		gal0Q.poly[u].Coeffs = gal0Q.poly[u].Coeffs[:level+1]

		gal0P.poly[u] = galKeys[j-1].key[0].poly[u].CopyNew()
		gal0P.poly[u].Coeffs = gal0P.poly[u].Coeffs[modulus:]

		gal1Q.poly[u] = galKeys[j-1].key[0].poly[u].CopyNew()
		gal1Q.poly[u].Coeffs = gal1Q.poly[u].Coeffs[:level+1]

		gal1P.poly[u] = galKeys[j-1].key[0].poly[u].CopyNew()
		gal1P.poly[u].Coeffs = gal1P.poly[u].Coeffs[modulus:]

	}

	return gal0Q, gal0P, gal1Q, gal1P
}

// SwitchKeysNew perform the key switch for a ciphertext involving a participant
func (eval *mkEvaluator) SwitchKeysNew(ct *MKCiphertext, switchingKey *rlwe.SwitchingKey) (ctOut *MKCiphertext) {
	if ct.ciphertexts.Degree() != 1 {
		panic("Key switch only work for degree 1 ciphertexts")
	}
	// TODO : if it is possible to perform key switch on multi-participant ciphertext
	//implement it by searching for corresponding cipherpart and puting peerID in switchingKey
	var reduce uint64

	level := ct.ciphertexts.Level()

	ringP := GetRingP(eval.params)

	ctOut = NewMKCiphertext(ct.peerIDs, eval.ringQ, eval.params, level)

	cipherQ, cipherP := GInverse(ct.ciphertexts.Value()[1], eval.params, level)

	evakey0Q := new(ring.Poly)
	evakey1Q := new(ring.Poly)
	evakey0P := new(ring.Poly)
	evakey1P := new(ring.Poly)

	pool2Q := eval.ringQ.NewPoly()
	pool3Q := eval.ringQ.NewPoly()
	pool2P := ringP.NewPoly()
	pool3P := ringP.NewPoly()

	reduce = 0

	for i := 0; i < int(eval.params.Beta()); i++ {

		c2QiQ := cipherQ.poly[i]
		c2QiP := cipherP.poly[i]

		//prepare switching key for split computation
		evakey0Q.Coeffs = switchingKey.Value[i][0].Coeffs[:level+1]
		evakey1Q.Coeffs = switchingKey.Value[i][1].Coeffs[:level+1]
		evakey0P.Coeffs = switchingKey.Value[i][0].Coeffs[len(eval.ringQ.Modulus):]
		evakey1P.Coeffs = switchingKey.Value[i][1].Coeffs[len(eval.ringQ.Modulus):]

		if i == 0 {
			eval.ringQ.MulCoeffsMontgomeryConstantLvl(level, evakey0Q, c2QiQ, pool2Q)
			eval.ringQ.MulCoeffsMontgomeryConstantLvl(level, evakey1Q, c2QiQ, pool3Q)
			ringP.MulCoeffsMontgomeryConstant(evakey0P, c2QiP, pool2P)
			ringP.MulCoeffsMontgomeryConstant(evakey1P, c2QiP, pool3P)
		} else {
			eval.ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(level, evakey0Q, c2QiQ, pool2Q)
			eval.ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(level, evakey1Q, c2QiQ, pool3Q)
			ringP.MulCoeffsMontgomeryConstantAndAddNoMod(evakey0P, c2QiP, pool2P)
			ringP.MulCoeffsMontgomeryConstantAndAddNoMod(evakey1P, c2QiP, pool3P)
		}

		//
		if reduce&3 == 3 {
			eval.ringQ.ReduceConstantLvl(level, pool2Q, pool2Q)
			eval.ringQ.ReduceConstantLvl(level, pool3Q, pool3Q)
			ringP.ReduceConstant(pool2P, pool2P)
			ringP.ReduceConstant(pool3P, pool3P)
		}

		reduce++
	}

	eval.ringQ.ReduceLvl(level, pool2Q, pool2Q)
	eval.ringQ.ReduceLvl(level, pool3Q, pool3Q)
	ringP.Reduce(pool2P, pool2P)
	ringP.Reduce(pool3P, pool3P)

	eval.convertor.ModDownSplitNTTPQ(level, pool2Q, pool2P, pool2Q)
	eval.convertor.ModDownSplitNTTPQ(level, pool3Q, pool3P, pool3Q)

	eval.ringQ.AddLvl(level, ct.ciphertexts.Value()[0], pool2Q, ctOut.ciphertexts.Value()[0])
	eval.ringQ.CopyLvl(level, pool3Q, ctOut.ciphertexts.Value()[1])

	return
}

// NewPlaintextFromValue returns a plaintext from the provided values
func (eval *mkEvaluator) NewPlaintextFromValue(values []complex128) *ckks.Plaintext {

	return eval.encoder.EncodeNTTAtLvlNew(eval.params.MaxLevel(), values, eval.params.LogSlots())
}
