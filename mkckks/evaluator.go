package mkckks

import (
	"sort"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
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
	RelinInPlace(ct *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey)
	Rescale(c *MKCiphertext, out *MKCiphertext)
	Rotate(c *MKCiphertext, n int, keys []*MKEvalGalKey) *MKCiphertext
	SwitchKeysNew(ct *MKCiphertext, switchingKey *MKSwitchingKey) (ctOut *MKCiphertext)
	NewPlaintextFromValue([]complex128) *ckks.Plaintext
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

	ringQ := GetRingQ(params)
	ringP := GetRingP(params)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	sampler := GetGaussianSampler(params, ringQ, prng)
	convertor := ring.NewFastBasisExtender(ringQ, ringP)

	return &mkEvaluator{
		ckksEval:        ckks.NewEvaluator(params, ckks.EvaluationKey{}),
		params:          params,
		ringQ:           ringQ,
		ringP:           ringP,
		samplerGaussian: sampler,
		convertor:       convertor,
		encoder:         ckks.NewEncoder(params),
	}
}

// Add adds the ciphertexts component wise and expend their list of involved peers. A new ciphertext is returned
func (eval *mkEvaluator) Add(c0 *MKCiphertext, c1 *MKCiphertext) *MKCiphertext {

	if c0 == nil || c1 == nil || c0.ciphertexts == nil || c1.ciphertexts == nil {
		panic("Uninitialized ciphertexts")
	}

	padded1, padded2 := PadCiphers(c0, c1, eval.params)

	out := NewMKCiphertext(padded1.peerIDs, eval.ringQ, eval.params, padded1.ciphertexts.Level())

	out.ciphertexts = eval.ckksEval.AddNew(padded1.ciphertexts, padded2.ciphertexts)
	return out
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

// Mul will compute the tensor product and output a ciphertext with degree (k+1)**2
func (eval *mkEvaluator) Mul(c1 *MKCiphertext, c2 *MKCiphertext) *MKCiphertext {

	padded1, padded2 := PadCiphers(c1, c2, eval.params)

	nbrElements := padded1.ciphertexts.Degree() + 1 // k+1

	outputDegree := nbrElements * nbrElements // (k+1)**2

	el1 := padded1.ciphertexts
	el2 := padded2.ciphertexts
	level := utils.MinUint64(el1.Level(), el2.Level())

	out := new(MKCiphertext)
	out.ciphertexts = ckks.NewCiphertext(eval.params, outputDegree-1, level, el1.Scale()*el2.Scale())
	out.peerIDs = padded1.peerIDs

	if !el1.IsNTT() {
		panic("cannot MulRelin: op0 must be in NTT")
	}

	if !el2.IsNTT() {
		panic("cannot MulRelin: op1 must be in NTT")
	}

	ringQ := eval.ringQ

	/*tmp1 := ringQ.NewPoly()
	tmp2 := ringQ.NewPoly()*/

	for i, v1 := range el1.Value() {

		//ringQ.MFormLvl(level, v1, tmp1)

		for j, v2 := range el2.Value() {
			//ringQ.MFormLvl(level, v2, tmp2)
			ringQ.MulCoeffs(v1, v2, out.ciphertexts.Value()[int(nbrElements)*i+j])
		}
	}

	return out
}

// Relinearize a ciphertext after a multiplication
func (eval *mkEvaluator) RelinInPlace(ct *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey) {

	sort.Slice(evalKeys, func(i, j int) bool { return evalKeys[i].peerID < evalKeys[j].peerID })
	sort.Slice(publicKeys, func(i, j int) bool { return publicKeys[i].peerID < publicKeys[j].peerID })

	checkParticipantsEvalKey(ct.peerIDs, evalKeys)
	checkParticipantsPubKey(ct.peerIDs, publicKeys)

	Relinearization(evalKeys, publicKeys, ct, eval.params)
}

// Rescale takes a ciphertext at level l reduces it until it reaches its original
// this function is the same as in ckks/evaluator.go
func (eval *mkEvaluator) Rescale(c *MKCiphertext, out *MKCiphertext) {

	eval.ckksEval.Rescale(c.ciphertexts, eval.params.Scale(), c.ciphertexts)
}

// Rotate rotate the columns of the ciphertext by n to the left and return the result in a new ciphertext
func (eval *mkEvaluator) Rotate(c *MKCiphertext, n int, keys []*MKEvalGalKey) *MKCiphertext {

	sort.Slice(keys, func(i, j int) bool { return keys[i].peerID < keys[j].peerID })

	checkParticipantsGalKey(c.peerIDs, keys)

	out := NewMKCiphertext(c.peerIDs, eval.ringQ, eval.params, c.ciphertexts.Level())

	galEl := eval.params.GaloisElementForColumnRotationBy(n)

	level := c.ciphertexts.Level()

	ringQP := GetRingQP(eval.params)

	k := uint64(len(c.peerIDs))

	res := make([]*ring.Poly, k+1)

	restmpQ := make([]*ring.Poly, k+1)
	restmpP := make([]*ring.Poly, k+1)

	for i := uint64(0); i < k+1; i++ {
		restmpQ[i] = eval.ringQ.NewPoly()
		restmpP[i] = eval.ringP.NewPoly()
		res[i] = eval.ringQ.NewPoly()
	}

	for i := uint64(1); i <= k; i++ {

		gal0Q, gal0P, gal1Q, gal1P := prepareGaloisEvaluationKey(i, level, uint64(len(eval.ringQ.Modulus)), eval.params.Beta(), keys)

		permutedCipher := eval.ringQ.NewPoly() // apply rotation to the ciphertext
		index := ring.PermuteNTTIndex(galEl, ringQP.N)
		ring.PermuteNTTWithIndexLvl(level, c.ciphertexts.Value()[i], index, permutedCipher)

		decomposedPermutedQ, decomposedPermutedP := GInverse(permutedCipher, eval.params, level)

		res0P := Dot(decomposedPermutedP, gal0P, eval.ringP) // dot product and add in c0''
		res0Q := DotLvl(level, decomposedPermutedQ, gal0Q, eval.ringQ)

		eval.ringP.Add(restmpP[0], res0P, restmpP[0])
		eval.ringQ.AddLvl(level, restmpQ[0], res0Q, restmpQ[0])

		restmpP[i] = Dot(decomposedPermutedP, gal1P, eval.ringP) // dot product and put in ci''
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
		eval.ringQ.CopyLvl(level, tmpModDown, res[i])
	}

	out.ciphertexts.Ciphertext().SetValue(res)

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

// SwitchKeysNew perform the key switch for a ciphertext involving one participant.
// After the key switch, the participants's secret Key must be update with SetSecretKey
func (eval *mkEvaluator) SwitchKeysNew(ct *MKCiphertext, switchingKey *MKSwitchingKey) (ctOut *MKCiphertext) {

	if ct.ciphertexts.Degree() != 1 {
		panic("Key switch only work for degree 1 ciphertexts")
	}

	index := getCiphertextIndex(switchingKey.peerID, ct)

	if index == 0 {
		panic("Participant not involved in ciphertext. Key switch impossible")
	}

	var reduce uint64

	level := ct.ciphertexts.Level()

	ctOut = NewMKCiphertext(ct.peerIDs, eval.ringQ, eval.params, level)

	cipherQ, cipherP := GInverse(ct.ciphertexts.Value()[index], eval.params, level)

	evakey0Q := new(ring.Poly)
	evakey1Q := new(ring.Poly)
	evakey0P := new(ring.Poly)
	evakey1P := new(ring.Poly)

	pool2Q := eval.ringQ.NewPoly()
	pool3Q := eval.ringQ.NewPoly()
	pool2P := eval.ringP.NewPoly()
	pool3P := eval.ringP.NewPoly()

	reduce = 0

	for i := 0; i < int(eval.params.Beta()); i++ {

		c2QiQ := cipherQ.poly[i]
		c2QiP := cipherP.poly[i]

		//prepare switching key for split computation
		evakey0Q.Coeffs = switchingKey.key[0].poly[i].Coeffs[:level+1]
		evakey1Q.Coeffs = switchingKey.key[1].poly[i].Coeffs[:level+1]
		evakey0P.Coeffs = switchingKey.key[0].poly[i].Coeffs[len(eval.ringQ.Modulus):]
		evakey1P.Coeffs = switchingKey.key[1].poly[i].Coeffs[len(eval.ringQ.Modulus):]

		if i == 0 {
			eval.ringQ.MulCoeffsMontgomeryConstantLvl(level, evakey0Q, c2QiQ, pool2Q)
			eval.ringQ.MulCoeffsMontgomeryConstantLvl(level, evakey1Q, c2QiQ, pool3Q)
			eval.ringP.MulCoeffsMontgomeryConstant(evakey0P, c2QiP, pool2P)
			eval.ringP.MulCoeffsMontgomeryConstant(evakey1P, c2QiP, pool3P)
		} else {
			eval.ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(level, evakey0Q, c2QiQ, pool2Q)
			eval.ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(level, evakey1Q, c2QiQ, pool3Q)
			eval.ringP.MulCoeffsMontgomeryConstantAndAddNoMod(evakey0P, c2QiP, pool2P)
			eval.ringP.MulCoeffsMontgomeryConstantAndAddNoMod(evakey1P, c2QiP, pool3P)
		}

		//
		if reduce&3 == 3 {
			eval.ringQ.ReduceConstantLvl(level, pool2Q, pool2Q)
			eval.ringQ.ReduceConstantLvl(level, pool3Q, pool3Q)
			eval.ringP.ReduceConstant(pool2P, pool2P)
			eval.ringP.ReduceConstant(pool3P, pool3P)
		}

		reduce++
	}

	eval.ringQ.ReduceLvl(level, pool2Q, pool2Q)
	eval.ringQ.ReduceLvl(level, pool3Q, pool3Q)
	eval.ringP.Reduce(pool2P, pool2P)
	eval.ringP.Reduce(pool3P, pool3P)

	eval.convertor.ModDownSplitNTTPQ(level, pool2Q, pool2P, pool2Q)
	eval.convertor.ModDownSplitNTTPQ(level, pool3Q, pool3P, pool3Q)

	eval.ringQ.AddLvl(level, ct.ciphertexts.Value()[0], pool2Q, ctOut.ciphertexts.Value()[0])
	eval.ringQ.CopyLvl(level, pool3Q, ctOut.ciphertexts.Value()[index])

	return
}

// NewPlaintextFromValue returns a plaintext from the provided values
func (eval *mkEvaluator) NewPlaintextFromValue(values []complex128) *ckks.Plaintext {

	return eval.encoder.EncodeNTTAtLvlNew(eval.params.MaxLevel(), values, eval.params.LogSlots())
}

func getCiphertextIndex(peerID uint64, ct *MKCiphertext) uint64 {

	for i, id := range ct.peerIDs {

		if id == peerID {
			return uint64(i + 1)
		}
	}

	return 0
}

func checkParticipantsEvalKey(peerID []uint64, evalKeys []*MKEvaluationKey) {

	for i, id := range peerID {
		if id != evalKeys[i].peerID {
			panic("Incorrect evaluation keys for the given ciphertexts")
		}
	}
}

func checkParticipantsPubKey(peerID []uint64, pubKeys []*MKPublicKey) {

	for i, id := range peerID {
		if id != pubKeys[i].peerID {
			panic("Incorrect public keys for the given ciphertexts")
		}
	}
}

func checkParticipantsGalKey(peerID []uint64, galKeys []*MKEvalGalKey) {

	for i, id := range peerID {
		if id != galKeys[i].peerID {
			panic("Incorrect galois evaluation keys for the given ciphertexts")
		}
	}
}
