package bfv

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
)

// Evaluator is a struct holding the necessary elements to operates the homomorphic operations between ciphertext and/or plaintexts.
// It also holds a small memory pool used to store intermediate computations.
type Evaluator struct {
	bfvcontext    *BfvContext
	basisextender *ring.BasisExtender
	complexscaler *ring.ComplexScaler
	polypool      [4]*ring.Poly
	ctxpool       [3]*Ciphertext
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the ciphertexts and/or plaintexts. It stores a small pool of polynomials
// and ciphertexts that will be used for intermediate values.
func (bfvcontext *BfvContext) NewEvaluator() (evaluator *Evaluator) {

	evaluator = new(Evaluator)
	evaluator.bfvcontext = bfvcontext

	evaluator.basisextender = ring.NewBasisExtender(bfvcontext.contextQ, bfvcontext.contextP)
	evaluator.complexscaler = ring.NewComplexScaler(bfvcontext.t, bfvcontext.contextQ, bfvcontext.contextP)

	for i := 0; i < 4; i++ {
		evaluator.polypool[i] = bfvcontext.contextQP.NewPoly()
	}

	evaluator.ctxpool[0] = bfvcontext.NewCiphertextBig(5)
	evaluator.ctxpool[1] = bfvcontext.NewCiphertextBig(5)
	evaluator.ctxpool[2] = bfvcontext.NewCiphertextBig(5)

	return evaluator
}

func (evaluator *Evaluator) getElemAndCheckBinary(op0, op1, opOut Operand, opOutMinDegree uint64) (el0, el1, elOut *bfvElement, err error) {
	if op0 == nil || op1 == nil || opOut == nil {
		return nil, nil, nil, errors.New("operands cannot be nil")
	}
	if opOut.Degree() < opOutMinDegree {
		return nil, nil, nil, errors.New("receiver operand degree is too small")
	}
	el0, el1, elOut = op0.Element(), op1.Element(), opOut.Element()
	return // TODO: more checks on elements
}

func (evaluator *Evaluator) getElemAndCheckUnary(op0, opOut Operand, opOutMinDegree uint64) (el0, elOut *bfvElement, err error) {
	if op0 == nil || opOut == nil {
		return nil, nil, errors.New("operand cannot be nil")
	}
	if opOut.Degree() < opOutMinDegree {
		return nil, nil, errors.New("receiver operand degree is too small")
	}
	el0, elOut = op0.Element(), opOut.Element()
	return // TODO: more checks on elements
}

// evaluateInPlaceBinary applies the provided function in place on c0 and c1 and returns the result in cOut
func evaluateInPlaceBinary(el0, el1, elOut *bfvElement, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {

	maxDegree := max(el0.Degree(), el1.Degree())
	minDegree := min(el0.Degree(), el1.Degree())

	for i := uint64(0); i < minDegree+1; i++ {
		evaluate(el0.value[i], el1.value[i], elOut.value[i])
	}

	// If the inputs degree differ, copies the remaining degree on the receiver
	var largest *bfvElement
	if el0.Degree() > el1.Degree() {
		largest = el0
	} else if el1.Degree() > el0.Degree() {
		largest = el1
	}
	if largest != nil && largest != elOut { // checks to avoid unnecessary work.
		for i := minDegree + 1; i < maxDegree+1; i++ {
			_ = largest.value[i].Copy(elOut.value[i])
		}
	}
}

// evaluateInPlaceUnary applies the provided function in place on c0 and c1 and returns the result in cOut
func evaluateInPlaceUnary(el0, elOut *bfvElement, evaluate func(*ring.Poly, *ring.Poly)) {
	for i := range el0.value {
		evaluate(el0.value[i], elOut.value[i])
	}
}

// Add adds c0 to c1 and returns the result on cOut.
func (evaluator *Evaluator) Add(op0, op1 Operand, ctOut *Ciphertext) (err error) {
	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max(op0.Degree(), op1.Degree()))
	if err != nil {
		return err
	}
	evaluateInPlaceBinary(el0, el1, elOut, evaluator.bfvcontext.contextQ.Add)
	return
}

// AddNew adds c0 to c1 and creates a new element to store the result.
func (evaluator *Evaluator) AddNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {
	ctOut = evaluator.bfvcontext.NewCiphertext(max(op0.Degree(), op1.Degree()))
	return ctOut, evaluator.Add(op0, op1, ctOut)
}

// AddNoMod adds c0 to c1 without modular reduction, and returns the result on cOut.
func (evaluator *Evaluator) AddNoMod(op0, op1 Operand, ctOut *Ciphertext) (err error) {
	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max(op0.Degree(), op1.Degree()))
	if err != nil {
		return err
	}
	evaluateInPlaceBinary(el0, el1, elOut, evaluator.bfvcontext.contextQ.AddNoMod)
	return nil
}

// AddNoModNew adds c0 to c1 without modular reduction and creates a new element to store the result.
func (evaluator *Evaluator) AddNoModNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {
	ctOut = evaluator.bfvcontext.NewCiphertext(max(op0.Degree(), op1.Degree()))
	return ctOut, evaluator.AddNoMod(op0, op1, ctOut)
}

// Sub subtracts c1 to c0 and returns the result on cOut.
func (evaluator *Evaluator) Sub(op0, op1 Operand, ctOut *Ciphertext) (err error) {
	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max(op0.Degree(), op1.Degree()))
	if err != nil {
		return err
	}
	evaluateInPlaceBinary(el0, el1, elOut, evaluator.bfvcontext.contextQ.Sub)
	return nil
}

// SubNew subtracts c1 to c0 and creates a new element to store the result.
func (evaluator *Evaluator) SubNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {
	ctOut = evaluator.bfvcontext.NewCiphertext(max(op0.Degree(), op1.Degree()))
	return ctOut, evaluator.Sub(op0, op1, ctOut)
}

// SubNoMod subtracts c1 to c0 without modular reduction and returns the result on cOut.
func (evaluator *Evaluator) SubNoMod(op0, op1 Operand, ctOut *Ciphertext) (err error) {
	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max(op0.Degree(), op1.Degree()))
	if err != nil {
		return err
	}
	evaluateInPlaceBinary(el0, el1, elOut, evaluator.bfvcontext.contextQ.SubNoMod)
	return nil
}

// SubNoModNew subtracts c1 to c0 without modular reduction and creates a new element to store the result.
func (evaluator *Evaluator) SubNoModNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {
	ctOut = evaluator.bfvcontext.NewCiphertext(max(op0.Degree(), op1.Degree()))
	return ctOut, evaluator.SubNoMod(op0, op1, ctOut)
}

// Neg negates c0 and returns the result on cOut.
func (evaluator *Evaluator) Neg(op Operand, ctOut *Ciphertext) error {
	el0, elOut, err := evaluator.getElemAndCheckUnary(op, ctOut, op.Degree())
	if err != nil {
		return err
	}
	evaluateInPlaceUnary(el0, elOut, evaluator.bfvcontext.contextQ.Neg)
	return nil
}

// Neg negates c0 and creates a new element to store the result.
func (evaluator *Evaluator) NegNew(op Operand) (ctOut *Ciphertext, err error) {
	ctOut = evaluator.bfvcontext.NewCiphertext(op.Degree())
	return ctOut, evaluator.Neg(op, ctOut)
}

// Reduce applies a modular reduction on c0 and returns the result on cOut.
func (evaluator *Evaluator) Reduce(op Operand, ctOut *Ciphertext) error {
	el0, elOut, err := evaluator.getElemAndCheckUnary(op, ctOut, op.Degree())
	if err != nil {
		return err
	}
	evaluateInPlaceUnary(el0, elOut, evaluator.bfvcontext.contextQ.Reduce)
	return nil
}

// Reduce applies a modular reduction on c0 and creates a new element to store the result.
func (evaluator *Evaluator) ReduceNew(op Operand) (ctOut *Ciphertext, err error) {
	ctOut = evaluator.bfvcontext.NewCiphertext(op.Degree())
	return ctOut, evaluator.Reduce(op, ctOut)
}

// MulScalar multiplies c0 by an uint64 scalar and returns the result on cOut.
func (evaluator *Evaluator) MulScalar(op Operand, scalar uint64, ctOut *Ciphertext) error {

	el0, elOut, err := evaluator.getElemAndCheckUnary(op, ctOut, op.Degree())
	if err != nil {
		return err
	}
	fun := func(el, elOut *ring.Poly) { evaluator.bfvcontext.contextQ.MulScalar(el, scalar, elOut) }
	evaluateInPlaceUnary(el0, elOut, fun)
	return nil
}

// MulScalarNew multiplies c0 by an uint64 scalar and creates a new element to store the result.
func (evaluator *Evaluator) MulScalarNew(op Operand, scalar uint64) (ctOut *Ciphertext, err error) {
	ctOut = evaluator.bfvcontext.NewCiphertext(op.Degree())
	return ctOut, evaluator.MulScalar(op, scalar, ctOut)
}

// tensorAndRescales computes (ct0 x ct1) * (t/Q) and stores the result on cOut.
func (evaluator *Evaluator) tensorAndRescale(ct0, ct1, cOut *bfvElement) {

	// Prepares the ciphertexts for the Tensoring by extending their
	// basis from Q to QP and transforming them in NTT form
	c0 := evaluator.ctxpool[0]
	c1 := evaluator.ctxpool[1]
	tmpCout := evaluator.ctxpool[2]

	if ct0 == ct1 {

		for i := range ct0.value {
			evaluator.basisextender.ExtendBasis(ct0.value[i], c0.value[i])
			evaluator.bfvcontext.contextQP.NTT(c0.value[i], c0.value[i])
		}

	} else {

		for i := range ct0.value {
			evaluator.basisextender.ExtendBasis(ct0.value[i], c0.value[i])
			evaluator.bfvcontext.contextQP.NTT(c0.value[i], c0.value[i])
		}

		for i := range ct1.value {
			evaluator.basisextender.ExtendBasis(ct1.value[i], c1.value[i])
			evaluator.bfvcontext.contextQP.NTT(c1.value[i], c1.value[i])
		}
	}

	// Tensoring : multiplies each elements of the ciphertexts together
	// and adds them to their correspongint position in the new ciphertext
	// based on their respective degree

	// Case where both BfvElements are of degree 1
	if ct0.Degree() == 1 && ct1.Degree() == 1 {

		c_00 := evaluator.polypool[0]
		c_01 := evaluator.polypool[1]

		d0 := tmpCout.value[0]
		d1 := tmpCout.value[1]
		d2 := tmpCout.value[2]

		evaluator.bfvcontext.contextQP.MForm(c0.value[0], c_00)
		evaluator.bfvcontext.contextQP.MForm(c0.value[1], c_01)

		// Squaring case
		if ct0 == ct1 {

			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00, c0.value[0], d0) // c0 = c0[0]*c0[0]
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00, c0.value[1], d1) // c1 = 2*c0[0]*0[1]
			evaluator.bfvcontext.contextQP.Add(d1, d1, d1)
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_01, c0.value[1], d2) // c2 = c0[1]*c0[1]

			// Normal case
		} else {

			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00, c1.value[0], d0) // c0 = c0[0]*c0[0]
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00, c1.value[1], d1)
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomeryAndAddNoMod(c_01, c1.value[0], d1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_01, c1.value[1], d2)            // c2 = c0[1]*c1[1]
		}

		// Case where both BfvElements are not of degree 1
	} else {

		for i := 0; i < len(ct0.Value())+len(ct1.Value()); i++ {
			tmpCout.value[i].Zero()
		}

		// Squaring case
		if ct0 == ct1 {

			c_00 := evaluator.ctxpool[1]

			for i := range ct0.value {
				evaluator.bfvcontext.contextQP.MForm(c0.value[i], c_00.value[i])
			}

			for i := uint64(0); i < ct0.Degree()+1; i++ {
				for j := i + 1; j < ct0.Degree()+1; j++ {
					evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00.value[i], c0.value[j], tmpCout.value[i+j])
					evaluator.bfvcontext.contextQP.Add(tmpCout.value[i+j], tmpCout.value[i+j], tmpCout.value[i+j])
				}
			}

			for i := uint64(0); i < ct0.Degree()+1; i++ {
				evaluator.bfvcontext.contextQP.MulCoeffsMontgomeryAndAdd(c_00.value[i], c0.value[i], tmpCout.value[i<<1])
			}

			// Normal case
		} else {
			for i := range ct0.value {
				evaluator.bfvcontext.contextQP.MForm(c0.value[i], c0.value[i])
				for j := range ct1.value {
					evaluator.bfvcontext.contextQP.MulCoeffsMontgomeryAndAdd(c0.value[i], c1.value[j], tmpCout.value[i+j])
				}
			}
		}
	}

	// Applies the inverse NTT to the ciphertext, scales the down ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range cOut.value {
		evaluator.bfvcontext.contextQP.InvNTT(tmpCout.value[i], tmpCout.value[i])
		evaluator.complexscaler.Scale(tmpCout.value[i], cOut.value[i])
	}
}

// Mul multiplies c0 by c1 and returns the result on cOut.
func (evaluator *Evaluator) Mul(op0, op1 Operand, ctOut *Ciphertext) (err error) {

	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, op0.Degree()+op1.Degree())
	if err != nil {
		return err
	}
	evaluator.tensorAndRescale(el0, el1, elOut)
	return nil
}

// MulNew multiplies c0 by c1 and creates a new element to store the result.
func (evaluator *Evaluator) MulNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.bfvcontext.NewCiphertext(op0.Degree() + op1.Degree())
	return ctOut, evaluator.Mul(op0, op1, ctOut)
}

// switchKeys compute cOut = [c0 + c2*evakey[0], c1 + c2*evakey[1]].
func (evaluator *Evaluator) switchKeys(c0, c1, c2 *ring.Poly, evakey *SwitchingKey, elOut *bfvElement) {

	var mask, reduce, bitLog uint64

	if c0 != elOut.value[0] {
		c0.Copy(elOut.value[0])
	}
	if c1 != elOut.value[1] {
		c1.Copy(elOut.value[1])
	}

	c2_qi_w := evaluator.polypool[3]

	mask = uint64((1 << evakey.bitDecomp) - 1)

	reduce = 0

	for i := range evaluator.bfvcontext.contextQ.Modulus {

		bitLog = uint64(len(evakey.evakey[i]))

		for j := uint64(0); j < bitLog; j++ {
			//c2_qi_w = (c2_qi_w >> (w*z)) & (w-1)
			for u := uint64(0); u < evaluator.bfvcontext.n; u++ {
				for v := range evaluator.bfvcontext.contextQ.Modulus {
					c2_qi_w.Coeffs[v][u] = (c2.Coeffs[i][u] >> (j * evakey.bitDecomp)) & mask
				}
			}

			evaluator.bfvcontext.contextQ.NTT(c2_qi_w, c2_qi_w)

			evaluator.bfvcontext.contextQ.MulCoeffsMontgomeryAndAddNoMod(evakey.evakey[i][j][0], c2_qi_w, elOut.value[0])
			evaluator.bfvcontext.contextQ.MulCoeffsMontgomeryAndAddNoMod(evakey.evakey[i][j][1], c2_qi_w, elOut.value[1])

			if reduce&7 == 7 {
				evaluator.bfvcontext.contextQ.Reduce(elOut.value[0], elOut.value[0])
				evaluator.bfvcontext.contextQ.Reduce(elOut.value[1], elOut.value[1])
			}

			reduce += 1
		}
	}

	if (reduce-1)&7 != 7 {
		evaluator.bfvcontext.contextQ.Reduce(elOut.value[0], elOut.value[0])
		evaluator.bfvcontext.contextQ.Reduce(elOut.value[1], elOut.value[1])
	}
}

// relinearize is a methode common to Relinearize and RelinearizeNew. It switches cIn out in the NTT domain, applies the keyswitch, and returns the result out of the NTT domain.
func (evaluator *Evaluator) relinearize(el *bfvElement, evakey *EvaluationKey, elOut *bfvElement) {

	evaluator.bfvcontext.contextQ.NTT(el.value[0], elOut.value[0])
	evaluator.bfvcontext.contextQ.NTT(el.value[1], elOut.value[1])

	for deg := uint64(el.Degree()); deg > 1; deg-- {
		evaluator.switchKeys(elOut.value[0], elOut.value[1], el.value[deg], evakey.evakey[deg-2], elOut)
	}

	if len(elOut.value) > 2 {
		elOut.SetValue(elOut.value[:2])
	}

	evaluator.bfvcontext.contextQ.InvNTT(elOut.value[0], elOut.value[0])
	evaluator.bfvcontext.contextQ.InvNTT(elOut.value[1], elOut.value[1])
}

// Relinearize relinearize the ciphertext cIn of degree > 1 until it is of degree 1 and returns the result on cOut.
//
// Requires a correct evaluation key as additional input :
//
// - it must match the secret-key that was used to create the public key under which the current cIn is encrypted.
//
// - it must be of degree high enough to relinearize the input ciphertext to degree 1 (ex. a ciphertext
//of degree 3 will require that the evaluation key stores the keys for both degree 3 and 2 ciphertexts).
func (evaluator *Evaluator) Relinearize(op Operand, evakey *EvaluationKey, ctOut *Ciphertext) error {

	el, elOut, err := evaluator.getElemAndCheckUnary(op, ctOut, 1)

	if err != nil {
		return err
	}

	if int(el.Degree()-1) > len(evakey.evakey) {
		return errors.New("error : ciphertext degree too large to allow relinearization")
	}

	if el.Degree() < 2 {
		return nil
	}

	evaluator.relinearize(el, evakey, elOut)
	return nil
}

// Relinearize relinearize the ciphertext cIn of degree > 1 until it is of degree 1 and creates a new ciphertext to store the result.
//
// Requires a correct evaluation key as additional input :
//
// - it must match the secret-key that was used to create the public key under which the current cIn is encrypted
//
// - it must be of degree high enough to relinearize the input ciphertext to degree 1 (ex. a ciphertext
// of degree 3 will require that the evaluation key stores the keys for both degree 3 and 2 ciphertexts).
func (evaluator *Evaluator) RelinearizeNew(op Operand, evakey *EvaluationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.bfvcontext.NewCiphertext(1)

	return ctOut, evaluator.Relinearize(op, evakey, ctOut)
}

// SwitchKeys applies the key-switching procedure to the ciphertext cIn and returns the result on cOut. It requires as an additional input a valide switching-key :
// it must encrypt the target key under the public key under which cIn is currently encrypted.
func (evaluator *Evaluator) SwitchKeys(op Operand, switchkey *SwitchingKey, ctOut *Ciphertext) (err error) {

	if op.Degree() != 1 {
		return errors.New("error : ciphertext must be of degree 1 to allow key switching")
	}

	el, elOut, err := evaluator.getElemAndCheckUnary(op, ctOut, 1)
	if err != nil {
		return err
	}

	var c2 *ring.Poly
	if el == elOut {
		c2 = evaluator.polypool[1]
		el.value[1].Copy(c2)
	} else {
		c2 = el.value[1]
	}
	_ = el.NTT(evaluator.bfvcontext, elOut)
	evaluator.switchKeys(elOut.value[0], elOut.value[1], c2, switchkey, elOut)
	_ = elOut.InvNTT(evaluator.bfvcontext, elOut)

	return nil
}

// SwitchKeys applies the key-switching procedure to the ciphertext cIn and creates a new ciphertext to store the result. It requires as an additional input a valide switching-key :
// it must encrypt the target key under the public key under which cIn is currently encrypted.
func (evaluator *Evaluator) SwitchKeysNew(op Operand, switchkey *SwitchingKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.bfvcontext.NewCiphertext(1)
	return ctOut, evaluator.SwitchKeys(op, switchkey, ctOut)
}

// RotateColumns rotates the columns of c0 by k position to the left and returns the result on c1. As an additional input it requires a rotationkeys :
//
// - it must either store all the left and right power of 2 rotations or the specific rotation that is asked.
//
// If only the power of two rotations are stored, the numbers k and n/2-k will be decomposed in base 2 and the rotation with the least
// hamming weight will be chosen, then the specific rotation will be computed as a sum of powers of two rotations.
func (evaluator *Evaluator) RotateColumns(op Operand, k uint64, evakey *RotationKeys, ctOut *Ciphertext) (err error) {

	el, elOut, err := evaluator.getElemAndCheckUnary(op, ctOut, 1)
	if err != nil {
		return err
	}

	k &= ((evaluator.bfvcontext.n >> 1) - 1)

	if k == 0 {
		if el != elOut {
			if err = el.Copy(elOut); err != nil {
				return err
			}
		}

		return nil
	}

	if elOut.Degree() != el.Degree() {
		return errors.New("cannot rotate -> receiver degree doesn't match input degree ")
	}

	if el.Degree() > 1 {
		return errors.New("cannot rotate -> input and or output degree not 0 or 1")
	}

	context := evaluator.bfvcontext.contextQ

	if el.Degree() == 0 {

		if elOut != el {

			if el.IsNTT() {
				ring.PermuteNTT(el.value[0], evaluator.bfvcontext.galElRotColLeft[k], elOut.value[0])
			} else {
				context.Permute(el.value[0], evaluator.bfvcontext.galElRotColLeft[k], elOut.value[0])
			}

		} else {

			if el.IsNTT() {
				ring.PermuteNTT(el.value[0], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[0])
			} else {
				context.Permute(el.value[0], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[0])
			}

			context.Copy(evaluator.polypool[0], elOut.value[0])
		}

		return nil

	} else {
		// Looks in the rotationkey if the corresponding rotation has been generated
		if evakey.evakey_rot_col_L[k] != nil {

			if el.IsNTT() {

				ring.PermuteNTT(el.value[0], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[0])
				ring.PermuteNTT(el.value[1], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[1])

				context.Copy(evaluator.polypool[0], elOut.value[0])
				context.Copy(evaluator.polypool[1], elOut.value[1])

				context.InvNTT(evaluator.polypool[1], evaluator.polypool[1])

				evaluator.switchKeys(elOut.value[0], elOut.value[1], evaluator.polypool[1], evakey.evakey_rot_col_L[k], elOut)

			} else {

				context.Permute(el.value[0], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[0])
				context.Permute(el.value[1], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[1])

				context.NTT(evaluator.polypool[0], elOut.value[0])
				context.NTT(evaluator.polypool[1], elOut.value[1])

				evaluator.switchKeys(elOut.value[0], elOut.value[1], evaluator.polypool[1], evakey.evakey_rot_col_L[k], elOut)

				context.InvNTT(elOut.value[0], elOut.value[0])
				context.InvNTT(elOut.value[1], elOut.value[1])
			}

			return nil

		} else {

			// If not looks if the left and right pow2 rotations have been generated
			has_pow2_rotations := true
			for i := uint64(1); i < evaluator.bfvcontext.n>>1; i <<= 1 {
				if evakey.evakey_rot_col_L[i] == nil || evakey.evakey_rot_col_R[i] == nil {
					has_pow2_rotations = false
					break
				}
			}

			// If yes, computes the least amount of rotation between k to the left and n/2 -k to the right required to apply the demanded rotation
			if has_pow2_rotations {

				if hammingWeight64(k) <= hammingWeight64((evaluator.bfvcontext.n>>1)-k) {
					evaluator.rotateColumnsLPow2(el, k, evakey, elOut)
				} else {
					evaluator.rotateColumnsRPow2(el, (evaluator.bfvcontext.n>>1)-k, evakey, elOut)
				}

				return nil

				// Else returns an error indicating that the keys have not been generated
			} else {
				return errors.New("error : specific rotation and pow2 rotations have not been generated")
			}
		}
	}
}

// rotateColumnsLPow2 applies the Galois Automorphism on the element, rotating the element by k positions to the left.
func (evaluator *Evaluator) rotateColumnsLPow2(c0 *bfvElement, k uint64, evakey *RotationKeys, c1 *bfvElement) {
	evaluator.rotateColumnsPow2(c0, evaluator.bfvcontext.gen, k, evakey.evakey_rot_col_L, c1)
}

// rotateColumnsRPow2 applies the Galois Endomorphism on the element, rotating the element by k positions to the right.
func (evaluator *Evaluator) rotateColumnsRPow2(c0 *bfvElement, k uint64, evakey *RotationKeys, c1 *bfvElement) {
	evaluator.rotateColumnsPow2(c0, evaluator.bfvcontext.genInv, k, evakey.evakey_rot_col_R, c1)
}

// rotateColumnsPow2 rotates c0 by k position (left or right depending on the input), decomposing k as a sum of power of 2 rotations, and returns the result on c1.
func (evaluator *Evaluator) rotateColumnsPow2(c0 *bfvElement, generator, k uint64, evakey_rot_col map[uint64]*SwitchingKey, c1 *bfvElement) {

	var mask, evakey_index uint64

	context := evaluator.bfvcontext.contextQ

	mask = (evaluator.bfvcontext.n << 1) - 1

	evakey_index = 1

	if c0.IsNTT() {
		c0.Copy(c1)
	} else {
		for i := range c0.value {
			context.NTT(c0.value[i], c1.value[i])
		}
	}

	// Applies the galois automorphism and the switching-key process
	for k > 0 {

		if k&1 == 1 {

			if c0.Degree() == 0 {

				ring.PermuteNTT(c1.value[0], generator, evaluator.polypool[0])
				context.Copy(evaluator.polypool[0], c1.value[0])

			} else {

				ring.PermuteNTT(c1.value[0], generator, evaluator.polypool[0])
				ring.PermuteNTT(c1.value[1], generator, evaluator.polypool[1])

				context.Copy(evaluator.polypool[0], c1.value[0])
				context.Copy(evaluator.polypool[1], c1.value[1])
				context.InvNTT(evaluator.polypool[1], evaluator.polypool[2])

				evaluator.switchKeys(evaluator.polypool[0], evaluator.polypool[1], evaluator.polypool[2], evakey_rot_col[evakey_index], c1)
			}

		}

		generator *= generator
		generator &= mask

		evakey_index <<= 1
		k >>= 1
	}

	if !c0.IsNTT() {
		for i := range c1.value {
			context.InvNTT(c1.value[i], c1.value[i])
		}
	}
}

// RotateRows swaps the rows of c0 and returns the result on c1.
func (evaluator *Evaluator) RotateRows(op Operand, evakey *RotationKeys, ctOut *Ciphertext) error {

	c0, c1 := op.Element(), ctOut.Element()

	if c1.Degree() != c0.Degree() {
		return errors.New("cannot rotate -> receiver degree doesn't match input degree ")
	}

	if c0.Degree() > 1 {
		return errors.New("cannot rotate -> input and or output degree not 0 or 1")
	}

	if evakey.evakey_rot_row == nil {
		return errors.New("error : rows rotation key not generated")
	}

	context := evaluator.bfvcontext.contextQ

	if c0.Degree() == 0 {

		if c0.IsNTT() {

			if c0 != c1 {

				ring.PermuteNTT(c0.value[0], evaluator.bfvcontext.galElRotRow, c1.value[0])

			} else {

				ring.PermuteNTT(c0.value[0], evaluator.bfvcontext.galElRotRow, evaluator.polypool[0])
				context.Copy(evaluator.polypool[0], c1.value[0])
			}

		} else {

			if c0 != c1 {

				context.Permute(c0.value[0], evaluator.bfvcontext.galElRotRow, c1.value[0])

			} else {

				context.Permute(c0.value[0], evaluator.bfvcontext.galElRotRow, evaluator.polypool[0])
				context.Copy(evaluator.polypool[0], c1.value[0])
			}
		}

	} else {

		if c0.IsNTT() {

			if c0 != c1 {

				ring.PermuteNTT(c0.value[0], evaluator.bfvcontext.galElRotRow, c1.value[0])
				ring.PermuteNTT(c0.value[1], evaluator.bfvcontext.galElRotRow, c1.value[1])

				context.InvNTT(c1.value[1], evaluator.polypool[1])

				evaluator.switchKeys(c1.value[0], c1.value[1], evaluator.polypool[1], evakey.evakey_rot_row, c1)

			} else {

				ring.PermuteNTT(c0.value[0], evaluator.bfvcontext.galElRotRow, evaluator.polypool[0])
				ring.PermuteNTT(c0.value[1], evaluator.bfvcontext.galElRotRow, evaluator.polypool[1])

				context.Copy(evaluator.polypool[0], c1.value[0])
				context.Copy(evaluator.polypool[1], c1.value[1])

				context.InvNTT(evaluator.polypool[1], evaluator.polypool[1])

				evaluator.switchKeys(c1.value[0], c1.value[1], evaluator.polypool[1], evakey.evakey_rot_row, c1)
			}

		} else {

			if c0 != c1 {

				context.Permute(c0.value[0], evaluator.bfvcontext.galElRotRow, c1.value[0])
				context.Permute(c0.value[1], evaluator.bfvcontext.galElRotRow, c1.value[1])

				context.Copy(c1.value[1], evaluator.polypool[1])

				context.NTT(c1.value[0], c1.value[0])
				context.NTT(c1.value[1], c1.value[1])

				evaluator.switchKeys(c1.value[0], c1.value[1], evaluator.polypool[1], evakey.evakey_rot_row, c1)

				context.InvNTT(c1.value[0], c1.value[0])
				context.InvNTT(c1.value[1], c1.value[1])

			} else {

				context.Permute(c0.value[0], evaluator.bfvcontext.galElRotRow, evaluator.polypool[0])
				context.Permute(c0.value[1], evaluator.bfvcontext.galElRotRow, evaluator.polypool[1])

				context.NTT(evaluator.polypool[0], c1.value[0])
				context.NTT(evaluator.polypool[1], c1.value[1])

				evaluator.switchKeys(c1.value[0], c1.value[1], evaluator.polypool[1], evakey.evakey_rot_row, c1)

				context.InvNTT(c1.value[0], c1.value[0])
				context.InvNTT(c1.value[1], c1.value[1])
			}
		}
	}

	return nil
}

// InnerSum computs the inner sum of c0 and returns the result on c1. It requires a rotation key storing all the left power of two rotations.
// The resulting vector will be of the form [sum, sum, .., sum, sum ].
func (evaluator *Evaluator) InnerSum(op Operand, evakey *RotationKeys, ctOut *Ciphertext) error {

	el, elOut, err := evaluator.getElemAndCheckUnary(op, ctOut, 1)
	if err != nil {
		return err
	}

	if el.Degree() != 1 {
		return errors.New("error : ciphertext must be of degree 1 to allow Galois Auotomorphism (required for inner sum)")
	}

	if elOut.Degree() != 1 {
		return errors.New("error : receiver ciphertext must be of degree 1 to allow Galois Automorphism (required for inner sum)")
	}

	cTmp := evaluator.bfvcontext.NewCiphertext(1)

	if el != elOut {
		if err := elOut.Copy(el); err != nil {
			return err
		}
	}

	for i := uint64(1); i < evaluator.bfvcontext.n>>1; i <<= 1 {
		if err := evaluator.RotateColumns(elOut, i, evakey, cTmp); err != nil {
			return err
		}
		evaluator.Add(cTmp.bfvElement, elOut, elOut.Ciphertext())
	}

	if err := evaluator.RotateRows(elOut, evakey, cTmp); err != nil {
		return err
	}
	evaluator.Add(elOut, cTmp.bfvElement, elOut.Ciphertext())

	return nil

}
