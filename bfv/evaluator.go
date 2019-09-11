package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
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
func (bfvcontext *BfvContext) NewEvaluator() (evaluator *Evaluator, err error) {

	evaluator = new(Evaluator)
	evaluator.bfvcontext = bfvcontext

	if evaluator.basisextender, err = ring.NewBasisExtender(bfvcontext.contextQ, bfvcontext.contextP); err != nil {
		return nil, err
	}

	if evaluator.complexscaler, err = ring.NewComplexScaler(bfvcontext.t, bfvcontext.contextQ, bfvcontext.contextP); err != nil {
		return nil, err
	}

	for i := 0; i < 4; i++ {
		evaluator.polypool[i] = bfvcontext.contextQP.NewPoly()
	}

	evaluator.ctxpool[0] = bfvcontext.NewCiphertextBig(5)
	evaluator.ctxpool[1] = bfvcontext.NewCiphertextBig(5)
	evaluator.ctxpool[2] = bfvcontext.NewCiphertextBig(5)

	return evaluator, nil
}

// Add adds c0 to c1 and returns the result on cOut.
func (evaluator *Evaluator) Add(c0, c1, cOut Operand) (err error) {

	if err = evaluateInPlace(c0.Element(), c1.Element(), cOut.Element(), evaluator.bfvcontext.contextQ.Add, evaluator.bfvcontext); err != nil {
		return err
	}

	return nil
}

// AddNoMod adds c0 to c1 without modular reduction, and returns the result on cOut.
func (evaluator *Evaluator) AddNoMod(c0, c1, cOut Operand) (err error) {

	if err = evaluateInPlace(c0.Element(), c1.Element(), cOut.Element(), evaluator.bfvcontext.contextQ.AddNoMod, evaluator.bfvcontext); err != nil {
		return err
	}

	return nil
}

// AddNew adds c0 to c1 and creates a new element to store the result.
func (evaluator *Evaluator) AddNew(c0, c1 Operand) (cOut *bfvElement, err error) {

	if cOut, err = evaluateNew(c0.Element(), c1.Element(), evaluator.bfvcontext.contextQ.Add, evaluator.bfvcontext); err != nil {
		return nil, err
	}

	return
}

// AddNoModNew adds c0 to c1 without modular reduction and creates a new element to store the result.
func (evaluator *Evaluator) AddNoModNew(c0, c1 Operand) (cOut *bfvElement, err error) {

	if cOut, err = evaluateNew(c0.Element(), c1.Element(), evaluator.bfvcontext.contextQ.AddNoMod, evaluator.bfvcontext); err != nil {
		return nil, err
	}

	return
}

// Sub subtracts c1 to c0 and returns the result on cOut.
func (evaluator *Evaluator) Sub(c0, c1, cOut Operand) (err error) {

	if err = evaluateInPlace(c0.Element(), c1.Element(), cOut.Element(), evaluator.bfvcontext.contextQ.Sub, evaluator.bfvcontext); err != nil {
		return err
	}

	return nil
}

// SubNoMod subtracts c1 to c0 without modular reduction and returns the result on cOut.
func (evaluator *Evaluator) SubNoMod(op0, op1, opOut Operand) (err error) {

	if err = evaluateInPlace(op0.Element(), op1.Element(), opOut.Element(), evaluator.bfvcontext.contextQ.SubNoMod, evaluator.bfvcontext); err != nil {
		return err
	}

	return nil
}

// SubNew subtracts c1 to c0 and creates a new element to store the result.
func (evaluator *Evaluator) SubNew(op0, op1 Operand) (cOut *bfvElement, err error) {

	c0, c1 := op0.Element(), op1.Element()
	if cOut, err = evaluateNew(c0, c1, evaluator.bfvcontext.contextQ.Sub, evaluator.bfvcontext); err != nil {
		return nil, err
	}

	return
}

// SubNoModNew subtracts c1 to c0 without modular reduction and creates a new element to store the result.
func (evaluator *Evaluator) SubNoModNew(op0, op1 Operand) (cOut *bfvElement, err error) {

	if cOut, err = evaluateNew(op0.Element(), op1.Element(), evaluator.bfvcontext.contextQ.SubNoMod, evaluator.bfvcontext); err != nil {
		return nil, err
	}

	return
}

// evaluateInPlace applies the provided function in place on c0 and c1 and returns the result in cOut
func evaluateInPlace(c0, c1, cOut *bfvElement, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly), bfvcontext *BfvContext) (err error) {

	maxDegree := max([]uint64{c0.Degree(), c1.Degree()})
	minDegree := min([]uint64{c0.Degree(), c1.Degree()})

	// Checks the validity of the receiver element
	if cOut.Degree() == 0 && cOut.Degree() < maxDegree {
		return errors.New("cannot evaluate(c0, c1 cOut) -> cOut is a plaintext (or an invalid ciphertext of degree 0) while c1 and/or c2 are ciphertexts of degree >= 1")
	} else {
		// Else resizes the receiver element
		cOut.Resize(bfvcontext, maxDegree)
	}

	for i := uint64(0); i < minDegree+1; i++ {
		evaluate(c0.Value()[i], c1.Value()[i], cOut.Value()[i])
	}

	// If the inputs degree differ, copies the remaining degree on the receiver
	// Also checks that the receiver is ont one of the inputs to avoid unnecessary work.

	if c0.Degree() > c1.Degree() && c0 != cOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			c0.Value()[i].Copy(cOut.Value()[i])
		}
	} else if c1.Degree() > c0.Degree() && c1 != cOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			c1.Value()[i].Copy(cOut.Value()[i])
		}
	}

	return nil
}

// evaluateNew applies the provided function on c0 and c1 and returns the result on a new element cOut.
func evaluateNew(c0, c1 *bfvElement, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly), bfvcontext *BfvContext) (cOut *bfvElement, err error) {

	if c0.Degree() >= c1.Degree() {

		cOut = c0.CopyNew()

		for i := range c1.Value() {
			evaluate(cOut.Value()[i], c1.Value()[i], cOut.Value()[i])
		}

	} else {

		cOut = c1.CopyNew()

		for i := range c0.Value() {
			evaluate(cOut.Value()[i], c0.Value()[i], cOut.Value()[i])
		}
	}

	return cOut, nil
}

// Neg negates c0 and returns the result on cOut.
func (evaluator *Evaluator) Neg(op0, opOut Operand) error {

	c0, cOut := op0.Element(), opOut.Element()

	if c0.Degree() != cOut.Degree() {
		return errors.New("error : invalid receiver ciphertext (degree not equal to input ciphertext")
	}

	for i := range c0.Value() {
		evaluator.bfvcontext.contextQ.Neg(c0.Value()[i], cOut.Value()[i])
	}

	return nil
}

// Neg negates c0 and creates a new element to store the result.
func (evaluator *Evaluator) NegNew(c0 *bfvElement) (cOut *bfvElement) {

	cOut = evaluator.bfvcontext.NewBfvElement(c0.Degree())

	for i := range c0.Value() {
		evaluator.bfvcontext.contextQ.Neg(c0.Value()[i], cOut.Value()[i])
	}

	return nil
}

// Reduce applies a modular reduction on c0 and returns the result on cOut.
func (evaluator *Evaluator) Reduce(op, opOut Operand) error {

	c0, cOut := op.Element(), opOut.Element()

	if c0.Degree() != cOut.Degree() {
		return errors.New("error : invalide ciphertext receiver (degree doesn't match c0.Degree")
	}

	for i := range c0.Value() {
		evaluator.bfvcontext.contextQ.Reduce(c0.Value()[i], cOut.Value()[i])
	}

	return nil
}

// Reduce applies a modular reduction on c0 and creates a new element to store the result.
func (evaluator *Evaluator) ReduceNew(c0 *bfvElement) (cOut *bfvElement) {

	cOut = evaluator.bfvcontext.NewBfvElement(c0.Degree())

	evaluator.Reduce(c0, cOut)

	return
}

// MulScalar multiplies c0 by an uint64 scalar and returns the result on cOut.
func (evaluator *Evaluator) MulScalar(c0 *bfvElement, scalar uint64, cOut *bfvElement) error {

	if c0.Degree() != cOut.Degree() {
		return errors.New("error : invalide ciphertext receiver (degree doesn't match c0.Degree")
	}

	for i := range c0.Value() {
		evaluator.bfvcontext.contextQ.MulScalar(c0.Value()[i], scalar, cOut.Value()[i])
	}

	return nil
}

// MulScalarNew multiplies c0 by an uint64 scalar and creates a new element to store the result.
func (evaluator *Evaluator) MulScalarNew(c0 *bfvElement, scalar uint64) (cOut *bfvElement) {

	cOut = evaluator.bfvcontext.NewBfvElement(c0.Degree())

	evaluator.MulScalar(c0, scalar, cOut)

	return
}

// Mul multiplies c0 by c1 and returns the result on cOut.
func (evaluator *Evaluator) Mul(op0, op1, opOut Operand) (err error) {

	el0, el1, elOut := op0.Element(), op1.Element(), opOut.Element()

	if elOut.Degree() < el0.Degree()+el1.Degree() {
		return errors.New("cannot Mul -> opOut (receiver) degree is to small to store the result")
	}

	tensorAndRescale(evaluator, el0, el1, elOut)

	return nil
}

// MulNew multiplies c0 by c1 and creates a new element to store the result.
func (evaluator *Evaluator) MulNew(op0, op1 Operand) (opOut *bfvElement) {

	el0, el1 := op0.Element(), op1.Element()

	opOut = evaluator.bfvcontext.NewBfvElement(el0.Degree() + el1.Degree())

	tensorAndRescale(evaluator, el0, el1, opOut.Element())

	return opOut
}

// tensorAndRescales computes (ct0 x ct1) * (t/Q) and stores the result on cOut.
func tensorAndRescale(evaluator *Evaluator, ct0, ct1, cOut *bfvElement) {

	// Prepares the ciphertexts for the Tensoring by extending their
	// basis from Q to QP and transforming them in NTT form
	c0 := evaluator.ctxpool[0]
	c1 := evaluator.ctxpool[1]
	tmpCout := evaluator.ctxpool[2]

	if ct0 == ct1 {

		for i := range ct0.Value() {
			evaluator.basisextender.ExtendBasis(ct0.Value()[i], c0.Value()[i])
			evaluator.bfvcontext.contextQP.NTT(c0.Value()[i], c0.Value()[i])
		}

	} else {

		for i := range ct0.Value() {
			evaluator.basisextender.ExtendBasis(ct0.Value()[i], c0.Value()[i])
			evaluator.bfvcontext.contextQP.NTT(c0.Value()[i], c0.Value()[i])
		}

		for i := range ct1.Value() {
			evaluator.basisextender.ExtendBasis(ct1.Value()[i], c1.Value()[i])
			evaluator.bfvcontext.contextQP.NTT(c1.Value()[i], c1.Value()[i])
		}
	}

	// Tensoring : multiplies each elements of the ciphertexts together
	// and adds them to their correspongint position in the new ciphertext
	// based on their respective degree

	// Case where both BfvElements are of degree 1
	if ct0.Degree() == 1 && ct1.Degree() == 1 {

		c_00 := evaluator.polypool[0]
		c_01 := evaluator.polypool[1]

		d0 := tmpCout.Value()[0]
		d1 := tmpCout.Value()[1]
		d2 := tmpCout.Value()[2]

		evaluator.bfvcontext.contextQP.MForm(c0.Value()[0], c_00)
		evaluator.bfvcontext.contextQP.MForm(c0.Value()[1], c_01)

		// Squaring case
		if ct0 == ct1 {

			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00, c0.Value()[0], d0) // c0 = c0[0]*c0[0]
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00, c0.Value()[1], d1) // c1 = 2*c0[0]*0[1]
			evaluator.bfvcontext.contextQP.Add(d1, d1, d1)
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_01, c0.Value()[1], d2) // c2 = c0[1]*c0[1]

			// Normal case
		} else {

			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00, c1.Value()[0], d0) // c0 = c0[0]*c0[0]
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00, c1.Value()[1], d1)
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomeryAndAddNoMod(c_01, c1.Value()[0], d1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
			evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_01, c1.Value()[1], d2)            // c2 = c0[1]*c1[1]
		}

		// Case where both BfvElements are not of degree 1
	} else {

		for i := 0; i < len(ct0.Value())+len(ct1.Value()); i++ {
			tmpCout.Value()[i].Zero()
		}

		// Squaring case
		if ct0 == ct1 {

			c_00 := evaluator.ctxpool[1]

			for i := range ct0.Value() {
				evaluator.bfvcontext.contextQP.MForm(c0.Value()[i], c_00.Value()[i])
			}

			for i := uint64(0); i < ct0.Degree()+1; i++ {
				for j := i + 1; j < ct0.Degree()+1; j++ {
					evaluator.bfvcontext.contextQP.MulCoeffsMontgomery(c_00.Value()[i], c0.Value()[j], tmpCout.Value()[i+j])
					evaluator.bfvcontext.contextQP.Add(tmpCout.Value()[i+j], tmpCout.Value()[i+j], tmpCout.Value()[i+j])
				}
			}

			for i := uint64(0); i < ct0.Degree()+1; i++ {
				evaluator.bfvcontext.contextQP.MulCoeffsMontgomeryAndAdd(c_00.Value()[i], c0.Value()[i], tmpCout.Value()[i<<1])
			}

			// Normal case
		} else {
			for i := range ct0.Value() {
				evaluator.bfvcontext.contextQP.MForm(c0.Value()[i], c0.Value()[i])
				for j := range ct1.Value() {
					evaluator.bfvcontext.contextQP.MulCoeffsMontgomeryAndAdd(c0.Value()[i], c1.Value()[j], tmpCout.Value()[i+j])
				}
			}
		}
	}

	// Applies the inverse NTT to the ciphertext, scales the down ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range cOut.Value() {
		evaluator.bfvcontext.contextQP.InvNTT(tmpCout.Value()[i], tmpCout.Value()[i])
		evaluator.complexscaler.Scale(tmpCout.Value()[i], cOut.Value()[i])
	}
}

// Relinearize relinearize the ciphertext cIn of degree > 1 until it is of degree 1 and returns the result on cOut.
//
// Requires a correct evaluation key as additional input :
//
// - it must match the secret-key that was used to create the public key under which the current cIn is encrypted.
//
// - it must be of degree high enough to relinearize the input ciphertext to degree 1 (ex. a ciphertext
//of degree 3 will require that the evaluation key stores the keys for both degree 3 and 2 ciphertexts).
func (evaluator *Evaluator) Relinearize(cIn *Ciphertext, evakey *EvaluationKey, cOut *Ciphertext) error {

	if int(cIn.Degree()-1) > len(evakey.evakey) {
		return errors.New("error : ciphertext degree too large to allow relinearization")
	}

	if cIn.Degree() < 2 {
		return errors.New("error : ciphertext is already of degree 1 or 0")
	}

	relinearize(evaluator, cIn, evakey, cOut)

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
func (evaluator *Evaluator) RelinearizeNew(cIn *Ciphertext, evakey *EvaluationKey) (cOut *Ciphertext, err error) {

	if int(cIn.Degree()-1) > len(evakey.evakey) {
		return nil, errors.New("error : ciphertext degree too large to allow relinearization")
	}

	if cIn.Degree() < 2 {
		return nil, errors.New("error : ciphertext is already of degree 1 or 0")
	}

	cOut = evaluator.bfvcontext.NewCiphertext(1)

	relinearize(evaluator, cIn, evakey, cOut)

	return cOut, nil
}

// relinearize is a methode common to Relinearize and RelinearizeNew. It switches cIn out in the NTT domain, applies the keyswitch, and returns the result out of the NTT domain.
func relinearize(evaluator *Evaluator, cIn *Ciphertext, evakey *EvaluationKey, cOut *Ciphertext) {

	evaluator.bfvcontext.contextQ.NTT(cIn.Value()[0], cOut.Value()[0])
	evaluator.bfvcontext.contextQ.NTT(cIn.Value()[1], cOut.Value()[1])

	for deg := uint64(cIn.Degree()); deg > 1; deg-- {
		switchKeys(evaluator, cOut.Value()[0], cOut.Value()[1], cIn.Value()[deg], evakey.evakey[deg-2], cOut.bfvElement)
	}

	if len(cOut.Value()) > 2 {
		cOut.SetValue(cOut.Value()[:2])
	}

	evaluator.bfvcontext.contextQ.InvNTT(cOut.Value()[0], cOut.Value()[0])
	evaluator.bfvcontext.contextQ.InvNTT(cOut.Value()[1], cOut.Value()[1])
}

// SwitchKeys applies the key-switching procedure to the ciphertext cIn and returns the result on cOut. It requires as an additional input a valide switching-key :
// it must encrypt the target key under the public key under which cIn is currently encrypted.
func (evaluator *Evaluator) SwitchKeys(cIn *Ciphertext, switchkey *SwitchingKey, cOut *Ciphertext) (err error) {

	if cIn.Degree() != 1 {
		return errors.New("error : ciphertext must be of degree 1 to allow key switching")
	}

	if cOut.Degree() != 1 {
		return errors.New("error : receiver ciphertext must be of degree 1 to allow key switching")
	}

	var c2 *ring.Poly
	if cIn == cOut {
		c2 = evaluator.polypool[1]
		cIn.Value()[1].Copy(c2)
	} else {
		c2 = cIn.Value()[1]
	}
	cIn.NTT(evaluator.bfvcontext, cOut.bfvElement)
	switchKeys(evaluator, cOut.Value()[0], cOut.Value()[1], c2, switchkey, cOut.bfvElement)
	cOut.InvNTT(evaluator.bfvcontext, cOut.bfvElement)

	return nil
}

// SwitchKeys applies the key-switching procedure to the ciphertext cIn and creates a new ciphertext to store the result. It requires as an additional input a valide switching-key :
// it must encrypt the target key under the public key under which cIn is currently encrypted.
func (evaluator *Evaluator) SwitchKeysNew(cIn *Ciphertext, switchkey *SwitchingKey) (cOut *Ciphertext, err error) {

	if cIn.Degree() > 1 {
		return nil, errors.New("error : ciphertext must be of degree 1 to allow key switching")
	}

	cOut = evaluator.bfvcontext.NewCiphertext(1)

	cIn.NTT(evaluator.bfvcontext, cOut.bfvElement)
	switchKeys(evaluator, cOut.Value()[0], cOut.Value()[1], cIn.Value()[1], switchkey, cOut.bfvElement)
	cOut.InvNTT(evaluator.bfvcontext, cOut.bfvElement)

	return cOut, nil
}

// RotateColumns rotates the columns of c0 by k position to the left and returns the result on c1. As an additional input it requires a rotationkeys :
//
// - it must either store all the left and right power of 2 rotations or the specific rotation that is asked.
//
// If only the power of two rotations are stored, the numbers k and n/2-k will be decomposed in base 2 and the rotation with the least
// hamming weight will be chosen, then the specific rotation will be computed as a sum of powers of two rotations.
func (evaluator *Evaluator) RotateColumns(op Operand, k uint64, evakey *RotationKeys, opOut Operand) (err error) {

	c0, c1 := op.Element(), opOut.Element()

	k &= ((evaluator.bfvcontext.n >> 1) - 1)

	if k == 0 {
		if c0 != c1 {
			if err = c0.Copy(c1); err != nil {
				return err
			}
		}

		return nil
	}

	if c1.Degree() != c0.Degree() {
		return errors.New("cannot rotate -> receiver degree doesn't match input degree ")
	}

	if c0.Degree() > 1 {
		return errors.New("cannot rotate -> input and or output degree not 0 or 1")
	}

	context := evaluator.bfvcontext.contextQ

	if c0.Degree() == 0 {

		if c1 != c0 {

			if c0.IsNTT() {
				ring.PermuteNTT(c0.Value()[0], evaluator.bfvcontext.galElRotColLeft[k], c1.Value()[0])
			} else {
				context.Permute(c0.Value()[0], evaluator.bfvcontext.galElRotColLeft[k], c1.Value()[0])
			}

		} else {

			if c0.IsNTT() {
				ring.PermuteNTT(c0.Value()[0], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[0])
			} else {
				context.Permute(c0.Value()[0], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[0])
			}

			context.Copy(evaluator.polypool[0], c1.Value()[0])
		}

		return nil

	} else {
		// Looks in the rotationkey if the corresponding rotation has been generated
		if evakey.evakey_rot_col_L[k] != nil {

			if c0.IsNTT() {

				ring.PermuteNTT(c0.Value()[0], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[0])
				ring.PermuteNTT(c0.Value()[1], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[1])

				context.Copy(evaluator.polypool[0], c1.Value()[0])
				context.Copy(evaluator.polypool[1], c1.Value()[1])

				context.InvNTT(evaluator.polypool[1], evaluator.polypool[1])

				switchKeys(evaluator, c1.Value()[0], c1.Value()[1], evaluator.polypool[1], evakey.evakey_rot_col_L[k], c1)

			} else {

				context.Permute(c0.Value()[0], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[0])
				context.Permute(c0.Value()[1], evaluator.bfvcontext.galElRotColLeft[k], evaluator.polypool[1])

				context.NTT(evaluator.polypool[0], c1.Value()[0])
				context.NTT(evaluator.polypool[1], c1.Value()[1])

				switchKeys(evaluator, c1.Value()[0], c1.Value()[1], evaluator.polypool[1], evakey.evakey_rot_col_L[k], c1)

				context.InvNTT(c1.Value()[0], c1.Value()[0])
				context.InvNTT(c1.Value()[1], c1.Value()[1])
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
					rotateColumnsLPow2(evaluator, c0, k, evakey, c1)
				} else {
					rotateColumnsRPow2(evaluator, c0, (evaluator.bfvcontext.n>>1)-k, evakey, c1)
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
func rotateColumnsLPow2(evaluator *Evaluator, c0 *bfvElement, k uint64, evakey *RotationKeys, c1 *bfvElement) {
	rotateColumnsPow2(evaluator, c0, evaluator.bfvcontext.gen, k, evakey.evakey_rot_col_L, c1)
}

// rotateColumnsRPow2 applies the Galois Endomorphism on the element, rotating the element by k positions to the right.
func rotateColumnsRPow2(evaluator *Evaluator, c0 *bfvElement, k uint64, evakey *RotationKeys, c1 *bfvElement) {
	rotateColumnsPow2(evaluator, c0, evaluator.bfvcontext.genInv, k, evakey.evakey_rot_col_R, c1)
}

// rotateColumnsPow2 rotates c0 by k position (left or right depending on the input), decomposing k as a sum of power of 2 rotations, and returns the result on c1.
func rotateColumnsPow2(evaluator *Evaluator, c0 *bfvElement, generator, k uint64, evakey_rot_col map[uint64]*SwitchingKey, c1 *bfvElement) {

	var mask, evakey_index uint64

	context := evaluator.bfvcontext.contextQ

	mask = (evaluator.bfvcontext.n << 1) - 1

	evakey_index = 1

	if c0.IsNTT() {
		c0.Copy(c1)
	} else {
		for i := range c0.Value() {
			context.NTT(c0.Value()[i], c1.Value()[i])
		}
	}

	// Applies the galois automorphism and the switching-key process
	for k > 0 {

		if k&1 == 1 {

			if c0.Degree() == 0 {

				ring.PermuteNTT(c1.Value()[0], generator, evaluator.polypool[0])
				context.Copy(evaluator.polypool[0], c1.Value()[0])

			} else {

				ring.PermuteNTT(c1.Value()[0], generator, evaluator.polypool[0])
				ring.PermuteNTT(c1.Value()[1], generator, evaluator.polypool[1])

				context.Copy(evaluator.polypool[0], c1.Value()[0])
				context.Copy(evaluator.polypool[1], c1.Value()[1])
				context.InvNTT(evaluator.polypool[1], evaluator.polypool[2])

				switchKeys(evaluator, evaluator.polypool[0], evaluator.polypool[1], evaluator.polypool[2], evakey_rot_col[evakey_index], c1)
			}

		}

		generator *= generator
		generator &= mask

		evakey_index <<= 1
		k >>= 1
	}

	if !c0.IsNTT() {
		for i := range c1.Value() {
			context.InvNTT(c1.Value()[i], c1.Value()[i])
		}
	}
}

// RotateRows swaps the rows of c0 and returns the result on c1.
func (evaluator *Evaluator) RotateRows(op Operand, evakey *RotationKeys, opOut Operand) error {

	c0, c1 := op.Element(), opOut.Element()

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

				ring.PermuteNTT(c0.Value()[0], evaluator.bfvcontext.galElRotRow, c1.Value()[0])

			} else {

				ring.PermuteNTT(c0.Value()[0], evaluator.bfvcontext.galElRotRow, evaluator.polypool[0])
				context.Copy(evaluator.polypool[0], c1.Value()[0])
			}

		} else {

			if c0 != c1 {

				context.Permute(c0.Value()[0], evaluator.bfvcontext.galElRotRow, c1.Value()[0])

			} else {

				context.Permute(c0.Value()[0], evaluator.bfvcontext.galElRotRow, evaluator.polypool[0])
				context.Copy(evaluator.polypool[0], c1.Value()[0])
			}
		}

	} else {

		if c0.IsNTT() {

			if c0 != c1 {

				ring.PermuteNTT(c0.Value()[0], evaluator.bfvcontext.galElRotRow, c1.Value()[0])
				ring.PermuteNTT(c0.Value()[1], evaluator.bfvcontext.galElRotRow, c1.Value()[1])

				context.InvNTT(c1.Value()[1], evaluator.polypool[1])

				switchKeys(evaluator, c1.Value()[0], c1.Value()[1], evaluator.polypool[1], evakey.evakey_rot_row, c1)

			} else {

				ring.PermuteNTT(c0.Value()[0], evaluator.bfvcontext.galElRotRow, evaluator.polypool[0])
				ring.PermuteNTT(c0.Value()[1], evaluator.bfvcontext.galElRotRow, evaluator.polypool[1])

				context.Copy(evaluator.polypool[0], c1.Value()[0])
				context.Copy(evaluator.polypool[1], c1.Value()[1])

				context.InvNTT(evaluator.polypool[1], evaluator.polypool[1])

				switchKeys(evaluator, c1.Value()[0], c1.Value()[1], evaluator.polypool[1], evakey.evakey_rot_row, c1)
			}

		} else {

			if c0 != c1 {

				context.Permute(c0.Value()[0], evaluator.bfvcontext.galElRotRow, c1.Value()[0])
				context.Permute(c0.Value()[1], evaluator.bfvcontext.galElRotRow, c1.Value()[1])

				context.Copy(c1.Value()[1], evaluator.polypool[1])

				context.NTT(c1.Value()[0], c1.Value()[0])
				context.NTT(c1.Value()[1], c1.Value()[1])

				switchKeys(evaluator, c1.Value()[0], c1.Value()[1], evaluator.polypool[1], evakey.evakey_rot_row, c1)

				context.InvNTT(c1.Value()[0], c1.Value()[0])
				context.InvNTT(c1.Value()[1], c1.Value()[1])

			} else {

				context.Permute(c0.Value()[0], evaluator.bfvcontext.galElRotRow, evaluator.polypool[0])
				context.Permute(c0.Value()[1], evaluator.bfvcontext.galElRotRow, evaluator.polypool[1])

				context.NTT(evaluator.polypool[0], c1.Value()[0])
				context.NTT(evaluator.polypool[1], c1.Value()[1])

				switchKeys(evaluator, c1.Value()[0], c1.Value()[1], evaluator.polypool[1], evakey.evakey_rot_row, c1)

				context.InvNTT(c1.Value()[0], c1.Value()[0])
				context.InvNTT(c1.Value()[1], c1.Value()[1])
			}
		}
	}

	return nil
}

// InnerSum computs the inner sum of c0 and returns the result on c1. It requires a rotation key storing all the left power of two rotations.
// The resulting vector will be of the form [sum, sum, .., sum, sum ].
func (evaluator *Evaluator) InnerSum(c0 *bfvElement, evakey *RotationKeys, c1 *bfvElement) error {
	if c0.Degree() != 1 {
		return errors.New("error : ciphertext must be of degree 1 to allow Galois Auotomorphism (required for inner sum)")
	}

	if c1.Degree() != 1 {
		return errors.New("error : receiver ciphertext must be of degree 1 to allow Galois Automorphism (required for inner sum)")
	}

	cTmp := evaluator.bfvcontext.NewCiphertext(1)

	if c0 != c1 {
		if err := c1.Copy(c0); err != nil {
			return err
		}
	}

	for i := uint64(1); i < evaluator.bfvcontext.n>>1; i <<= 1 {
		if err := evaluator.RotateColumns(c1, i, evakey, cTmp.bfvElement); err != nil {
			return err
		}
		evaluator.Add(cTmp.bfvElement, c1, c1)
	}

	if err := evaluator.RotateRows(c1, evakey, cTmp.bfvElement); err != nil {
		return err
	}
	evaluator.Add(c1, cTmp.bfvElement, c1)

	return nil

}

// switchKeys compute cOut = [c0 + c2*evakey[0], c1 + c2*evakey[1]].
func switchKeys(evaluator *Evaluator, c0, c1, c2 *ring.Poly, evakey *SwitchingKey, cOut *bfvElement) {

	var mask, reduce, bitLog uint64

	c2_qi_w := evaluator.polypool[3]

	mask = uint64(((1 << evakey.bitDecomp) - 1))

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

			evaluator.bfvcontext.contextQ.MulCoeffsMontgomeryAndAddNoMod(evakey.evakey[i][j][0], c2_qi_w, cOut.Value()[0])
			evaluator.bfvcontext.contextQ.MulCoeffsMontgomeryAndAddNoMod(evakey.evakey[i][j][1], c2_qi_w, cOut.Value()[1])

			if reduce&7 == 7 {
				evaluator.bfvcontext.contextQ.Reduce(cOut.Value()[0], cOut.Value()[0])
				evaluator.bfvcontext.contextQ.Reduce(cOut.Value()[1], cOut.Value()[1])
			}

			reduce += 1
		}
	}

	if (reduce-1)&7 != 7 {
		evaluator.bfvcontext.contextQ.Reduce(cOut.Value()[0], cOut.Value()[0])
		evaluator.bfvcontext.contextQ.Reduce(cOut.Value()[1], cOut.Value()[1])
	}
}
