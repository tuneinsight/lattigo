package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// Ciphertext is a generic type for RGSW ciphertext.
type RGSWCiphertext struct {
	Value [2]GadgetCiphertext
}

// LevelQ returns the level of the modulus Q of the target.
func (ct *RGSWCiphertext) LevelQ() int {
	return ct.Value[0].LevelQ()
}

// LevelP returns the level of the modulus P of the target.
func (ct *RGSWCiphertext) LevelP() int {
	return ct.Value[0].LevelP()
}

// NewCiphertext allocates a new RGSW ciphertext in the NTT domain.
func NewRGSWCiphertext(levelQ, levelP, decompRNS, decompBit int, ringQP ringqp.Ring) (ct *RGSWCiphertext) {
	return &RGSWCiphertext{
		Value: [2]GadgetCiphertext{
			*NewGadgetCiphertext(levelQ, levelP, decompRNS, decompBit, ringQP),
			*NewGadgetCiphertext(levelQ, levelP, decompRNS, decompBit, ringQP),
		},
	}
}

// Plaintext stores an RGSW plaintext value.
type RGSWPlaintext GadgetPlaintext

// NewPlaintext creates a new RGSW plaintext from value, which can be either uint64, int64 or *ring.Poly.
// Plaintext is returned in the NTT and Mongtomery domain.
func NewRGSWPlaintext(value interface{}, levelQ, levelP, logBase2, decompBIT int, ringQP ringqp.Ring) (pt *RGSWPlaintext) {
	return &RGSWPlaintext{Value: NewGadgetPlaintext(value, levelQ, levelP, logBase2, decompBIT, ringQP).Value}
}

// AddNoModLvl adds op to ctOut, without modular reduction.
func AddNoModLvl(levelQ, levelP int, op interface{}, ringQP ringqp.Ring, ctOut *RGSWCiphertext) {
	switch el := op.(type) {
	case *RGSWPlaintext:

		nQ := levelQ + 1
		nP := levelP + 1

		if nP == 0 {
			nP = 1
		}

		for i := range ctOut.Value[0].Value {
			for j := range ctOut.Value[0].Value[i] {
				start, end := i*nP, (i+1)*nP
				if end > nQ {
					end = nQ
				}
				for k := start; k < end; k++ {
					ring.AddVecNoMod(ctOut.Value[0].Value[i][j].Value[0].Q.Coeffs[k], el.Value[j].Coeffs[k], ctOut.Value[0].Value[i][j].Value[0].Q.Coeffs[k])
					ring.AddVecNoMod(ctOut.Value[1].Value[i][j].Value[1].Q.Coeffs[k], el.Value[j].Coeffs[k], ctOut.Value[1].Value[i][j].Value[1].Q.Coeffs[k])
				}
			}
		}
	case *RGSWCiphertext:
		for i := range el.Value[0].Value {
			for j := range el.Value[0].Value[i] {
				ringQP.AddNoModLvl(levelQ, levelP, ctOut.Value[0].Value[i][j].Value[0], el.Value[0].Value[i][j].Value[0], ctOut.Value[0].Value[i][j].Value[0])
				ringQP.AddNoModLvl(levelQ, levelP, ctOut.Value[0].Value[i][j].Value[1], el.Value[0].Value[i][j].Value[1], ctOut.Value[0].Value[i][j].Value[1])
				ringQP.AddNoModLvl(levelQ, levelP, ctOut.Value[1].Value[i][j].Value[0], el.Value[1].Value[i][j].Value[0], ctOut.Value[1].Value[i][j].Value[0])
				ringQP.AddNoModLvl(levelQ, levelP, ctOut.Value[1].Value[i][j].Value[1], el.Value[1].Value[i][j].Value[1], ctOut.Value[1].Value[i][j].Value[1])
			}
		}
	default:
		panic("unsuported op.(type), must be either *rgsw.Plaintext or *rgsw.Ciphertext")
	}
}

// ReduceLvl applies the modular reduction on ctIn and returns the result on ctOut.
func ReduceLvl(levelQ, levelP int, ctIn *RGSWCiphertext, ringQP ringqp.Ring, ctOut *RGSWCiphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.ReduceLvl(levelQ, levelP, ctIn.Value[0].Value[i][j].Value[0], ctOut.Value[0].Value[i][j].Value[0])
			ringQP.ReduceLvl(levelQ, levelP, ctIn.Value[0].Value[i][j].Value[1], ctOut.Value[0].Value[i][j].Value[1])
			ringQP.ReduceLvl(levelQ, levelP, ctIn.Value[1].Value[i][j].Value[0], ctOut.Value[1].Value[i][j].Value[0])
			ringQP.ReduceLvl(levelQ, levelP, ctIn.Value[1].Value[i][j].Value[1], ctOut.Value[1].Value[i][j].Value[1])
		}
	}
}

// MulByXPowAlphaMinusOneConstantLvl multiplies ctOut by (X^alpha - 1) and returns the result on ctOut.
func MulByXPowAlphaMinusOneConstantLvl(levelQ, levelP int, ctIn *RGSWCiphertext, powXMinusOne ringqp.Poly, ringQP ringqp.Ring, ctOut *RGSWCiphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ctIn.Value[0].Value[i][j].Value[0], powXMinusOne, ctOut.Value[0].Value[i][j].Value[0])
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ctIn.Value[0].Value[i][j].Value[1], powXMinusOne, ctOut.Value[0].Value[i][j].Value[1])
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ctIn.Value[1].Value[i][j].Value[0], powXMinusOne, ctOut.Value[1].Value[i][j].Value[0])
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ctIn.Value[1].Value[i][j].Value[1], powXMinusOne, ctOut.Value[1].Value[i][j].Value[1])
		}
	}
}

// MulByXPowAlphaMinusOneAndAddNoModLvl multiplies ctOut by (X^alpha - 1) and adds the result on ctOut.
func MulByXPowAlphaMinusOneAndAddNoModLvl(levelQ, levelP int, ctIn *RGSWCiphertext, powXMinusOne ringqp.Poly, ringQP ringqp.Ring, ctOut *RGSWCiphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ctIn.Value[0].Value[i][j].Value[0], powXMinusOne, ctOut.Value[0].Value[i][j].Value[0])
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ctIn.Value[0].Value[i][j].Value[1], powXMinusOne, ctOut.Value[0].Value[i][j].Value[1])
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ctIn.Value[1].Value[i][j].Value[0], powXMinusOne, ctOut.Value[1].Value[i][j].Value[0])
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ctIn.Value[1].Value[i][j].Value[1], powXMinusOne, ctOut.Value[1].Value[i][j].Value[1])
		}
	}
}
