package rgsw

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// AddNoModLvl adds op to ctOut, without modular reduction.
func AddNoModLvl(levelQ, levelP int, op interface{}, ringQP ringqp.Ring, ctOut *Ciphertext) {
	switch el := op.(type) {
	case *Plaintext:

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
					ring.AddVecNoMod(ctOut.Value[0].Value[i][j][0].Q.Coeffs[k], el.Value[j].Coeffs[k], ctOut.Value[0].Value[i][j][0].Q.Coeffs[k])
					ring.AddVecNoMod(ctOut.Value[1].Value[i][j][1].Q.Coeffs[k], el.Value[j].Coeffs[k], ctOut.Value[1].Value[i][j][1].Q.Coeffs[k])
				}
			}
		}
	case *Ciphertext:
		for i := range el.Value[0].Value {
			for j := range el.Value[0].Value[i] {
				ringQP.AddNoModLvl(levelQ, levelP, ctOut.Value[0].Value[i][j][0], el.Value[0].Value[i][j][0], ctOut.Value[0].Value[i][j][0])
				ringQP.AddNoModLvl(levelQ, levelP, ctOut.Value[0].Value[i][j][1], el.Value[0].Value[i][j][1], ctOut.Value[0].Value[i][j][1])
				ringQP.AddNoModLvl(levelQ, levelP, ctOut.Value[1].Value[i][j][0], el.Value[1].Value[i][j][0], ctOut.Value[1].Value[i][j][0])
				ringQP.AddNoModLvl(levelQ, levelP, ctOut.Value[1].Value[i][j][1], el.Value[1].Value[i][j][1], ctOut.Value[1].Value[i][j][1])
			}
		}
	default:
		panic("unsuported op.(type), must be either *rgsw.Plaintext or *rgsw.Ciphertext")
	}
}

// ReduceLvl applies the modular reduction on ctIn and returns the result on ctOut.
func ReduceLvl(levelQ, levelP int, ctIn *Ciphertext, ringQP ringqp.Ring, ctOut *Ciphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.ReduceLvl(levelQ, levelP, ctIn.Value[0].Value[i][j][0], ctOut.Value[0].Value[i][j][0])
			ringQP.ReduceLvl(levelQ, levelP, ctIn.Value[0].Value[i][j][1], ctOut.Value[0].Value[i][j][1])
			ringQP.ReduceLvl(levelQ, levelP, ctIn.Value[1].Value[i][j][0], ctOut.Value[1].Value[i][j][0])
			ringQP.ReduceLvl(levelQ, levelP, ctIn.Value[1].Value[i][j][1], ctOut.Value[1].Value[i][j][1])
		}
	}
}

// MulByXPowAlphaMinusOneConstantLvl multiplies ctOut by (X^alpha - 1) and returns the result on ctOut.
func MulByXPowAlphaMinusOneConstantLvl(levelQ, levelP int, ctIn *Ciphertext, powXMinusOne ringqp.Poly, ringQP ringqp.Ring, ctOut *Ciphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ctIn.Value[0].Value[i][j][0], powXMinusOne, ctOut.Value[0].Value[i][j][0])
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ctIn.Value[0].Value[i][j][1], powXMinusOne, ctOut.Value[0].Value[i][j][1])
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ctIn.Value[1].Value[i][j][0], powXMinusOne, ctOut.Value[1].Value[i][j][0])
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ctIn.Value[1].Value[i][j][1], powXMinusOne, ctOut.Value[1].Value[i][j][1])
		}
	}
}

// MulByXPowAlphaMinusOneAndAddNoModLvl multiplies ctOut by (X^alpha - 1) and adds the result on ctOut.
func MulByXPowAlphaMinusOneAndAddNoModLvl(levelQ, levelP int, ctIn *Ciphertext, powXMinusOne ringqp.Poly, ringQP ringqp.Ring, ctOut *Ciphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ctIn.Value[0].Value[i][j][0], powXMinusOne, ctOut.Value[0].Value[i][j][0])
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ctIn.Value[0].Value[i][j][1], powXMinusOne, ctOut.Value[0].Value[i][j][1])
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ctIn.Value[1].Value[i][j][0], powXMinusOne, ctOut.Value[1].Value[i][j][0])
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ctIn.Value[1].Value[i][j][1], powXMinusOne, ctOut.Value[1].Value[i][j][1])
		}
	}
}
