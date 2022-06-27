package rlwe

import (
	"math"

	"github.com/tuneinsight/lattigo/v3/ring"
)

// KeySwitcher is a struct for RLWE key-switching.
type KeySwitcher struct {
	*Parameters
	*keySwitcherBuffer
	BasisExtender *ring.BasisExtender
	Decomposer    *ring.Decomposer
}

type keySwitcherBuffer struct {
	// BuffQ[0]/BuffP[0] : on the fly decomp(c2)
	// BuffQ[1-5]/BuffP[1-5] : available
	BuffQP       [6]PolyQP
	BuffInvNTT   *ring.Poly
	BuffDecompQP []PolyQP // Memory Buff for the basis extension in hoisting
}

func newKeySwitcherBuffer(params Parameters) *keySwitcherBuffer {

	buff := new(keySwitcherBuffer)
	beta := params.Beta()
	ringQP := params.RingQP()

	buff.BuffQP = [6]PolyQP{ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly()}

	buff.BuffInvNTT = params.RingQ().NewPoly()

	buff.BuffDecompQP = make([]PolyQP, beta)
	for i := 0; i < beta; i++ {
		buff.BuffDecompQP[i] = ringQP.NewPoly()
	}

	return buff
}

// NewKeySwitcher creates a new KeySwitcher.
func NewKeySwitcher(params Parameters) *KeySwitcher {
	ks := new(KeySwitcher)
	ks.Parameters = &params
	ks.BasisExtender = ring.NewBasisExtender(params.RingQ(), params.RingP())
	ks.Decomposer = ring.NewDecomposer(params.RingQ(), params.RingP())
	ks.keySwitcherBuffer = newKeySwitcherBuffer(params)
	return ks
}

// ShallowCopy creates a copy of a KeySwitcher, only reallocating the memory Buff.
func (ks *KeySwitcher) ShallowCopy() *KeySwitcher {
	return &KeySwitcher{
		Parameters:        ks.Parameters,
		Decomposer:        ks.Decomposer,
		keySwitcherBuffer: newKeySwitcherBuffer(*ks.Parameters),
		BasisExtender:     ks.BasisExtender.ShallowCopy(),
	}
}

// SwitchKeysInPlace applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
// Will return the result in the same NTT domain as the input cx.
func (ks *KeySwitcher) SwitchKeysInPlace(levelQ int, cx *ring.Poly, evakey *SwitchingKey, p0, p1 *ring.Poly) {
	ks.SwitchKeysInPlaceNoModDown(levelQ, cx, evakey, p0, ks.BuffQP[1].P, p1, ks.BuffQP[2].P)

	levelP := len(evakey.Value[0][0].P.Coeffs) - 1

	if cx.IsNTT {
		ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p0, ks.BuffQP[1].P, p0)
		ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p1, ks.BuffQP[2].P, p1)
	} else {

		ks.ringQ.InvNTTLazyLvl(levelQ, p0, p0)
		ks.ringQ.InvNTTLazyLvl(levelQ, p1, p1)
		ks.ringP.InvNTTLazyLvl(levelP, ks.BuffQP[1].P, ks.BuffQP[1].P)
		ks.ringP.InvNTTLazyLvl(levelP, ks.BuffQP[2].P, ks.BuffQP[2].P)

		ks.BasisExtender.ModDownQPtoQ(levelQ, levelP, p0, ks.BuffQP[1].P, p0)
		ks.BasisExtender.ModDownQPtoQ(levelQ, levelP, p1, ks.BuffQP[2].P, p1)
	}
}

// DecomposeNTT applies the full RNS basis decomposition for all q_alpha_i on c2.
// Expects the IsNTT flag of c2 to correctly reflect the domain of c2.
// BuffDecompQ and BuffDecompQ are vectors of polynomials (mod Q and mod P) that store the
// special RNS decomposition of c2 (in the NTT domain)
func (ks *KeySwitcher) DecomposeNTT(levelQ, levelP, alpha int, c2 *ring.Poly, BuffDecomp []PolyQP) {

	ringQ := ks.RingQ()

	var polyNTT, polyInvNTT *ring.Poly

	if c2.IsNTT {
		polyNTT = c2
		polyInvNTT = ks.BuffInvNTT
		ringQ.InvNTTLvl(levelQ, polyNTT, polyInvNTT)
	} else {
		polyNTT = ks.BuffInvNTT
		polyInvNTT = c2
		ringQ.NTTLvl(levelQ, polyInvNTT, polyNTT)
	}

	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	for i := 0; i < beta; i++ {
		ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, polyNTT, polyInvNTT, BuffDecomp[i].Q, BuffDecomp[i].P)
	}
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo q_alpha_beta, and returns the result on c2QiQ are c2QiP the receiver polynomials
// respectively mod Q and mod P (in the NTT domain)
func (ks *KeySwitcher) DecomposeSingleNTT(levelQ, levelP, alpha, beta int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := ks.RingQ()
	ringP := ks.RingP()

	ks.Decomposer.DecomposeAndSplit(levelQ, levelP, alpha, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * (levelP + 1)
	p0idxed := p0idxst + 1

	// c2_qi = cx mod qi mod qi
	for x := 0; x < levelQ+1; x++ {
		if p0idxst <= x && x < p0idxed {
			copy(c2QiQ.Coeffs[x], c2NTT.Coeffs[x])
		} else {
			ringQ.NTTSingleLazy(x, c2QiQ.Coeffs[x], c2QiQ.Coeffs[x])
		}
	}
	// c2QiP = c2 mod qi mod pj
	ringP.NTTLazyLvl(levelP, c2QiP, c2QiP)
}

// SwitchKeysInPlaceNoModDown applies the key-switch to the polynomial cx :
//
// Buff2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// Buff3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (ks *KeySwitcher) SwitchKeysInPlaceNoModDown(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

	var reduce int

	ringQ := ks.RingQ()
	ringP := ks.RingP()
	ringQP := ks.RingQP()

	c2QP := ks.BuffQP[0]

	var cxNTT, cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxNTT = cx
		cxInvNTT = ks.BuffInvNTT
		ringQ.InvNTTLvl(levelQ, cxNTT, cxInvNTT)
	} else {
		cxNTT = ks.BuffInvNTT
		cxInvNTT = cx
		ringQ.NTTLvl(levelQ, cxInvNTT, cxNTT)
	}

	c0QP := PolyQP{c0Q, c0P}
	c1QP := PolyQP{c1Q, c1P}

	reduce = 0

	alpha := len(evakey.Value[0][0].P.Coeffs)
	levelP := alpha - 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	for i := 0; i < beta; i++ {

		ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][1], c2QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][1], c2QP, c1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}
}

// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (BuffDecompQ and BuffDecompP)
// and divides the result by P, reducing the basis from QP to Q.
//
// Buff2 = dot(BuffDecompQ||BuffDecompP * evakey[0]) mod Q
// Buff3 = dot(BuffDecompQ||BuffDecompP * evakey[1]) mod Q
func (ks *KeySwitcher) KeyswitchHoisted(levelQ int, BuffDecompQP []PolyQP, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	ks.KeyswitchHoistedNoModDown(levelQ, BuffDecompQP, evakey, c0Q, c1Q, c0P, c1P)

	levelP := len(evakey.Value[0][0].P.Coeffs) - 1

	// Computes c0Q = c0Q/c0P and c1Q = c1Q/c1P
	ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0Q, c0P, c0Q)
	ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1Q, c1P, c1Q)
}

// KeyswitchHoistedNoModDown applies the key-switch to the decomposed polynomial c2 mod QP (BuffDecompQ and BuffDecompP)
//
// Buff2 = dot(BuffDecompQ||BuffDecompP * evakey[0]) mod QP
// Buff3 = dot(BuffDecompQ||BuffDecompP * evakey[1]) mod QP
func (ks *KeySwitcher) KeyswitchHoistedNoModDown(levelQ int, BuffDecompQP []PolyQP, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	ringQ := ks.RingQ()
	ringP := ks.RingP()
	ringQP := ks.RingQP()

	c0QP := PolyQP{c0Q, c0P}
	c1QP := PolyQP{c1Q, c1P}

	alpha := len(evakey.Value[0][0].P.Coeffs)
	levelP := alpha - 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(alpha)))

	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < beta; i++ {

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0], BuffDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][1], BuffDecompQP[i], c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0], BuffDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][1], BuffDecompQP[i], c1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}
}
