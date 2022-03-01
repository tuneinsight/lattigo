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
<<<<<<< dev_bfv_poly
	// BuffQ[0]/BuffP[0] : on the fly decomp(c2)
	// BuffQ[1-5]/BuffP[1-5] : available
	BuffQP       [6]PolyQP
	BuffInvNTT   *ring.Poly
	BuffDecompQP []PolyQP // Memory Buff for the basis extension in hoisting
=======
	// PoolQ[0]/PoolP[0] : on the fly decomp(c2)
	// PoolQ[1-5]/PoolP[1-5] : available
	Pool         [6]PolyQP
	PoolInvNTT   *ring.Poly
	PoolDecompQP []PolyQP // Memory pool for the basis extension in hoisting
	PoolBitDecomp []uint64
>>>>>>> wip
}

func newKeySwitcherBuffer(params Parameters) *keySwitcherBuffer {

	buff := new(keySwitcherBuffer)
	decompRNS := params.DecompRNS(params.QCount()-1, params.PCount()-1)
	ringQP := params.RingQP()

	buff.BuffQP = [6]PolyQP{ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly()}

	buff.BuffInvNTT = params.RingQ().NewPoly()

<<<<<<< dev_bfv_poly
	buff.BuffDecompQP = make([]PolyQP, beta)
	for i := 0; i < beta; i++ {
		buff.BuffDecompQP[i] = ringQP.NewPoly()
=======
	buff.PoolDecompQP = make([]PolyQP, decompRNS)
	for i := 0; i < decompRNS; i++ {
		buff.PoolDecompQP[i] = ringQP.NewPoly()
>>>>>>> First step for adding bit-decomp
	}

	buff.PoolBitDecomp = make([]uint64, params.RingQ().N)

	return buff
}

// NewKeySwitcher creates a new KeySwitcher.
func NewKeySwitcher(params Parameters) *KeySwitcher {
	ks := new(KeySwitcher)
	ks.Parameters = &params
	if params.RingP() != nil{
		ks.BasisExtender = ring.NewBasisExtender(params.RingQ(), params.RingP())
		ks.Decomposer = ring.NewDecomposer(params.RingQ(), params.RingP())
	}
	
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
<<<<<<< dev_bfv_poly
	ks.SwitchKeysInPlaceNoModDown(levelQ, cx, evakey, p0, ks.BuffQP[1].P, p1, ks.BuffQP[2].P)
=======
>>>>>>> wip

	levelP := evakey.LevelP()

<<<<<<< dev_bfv_poly
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
=======
	if levelP > 0{
		ks.SwitchKeysInPlaceNoModDown(levelQ, cx, evakey, p0, ks.Pool[1].P, p1, ks.Pool[2].P)
	}else{
		ks.SwitchKeyInPlaceSinglePAndBitDecomp(levelQ, cx, evakey, p0, ks.Pool[1].P, p1, ks.Pool[2].P)
	}

	if cx.IsNTT && levelP != -1{
		ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p0, ks.Pool[1].P, p0)
		ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p1, ks.Pool[2].P, p1)
	} else if !cx.IsNTT {
		ks.ringQ.InvNTTLazyLvl(levelQ, p0, p0)
		ks.ringQ.InvNTTLazyLvl(levelQ, p1, p1)

		if levelP != -1{
			ks.ringP.InvNTTLazyLvl(levelP, ks.Pool[1].P, ks.Pool[1].P)
			ks.ringP.InvNTTLazyLvl(levelP, ks.Pool[2].P, ks.Pool[2].P)
			ks.BasisExtender.ModDownQPtoQ(levelQ, levelP, p0, ks.Pool[1].P, p0)
			ks.BasisExtender.ModDownQPtoQ(levelQ, levelP, p1, ks.Pool[2].P, p1)
		}
>>>>>>> wip
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

	p0idxst := beta * alpha
	p0idxed := p0idxst + 1

	// c2_qi = cx mod qi mod qi
	for x := 0; x < levelQ+1; x++ {
		if p0idxst <= x && x < p0idxed {
			copy(c2QiQ.Coeffs[x], c2NTT.Coeffs[x])
		} else {
			ringQ.NTTSingle(x, c2QiQ.Coeffs[x], c2QiQ.Coeffs[x])
		}
	}

	if ringP != nil{
		// c2QiP = c2 mod qi mod pj
		ringP.NTTLvl(levelP, c2QiP, c2QiP)
	}
}

// SwitchKeysInPlaceNoModDown applies the key-switch to the polynomial cx :
//
// Buff2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// Buff3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (ks *KeySwitcher) SwitchKeysInPlaceNoModDown(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

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

	levelP := evakey.Value[0][0][0].P.Level()
	decompRNS := (levelQ + 1 + levelP) / (levelP + 1)

	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		ks.DecomposeSingleNTT(levelQ, levelP, levelP+1, i, cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][1], c2QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][1], c2QP, c1QP)
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

<<<<<<< dev_bfv_poly
// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (BuffDecompQ and BuffDecompP)
=======
// SwitchKeyInPlaceSinglePAndBitDecomp applies the key-switch to the polynomial cx :
//
// pool2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// pool3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (ks *KeySwitcher) SwitchKeyInPlaceSinglePAndBitDecomp(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

	ringQ := ks.RingQ()
	ringP := ks.RingP()

	var cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxInvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, cx, cxInvNTT)
	} else {
		cxInvNTT = cx
	}

	c0QP := PolyQP{c0Q, c0P}
	c1QP := PolyQP{c1Q, c1P}

	var levelP int
	if evakey.Value[0][0][0].P != nil{
		levelP = evakey.Value[0][0][0].P.Level()
	}else{
		levelP = -1
	}
	
	decompRNS := ks.DecompRNS(levelQ, levelP)
	decompBIT := ks.DecompBIT(levelQ, levelP)

	lb2 := ks.logbase2

	mask := uint64(((1<<lb2)-1))

	if mask == 0{
		mask = 0xFFFFFFFFFFFFFFFF
	}

	cw := ks.Pool[0].Q.Coeffs[0]
	cwNTT := ks.PoolBitDecomp

	QiOverF := ks.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {
		for j := 0; j < decompBIT; j++{

			ring.MaskVec(cxInvNTT.Coeffs[i], cw, j*lb2, mask)

			if i == 0 && j == 0{
				for u := 0; u < levelQ+1; u++{
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++{
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			} else {
				for u := 0; u < levelQ+1; u++{
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++{
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
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




// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (PoolDecompQ and PoolDecompP)
>>>>>>> wip
// and divides the result by P, reducing the basis from QP to Q.
//
// Buff2 = dot(BuffDecompQ||BuffDecompP * evakey[0]) mod Q
// Buff3 = dot(BuffDecompQ||BuffDecompP * evakey[1]) mod Q
func (ks *KeySwitcher) KeyswitchHoisted(levelQ int, BuffDecompQP []PolyQP, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	ks.KeyswitchHoistedNoModDown(levelQ, BuffDecompQP, evakey, c0Q, c1Q, c0P, c1P)

	levelP := evakey.Value[0][0][0].P.Level()

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

	levelP := evakey.Value[0][0][0].P.Level()
	decompRNS := (levelQ + 1 + levelP) / (levelP + 1)

	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		if i == 0 {
<<<<<<< dev_bfv_poly
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0], BuffDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][1], BuffDecompQP[i], c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0], BuffDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][1], BuffDecompQP[i], c1QP)
=======
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][0], PoolDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][1], PoolDecompQP[i], c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][0], PoolDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][1], PoolDecompQP[i], c1QP)
>>>>>>> First step for adding bit-decomp
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
