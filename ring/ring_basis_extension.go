package ring

import (
	"math"
	"math/big"
	"unsafe"
)

// FastBasisExtender stores the necessary parameters for RNS basis extension.
// The algorithm is from https://eprint.iacr.org/2018/117.pdf.
type FastBasisExtender struct {
	contextQ        *Context
	contextP        *Context
	paramsQP        *modupParams
	paramsPQ        *modupParams
	modDownParamsPQ []uint64
	modDownParamsQP []uint64
	polypoolQ       *Poly
	polypoolP       *Poly
}

type modupParams struct {
	Q []uint64
	P []uint64

	//Parameters for basis extension from Q to P
	// (Q/Qi)^-1) (mod each Qi) (in Montgomery form)
	qibMont []uint64
	// Q/qi (mod each Pj) (in Montgomery form)
	qispjMont [][]uint64
	// Q*v (mod each Pj) for v in [1,...,k] where k is the number of Pj moduli
	qpjInv [][]uint64

	bredParamsQ [][]uint64
	mredParamsQ []uint64

	bredParamsP [][]uint64
	mredParamsP []uint64
}

func genModDownParams(contextP, contextQ *Context) (params []uint64) {

	params = make([]uint64, len(contextP.Modulus))

	bredParams := contextP.GetBredParams()
	tmp := new(big.Int)
	for i, Qi := range contextP.Modulus {

		params[i] = tmp.Mod(contextQ.ModulusBigint, NewUint(Qi)).Uint64()
		params[i] = ModExp(params[i], Qi-2, Qi)
		params[i] = MForm(params[i], Qi, bredParams[i])
	}

	return
}

// NewFastBasisExtender creates a new FastBasisExtender, enabling RNS basis extension from Q to P and P to Q.
func NewFastBasisExtender(contextQ, contextP *Context) *FastBasisExtender {

	newParams := new(FastBasisExtender)

	newParams.contextQ = contextQ
	newParams.contextP = contextP

	newParams.paramsQP = basisextenderparameters(contextQ.Modulus, contextP.Modulus)
	newParams.paramsPQ = basisextenderparameters(contextP.Modulus, contextQ.Modulus)

	newParams.modDownParamsPQ = genModDownParams(contextQ, contextP)
	newParams.modDownParamsQP = genModDownParams(contextP, contextQ)

	newParams.polypoolQ = contextQ.NewPoly()
	newParams.polypoolP = contextP.NewPoly()

	return newParams
}

func basisextenderparameters(Q, P []uint64) (params *modupParams) {

	params = new(modupParams)

	params.Q = make([]uint64, len(Q))
	params.bredParamsQ = make([][]uint64, len(Q))
	params.mredParamsQ = make([]uint64, len(Q))
	for i, qi := range Q {
		params.Q[i] = Q[i]
		params.bredParamsQ[i] = BRedParams(qi)
		params.mredParamsQ[i] = MRedParams(qi)
	}

	params.P = make([]uint64, len(P))
	params.bredParamsP = make([][]uint64, len(P))
	params.mredParamsP = make([]uint64, len(P))
	for i, pj := range P {
		params.P[i] = P[i]
		params.bredParamsP[i] = BRedParams(pj)
		params.mredParamsP[i] = MRedParams(pj)
	}

	tmp := new(big.Int)
	QiB := new(big.Int)
	QiStar := new(big.Int)
	QiBarre := new(big.Int)

	modulusbigint := NewUint(1)
	for _, qi := range Q {
		modulusbigint.Mul(modulusbigint, NewUint(qi))
	}

	params.qibMont = make([]uint64, len(Q))
	params.qispjMont = make([][]uint64, len(P))

	for i := range P {
		params.qispjMont[i] = make([]uint64, len(Q))
	}

	for i, qi := range Q {

		QiB.SetUint64(qi)
		QiStar.Quo(modulusbigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)

		// (Q/Qi)^-1) * r (mod Qi) (in Montgomery form)
		params.qibMont[i] = MForm(QiBarre.Uint64(), qi, params.bredParamsQ[i])

		for j, pj := range P {
			// (Q/qi * r) (mod Pj) (in Montgomery form)
			params.qispjMont[j][i] = MForm(tmp.Mod(QiStar, NewUint(pj)).Uint64(), pj, params.bredParamsP[j])
		}
	}

	params.qpjInv = make([][]uint64, len(P))
	for j, pj := range P {
		params.qpjInv[j] = make([]uint64, len(Q)+1)
		// Correction Term (v*Q) mod each Pj
		v := pj - tmp.Mod(modulusbigint, NewUint(pj)).Uint64()
		params.qpjInv[j][0] = 0
		for i := 1; i < len(Q)+1; i++ {
			params.qpjInv[j][i] = CRed(params.qpjInv[j][i-1]+v, pj)
		}
	}

	return
}

// ModUpSplitQP extends the RNS basis of a polynomial from Q to QP.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel}
// Extends its basis from {Q0,Q1....Qlevel} to {Q0,Q1....Qlevel,P0,P1...Pj}
func (basisextender *FastBasisExtender) ModUpSplitQP(level uint64, p1, p2 *Poly) {
	modUpExact(p1.Coeffs[:level+1], p2.Coeffs[:uint64(len(basisextender.paramsQP.P))], basisextender.paramsQP)
}

// ModUpSplitPQ extends the RNS basis of a polynomial from P to PQ.
// Given a polynomial with coefficients in basis {P0,P1....Plevel}
// Extends its basis from {P0,P1....Plevel} to {Q0,Q1...Qj}
func (basisextender *FastBasisExtender) ModUpSplitPQ(level uint64, p1, p2 *Poly) {
	modUpExact(p1.Coeffs[:level+1], p2.Coeffs[:uint64(len(basisextender.paramsPQ.P))], basisextender.paramsPQ)
}

// ModDownNTTPQ reduces the basis RNS of a polynomial in the NTT domain
//  from QP to Q and divides its coefficients by P.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel,P0,P1...Pj}
// Reduces its basis from {Q0,Q1....Qlevel,P0,P1...Pj} to {Q0,Q1....Qlevel}
// and does a runded integer division of the result by P.
// Inputs must be in the NTT domain.
func (basisextender *FastBasisExtender) ModDownNTTPQ(level uint64, p1, p2 *Poly) {

	contextQ := basisextender.contextQ
	contextP := basisextender.contextP
	modDownParams := basisextender.modDownParamsPQ
	polypool := basisextender.polypoolQ
	nQi := len(contextQ.Modulus)
	nPj := len(contextP.Modulus)

	// First we get the P basis part of p1 out of the NTT domain
	for j := 0; j < nPj; j++ {
		InvNTT(p1.Coeffs[nQi+j], p1.Coeffs[nQi+j], contextP.N, contextP.GetNttPsiInv()[j], contextP.GetNttNInv()[j], contextP.Modulus[j], contextP.GetMredParams()[j])
	}

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	modUpExact(p1.Coeffs[nQi:nQi+nPj], polypool.Coeffs[:level+1], basisextender.paramsPQ)

	// Finaly, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := contextQ.Modulus[i]
		p1tmp := p1.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := modDownParams[i]
		mredParams := contextQ.mredParams[i]
		bredParams := contextQ.bredParams[i]

		// First we switch back the relevant polypool CRT array back to the NTT domain
		NTT(p3tmp, p3tmp, contextQ.N, contextQ.GetNttPsi()[i], qi, mredParams, bredParams)

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < contextQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(x[0]+(qi-y[0]), params, qi, mredParams)
			z[1] = MRed(x[1]+(qi-y[1]), params, qi, mredParams)
			z[2] = MRed(x[2]+(qi-y[2]), params, qi, mredParams)
			z[3] = MRed(x[3]+(qi-y[3]), params, qi, mredParams)
			z[4] = MRed(x[4]+(qi-y[4]), params, qi, mredParams)
			z[5] = MRed(x[5]+(qi-y[5]), params, qi, mredParams)
			z[6] = MRed(x[6]+(qi-y[6]), params, qi, mredParams)
			z[7] = MRed(x[7]+(qi-y[7]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownSplitedNTTPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi} and {P0,P1...Pj}
// Reduces its basis from {Q0,Q1....Qi} and {P0,P1...Pj} to {Q0,Q1....Qi}
// and does a runded integer division of the result by P.
// Inputs must be in the NTT domain.
func (basisextender *FastBasisExtender) ModDownSplitedNTTPQ(level uint64, p1Q, p1P, p2 *Poly) {

	contextQ := basisextender.contextQ
	contextP := basisextender.contextP
	modDownParams := basisextender.modDownParamsPQ
	polypool := basisextender.polypoolQ

	// First we get the P basis part of p1 out of the NTT domain
	contextP.InvNTT(p1P, p1P)

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	modUpExact(p1P.Coeffs, polypool.Coeffs[:level+1], basisextender.paramsPQ)

	// Finaly, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := contextQ.Modulus[i]
		p1tmp := p1Q.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := modDownParams[i]
		mredParams := contextQ.mredParams[i]
		bredParams := contextQ.bredParams[i]

		// First we switch back the relevant polypool CRT array back to the NTT domain
		NTT(p3tmp, p3tmp, contextQ.N, contextQ.GetNttPsi()[i], contextQ.Modulus[i], mredParams, bredParams)

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < contextQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(x[0]+(qi-y[0]), params, qi, mredParams)
			z[1] = MRed(x[1]+(qi-y[1]), params, qi, mredParams)
			z[2] = MRed(x[2]+(qi-y[2]), params, qi, mredParams)
			z[3] = MRed(x[3]+(qi-y[3]), params, qi, mredParams)
			z[4] = MRed(x[4]+(qi-y[4]), params, qi, mredParams)
			z[5] = MRed(x[5]+(qi-y[5]), params, qi, mredParams)
			z[6] = MRed(x[6]+(qi-y[6]), params, qi, mredParams)
			z[7] = MRed(x[7]+(qi-y[7]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel,P0,P1...Pj}
// Reduces its basis from {Q0,Q1....Qlevel,P0,P1...Pj} to {Q0,Q1....Qlevel}
// and does a runded integer division of the result by P.
func (basisextender *FastBasisExtender) ModDownPQ(level uint64, p1, p2 *Poly) {

	context := basisextender.contextQ
	modDownParams := basisextender.modDownParamsPQ
	polypool := basisextender.polypoolQ
	nPi := uint64(len(basisextender.paramsQP.P))

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	modUpExact(p1.Coeffs[level+1:level+1+nPi], polypool.Coeffs[:level+1], basisextender.paramsPQ)

	// Finaly, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := context.Modulus[i]
		p1tmp := p1.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := modDownParams[i]
		mredParams := context.mredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < context.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(x[0]+(qi-y[0]), params, qi, mredParams)
			z[1] = MRed(x[1]+(qi-y[1]), params, qi, mredParams)
			z[2] = MRed(x[2]+(qi-y[2]), params, qi, mredParams)
			z[3] = MRed(x[3]+(qi-y[3]), params, qi, mredParams)
			z[4] = MRed(x[4]+(qi-y[4]), params, qi, mredParams)
			z[5] = MRed(x[5]+(qi-y[5]), params, qi, mredParams)
			z[6] = MRed(x[6]+(qi-y[6]), params, qi, mredParams)
			z[7] = MRed(x[7]+(qi-y[7]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownSplitedPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel} and {P0,P1...Pj}
// Reduces its basis from {Q0,Q1....Qlevel} and {P0,P1...Pj} to {Q0,Q1....Qlevel}
// and does a runded integer division of the result by P.
func (basisextender *FastBasisExtender) ModDownSplitedPQ(level uint64, p1Q, p1P, p2 *Poly) {

	contextQ := basisextender.contextQ
	modDownParams := basisextender.modDownParamsPQ
	polypool := basisextender.polypoolQ

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	modUpExact(p1P.Coeffs, polypool.Coeffs[:level+1], basisextender.paramsPQ)

	// Finaly, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := contextQ.Modulus[i]
		p1tmp := p1Q.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := modDownParams[i]
		mredParams := contextQ.mredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < contextQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(x[0]+(qi-y[0]), params, qi, mredParams)
			z[1] = MRed(x[1]+(qi-y[1]), params, qi, mredParams)
			z[2] = MRed(x[2]+(qi-y[2]), params, qi, mredParams)
			z[3] = MRed(x[3]+(qi-y[3]), params, qi, mredParams)
			z[4] = MRed(x[4]+(qi-y[4]), params, qi, mredParams)
			z[5] = MRed(x[5]+(qi-y[5]), params, qi, mredParams)
			z[6] = MRed(x[6]+(qi-y[6]), params, qi, mredParams)
			z[7] = MRed(x[7]+(qi-y[7]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownSplitedQP reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....QlevelQ} and {P0,P1...PlevelP}
// Reduces its basis from {Q0,Q1....QlevelQ} and {P0,P1...PlevelP} to {P0,P1...PlevelP}
// and does a floored integer division of the result by Q.
func (basisextender *FastBasisExtender) ModDownSplitedQP(levelQ, levelP uint64, p1Q, p1P, p2 *Poly) {

	contextP := basisextender.contextP
	//contextQ := basisextender.contextQ
	modDownParams := basisextender.modDownParamsQP
	polypool := basisextender.polypoolP

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	basisextender.ModUpSplitQP(levelQ, p1Q, polypool)

	// Finaly, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < levelP+1; i++ {

		qi := contextP.Modulus[i]
		p1tmp := p1P.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := modDownParams[i]
		mredParams := contextP.mredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < contextP.N; j++ {
			p2tmp[j] = MRed(p1tmp[j]+(qi-p3tmp[j]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

func modUpExact(p1, p2 [][]uint64, params *modupParams) {

	var v0, v1, v2, v3, v4, v5, v6, v7 uint64
	var vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7 float64
	var xpj0, xpj1, xpj2, xpj3, xpj4, xpj5, xpj6, xpj7 uint64

	y0 := make([]uint64, len(p1), len(p1))
	y1 := make([]uint64, len(p1), len(p1))
	y2 := make([]uint64, len(p1), len(p1))
	y3 := make([]uint64, len(p1), len(p1))
	y4 := make([]uint64, len(p1), len(p1))
	y5 := make([]uint64, len(p1), len(p1))
	y6 := make([]uint64, len(p1), len(p1))
	y7 := make([]uint64, len(p1), len(p1))

	var qibMont, qi, pj, mredParams uint64
	var qif float64

	//We loop over each coefficient and apply the basis extension
	for x := uint64(0); x < uint64(len(p1[0])); x = x + 8 {

		vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7 = 0, 0, 0, 0, 0, 0, 0, 0

		for i := 0; i < len(p1); i++ {

			qibMont = params.qibMont[i]
			qi = params.Q[i]
			mredParams = params.mredParamsQ[i]
			qif = float64(qi)

			y0[i] = MRed(p1[i][x], qibMont, qi, mredParams)
			y1[i] = MRed(p1[i][x+1], qibMont, qi, mredParams)
			y2[i] = MRed(p1[i][x+2], qibMont, qi, mredParams)
			y3[i] = MRed(p1[i][x+3], qibMont, qi, mredParams)
			y4[i] = MRed(p1[i][x+4], qibMont, qi, mredParams)
			y5[i] = MRed(p1[i][x+5], qibMont, qi, mredParams)
			y6[i] = MRed(p1[i][x+6], qibMont, qi, mredParams)
			y7[i] = MRed(p1[i][x+7], qibMont, qi, mredParams)

			// Computation of the correction term v * Q%pi
			vi0 += float64(y0[i]) / qif
			vi1 += float64(y1[i]) / qif
			vi2 += float64(y2[i]) / qif
			vi3 += float64(y3[i]) / qif
			vi4 += float64(y4[i]) / qif
			vi5 += float64(y5[i]) / qif
			vi6 += float64(y6[i]) / qif
			vi7 += float64(y7[i]) / qif
		}

		// Index of the correction term
		v0, v1, v2, v3, v4, v5, v6, v7 = uint64(vi0), uint64(vi1), uint64(vi2), uint64(vi3), uint64(vi4), uint64(vi5), uint64(vi6), uint64(vi7)

		for j := 0; j < len(p2); j++ {

			xpj0, xpj1, xpj2, xpj3, xpj4, xpj5, xpj6, xpj7 = 0, 0, 0, 0, 0, 0, 0, 0

			pj = params.P[j]
			mredParams = params.mredParamsP[j]
			bredParams := params.bredParamsP[j]
			qpjInv := params.qpjInv[j]
			qispjMont := params.qispjMont[j]
			res := p2[j]

			for i := 0; i < len(p1); i++ {

				xpj0 += MRed(y0[i], qispjMont[i], pj, mredParams)
				xpj1 += MRed(y1[i], qispjMont[i], pj, mredParams)
				xpj2 += MRed(y2[i], qispjMont[i], pj, mredParams)
				xpj3 += MRed(y3[i], qispjMont[i], pj, mredParams)
				xpj4 += MRed(y4[i], qispjMont[i], pj, mredParams)
				xpj5 += MRed(y5[i], qispjMont[i], pj, mredParams)
				xpj6 += MRed(y6[i], qispjMont[i], pj, mredParams)
				xpj7 += MRed(y7[i], qispjMont[i], pj, mredParams)

				if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
					xpj0 = BRedAdd(xpj0, pj, bredParams)
					xpj1 = BRedAdd(xpj1, pj, bredParams)
					xpj2 = BRedAdd(xpj2, pj, bredParams)
					xpj3 = BRedAdd(xpj3, pj, bredParams)
					xpj4 = BRedAdd(xpj4, pj, bredParams)
					xpj5 = BRedAdd(xpj5, pj, bredParams)
					xpj6 = BRedAdd(xpj6, pj, bredParams)
					xpj7 = BRedAdd(xpj7, pj, bredParams)
				}
			}

			res[x+0] = BRedAdd(xpj0+qpjInv[v0], pj, bredParams)
			res[x+1] = BRedAdd(xpj1+qpjInv[v1], pj, bredParams)
			res[x+2] = BRedAdd(xpj2+qpjInv[v2], pj, bredParams)
			res[x+3] = BRedAdd(xpj3+qpjInv[v3], pj, bredParams)
			res[x+4] = BRedAdd(xpj4+qpjInv[v4], pj, bredParams)
			res[x+5] = BRedAdd(xpj5+qpjInv[v5], pj, bredParams)
			res[x+6] = BRedAdd(xpj6+qpjInv[v6], pj, bredParams)
			res[x+7] = BRedAdd(xpj7+qpjInv[v7], pj, bredParams)

		}
	}
}

// Decomposer is a structure storing the parameters of the arbitrary decomposer.
// This decomposer takes a p(x)_Q (in basis Q) and returns p(x) mod qi in basis QP. Where
// qi = prod(Q_i) for 0<=i<=L where L is the number of factors in P.
type Decomposer struct {
	nQprimes    uint64
	nPprimes    uint64
	alpha       uint64
	beta        uint64
	xalpha      []uint64
	modUpParams [][]*modupParams
	QInt        *big.Int
	PInt        *big.Int
}

// Xalpha returns a slice containing all the values of #Qi/#Pi.
func (decomposer *Decomposer) Xalpha() (xalpha []uint64) {
	return decomposer.xalpha
}

// NewDecomposer creates a new Decomposer.
func NewDecomposer(Q, P []uint64) (decomposer *Decomposer) {
	decomposer = new(Decomposer)

	decomposer.nQprimes = uint64(len(Q))
	decomposer.nPprimes = uint64(len(P))

	decomposer.QInt = NewUint(1)
	for i := range Q {
		decomposer.QInt.Mul(decomposer.QInt, NewUint(Q[i]))
	}

	decomposer.PInt = NewUint(1)
	for i := range P {
		decomposer.PInt.Mul(decomposer.PInt, NewUint(P[i]))
	}

	decomposer.alpha = uint64(len(P))
	decomposer.beta = uint64(math.Ceil(float64(len(Q)) / float64(decomposer.alpha)))

	decomposer.xalpha = make([]uint64, decomposer.beta)
	for i := range decomposer.xalpha {
		decomposer.xalpha[i] = decomposer.alpha
	}

	if uint64(len(Q))%decomposer.alpha != 0 {
		decomposer.xalpha[decomposer.beta-1] = uint64(len(Q)) % decomposer.alpha
	}

	decomposer.modUpParams = make([][]*modupParams, decomposer.beta)

	// Creates a basis extension for each possible combination of [Qi,Pj] according to xalpha
	for i := uint64(0); i < decomposer.beta; i++ {

		decomposer.modUpParams[i] = make([]*modupParams, decomposer.xalpha[i]-1)

		for j := uint64(0); j < decomposer.xalpha[i]-1; j++ {

			Qi := make([]uint64, j+2)
			Pi := make([]uint64, len(Q)+len(P))

			for k := uint64(0); k < j+2; k++ {
				Qi[k] = Q[i*decomposer.alpha+k]
			}

			for k := 0; k < len(Q); k++ {
				Pi[k] = Q[k]
			}

			for k := len(Q); k < len(Q)+len(P); k++ {
				Pi[k] = P[k-len(Q)]
			}

			decomposer.modUpParams[i][j] = basisextenderparameters(Qi, Pi)
		}
	}

	return
}

// Decompose decomposes takes a polynomial p(x) in basis Q, reduces it modulo qi, and returns
// the result in basis QP.
func (decomposer *Decomposer) Decompose(level, crtDecompLevel uint64, p0, p1 *Poly) {

	alphai := decomposer.xalpha[crtDecompLevel]

	p0idxst := crtDecompLevel * decomposer.alpha
	p0idxed := p0idxst + alphai

	// First we check if the vector can simply by coping and rearanging elements (the case where no reconstruction is needed)
	if (p0idxed > level+1 && (level+1)%decomposer.nPprimes == 1) || alphai == 1 {

		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x = x + 8 {

			tmp := p0.Coeffs[p0idxst]

			for j := uint64(0); j < level+decomposer.nPprimes+1; j++ {

				p1.Coeffs[j][x+0] = tmp[x+0]
				p1.Coeffs[j][x+1] = tmp[x+1]
				p1.Coeffs[j][x+2] = tmp[x+2]
				p1.Coeffs[j][x+3] = tmp[x+3]
				p1.Coeffs[j][x+4] = tmp[x+4]
				p1.Coeffs[j][x+5] = tmp[x+5]
				p1.Coeffs[j][x+6] = tmp[x+6]
				p1.Coeffs[j][x+7] = tmp[x+7]
			}
		}

		// Else we apply a fast exact base conversion for the reconstruction
	} else {

		var index uint64
		if level >= alphai+crtDecompLevel*decomposer.alpha {
			//fmt.Println("A")
			index = decomposer.xalpha[crtDecompLevel] - 2
		} else {
			//fmt.Println("B")
			index = (level - 1) % decomposer.alpha
		}

		params := decomposer.modUpParams[crtDecompLevel][index]

		v := make([]uint64, 8, 8)
		vi := make([]float64, 8, 8)
		xpj := make([]uint64, 8, 8)

		y0 := make([]uint64, index+2, index+2)
		y1 := make([]uint64, index+2, index+2)
		y2 := make([]uint64, index+2, index+2)
		y3 := make([]uint64, index+2, index+2)
		y4 := make([]uint64, index+2, index+2)
		y5 := make([]uint64, index+2, index+2)
		y6 := make([]uint64, index+2, index+2)
		y7 := make([]uint64, index+2, index+2)

		var qibMont, qi, pj, mredParams uint64
		var qif float64

		//We loop over each coefficient and apply the basis extension
		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x = x + 8 {

			vi[0], vi[1], vi[2], vi[3], vi[4], vi[5], vi[6], vi[7] = 0, 0, 0, 0, 0, 0, 0, 0

			// Coefficients to be decomposed
			for i, j := uint64(0), p0idxst; i < index+2; i, j = i+1, j+1 {

				qibMont = params.qibMont[i]
				qi = params.Q[i]
				mredParams = params.mredParamsQ[i]
				qif = float64(qi)

				px := p0.Coeffs[j]
				py := p1.Coeffs[j]

				// For the coefficients to be decomposed, we can simplly copy them
				py[x+0] = px[x+0]
				py[x+1] = px[x+1]
				py[x+2] = px[x+2]
				py[x+3] = px[x+3]
				py[x+4] = px[x+4]
				py[x+5] = px[x+5]
				py[x+6] = px[x+6]
				py[x+7] = px[x+7]

				y0[i] = MRed(px[x+0], qibMont, qi, mredParams)
				y1[i] = MRed(px[x+1], qibMont, qi, mredParams)
				y2[i] = MRed(px[x+2], qibMont, qi, mredParams)
				y3[i] = MRed(px[x+3], qibMont, qi, mredParams)
				y4[i] = MRed(px[x+4], qibMont, qi, mredParams)
				y5[i] = MRed(px[x+5], qibMont, qi, mredParams)
				y6[i] = MRed(px[x+6], qibMont, qi, mredParams)
				y7[i] = MRed(px[x+7], qibMont, qi, mredParams)

				// Computation of the correction term v * Q%pi
				vi[0] += float64(y0[i]) / qif
				vi[1] += float64(y1[i]) / qif
				vi[2] += float64(y2[i]) / qif
				vi[3] += float64(y3[i]) / qif
				vi[4] += float64(y4[i]) / qif
				vi[5] += float64(y5[i]) / qif
				vi[6] += float64(y6[i]) / qif
				vi[7] += float64(y7[i]) / qif
			}

			// Index of the correction term
			v[0] = uint64(vi[0])
			v[1] = uint64(vi[1])
			v[2] = uint64(vi[2])
			v[3] = uint64(vi[3])
			v[4] = uint64(vi[4])
			v[5] = uint64(vi[5])
			v[6] = uint64(vi[6])
			v[7] = uint64(vi[7])

			// Coefficients of index smaller than the ones to be decomposer
			for j := uint64(0); j < p0idxst; j++ {

				xpj[0], xpj[1], xpj[2], xpj[3], xpj[4], xpj[5], xpj[6], xpj[7] = 0, 0, 0, 0, 0, 0, 0, 0

				pj = params.P[j]
				mredParams = params.mredParamsP[j]
				bredParams := params.bredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]
				res := p1.Coeffs[j]

				for i := uint64(0); i < index+2; i++ {

					xpj[0] += MRed(y0[i], qispjMont[i], pj, mredParams)
					xpj[1] += MRed(y1[i], qispjMont[i], pj, mredParams)
					xpj[2] += MRed(y2[i], qispjMont[i], pj, mredParams)
					xpj[3] += MRed(y3[i], qispjMont[i], pj, mredParams)
					xpj[4] += MRed(y4[i], qispjMont[i], pj, mredParams)
					xpj[5] += MRed(y5[i], qispjMont[i], pj, mredParams)
					xpj[6] += MRed(y6[i], qispjMont[i], pj, mredParams)
					xpj[7] += MRed(y7[i], qispjMont[i], pj, mredParams)

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj[0] = BRedAdd(xpj[0], pj, bredParams)
						xpj[1] = BRedAdd(xpj[1], pj, bredParams)
						xpj[2] = BRedAdd(xpj[2], pj, bredParams)
						xpj[3] = BRedAdd(xpj[3], pj, bredParams)
						xpj[4] = BRedAdd(xpj[4], pj, bredParams)
						xpj[5] = BRedAdd(xpj[5], pj, bredParams)
						xpj[6] = BRedAdd(xpj[6], pj, bredParams)
						xpj[7] = BRedAdd(xpj[7], pj, bredParams)
					}
				}

				res[x+0] = BRedAdd(xpj[0]+qpjInv[v[0]], pj, bredParams)
				res[x+1] = BRedAdd(xpj[1]+qpjInv[v[1]], pj, bredParams)
				res[x+2] = BRedAdd(xpj[2]+qpjInv[v[2]], pj, bredParams)
				res[x+3] = BRedAdd(xpj[3]+qpjInv[v[3]], pj, bredParams)
				res[x+4] = BRedAdd(xpj[4]+qpjInv[v[4]], pj, bredParams)
				res[x+5] = BRedAdd(xpj[5]+qpjInv[v[5]], pj, bredParams)
				res[x+6] = BRedAdd(xpj[6]+qpjInv[v[6]], pj, bredParams)
				res[x+7] = BRedAdd(xpj[7]+qpjInv[v[7]], pj, bredParams)

			}

			// Coefficients of index greater than the ones to be decomposer
			for j := decomposer.alpha * crtDecompLevel; j < level+1; j = j + 1 {

				xpj[0], xpj[1], xpj[2], xpj[3], xpj[4], xpj[5], xpj[6], xpj[7] = 0, 0, 0, 0, 0, 0, 0, 0

				pj = params.P[j]
				mredParams = params.mredParamsP[j]
				bredParams := params.bredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]
				res := p1.Coeffs[j]

				for i := uint64(0); i < index+2; i++ {

					xpj[0] += MRed(y0[i], qispjMont[i], pj, mredParams)
					xpj[1] += MRed(y1[i], qispjMont[i], pj, mredParams)
					xpj[2] += MRed(y2[i], qispjMont[i], pj, mredParams)
					xpj[3] += MRed(y3[i], qispjMont[i], pj, mredParams)
					xpj[4] += MRed(y4[i], qispjMont[i], pj, mredParams)
					xpj[5] += MRed(y5[i], qispjMont[i], pj, mredParams)
					xpj[6] += MRed(y6[i], qispjMont[i], pj, mredParams)
					xpj[7] += MRed(y7[i], qispjMont[i], pj, mredParams)

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj[0] = BRedAdd(xpj[0], pj, bredParams)
						xpj[1] = BRedAdd(xpj[1], pj, bredParams)
						xpj[2] = BRedAdd(xpj[2], pj, bredParams)
						xpj[3] = BRedAdd(xpj[3], pj, bredParams)
						xpj[4] = BRedAdd(xpj[4], pj, bredParams)
						xpj[5] = BRedAdd(xpj[5], pj, bredParams)
						xpj[6] = BRedAdd(xpj[6], pj, bredParams)
						xpj[7] = BRedAdd(xpj[7], pj, bredParams)
					}
				}

				res[x+0] = BRedAdd(xpj[0]+qpjInv[v[0]], pj, bredParams)
				res[x+1] = BRedAdd(xpj[1]+qpjInv[v[1]], pj, bredParams)
				res[x+2] = BRedAdd(xpj[2]+qpjInv[v[2]], pj, bredParams)
				res[x+3] = BRedAdd(xpj[3]+qpjInv[v[3]], pj, bredParams)
				res[x+4] = BRedAdd(xpj[4]+qpjInv[v[4]], pj, bredParams)
				res[x+5] = BRedAdd(xpj[5]+qpjInv[v[5]], pj, bredParams)
				res[x+6] = BRedAdd(xpj[6]+qpjInv[v[6]], pj, bredParams)
				res[x+7] = BRedAdd(xpj[7]+qpjInv[v[7]], pj, bredParams)

			}

			// Coefficients of the special primes
			for j, u := level+1, decomposer.nQprimes; j < level+1+decomposer.nPprimes; j, u = u+1, j+1 {

				xpj[0], xpj[1], xpj[2], xpj[3], xpj[4], xpj[5], xpj[6], xpj[7] = 0, 0, 0, 0, 0, 0, 0, 0

				pj = params.P[j]
				mredParams = params.mredParamsP[j]
				bredParams := params.bredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]
				res := p1.Coeffs[u]

				for i := uint64(0); i < index+2; i++ {

					xpj[0] += MRed(y0[i], qispjMont[i], pj, mredParams)
					xpj[1] += MRed(y1[i], qispjMont[i], pj, mredParams)
					xpj[2] += MRed(y2[i], qispjMont[i], pj, mredParams)
					xpj[3] += MRed(y3[i], qispjMont[i], pj, mredParams)
					xpj[4] += MRed(y4[i], qispjMont[i], pj, mredParams)
					xpj[5] += MRed(y5[i], qispjMont[i], pj, mredParams)
					xpj[6] += MRed(y6[i], qispjMont[i], pj, mredParams)
					xpj[7] += MRed(y7[i], qispjMont[i], pj, mredParams)

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj[0] = BRedAdd(xpj[0], pj, bredParams)
						xpj[1] = BRedAdd(xpj[1], pj, bredParams)
						xpj[2] = BRedAdd(xpj[2], pj, bredParams)
						xpj[3] = BRedAdd(xpj[3], pj, bredParams)
						xpj[4] = BRedAdd(xpj[4], pj, bredParams)
						xpj[5] = BRedAdd(xpj[5], pj, bredParams)
						xpj[6] = BRedAdd(xpj[6], pj, bredParams)
						xpj[7] = BRedAdd(xpj[7], pj, bredParams)
					}
				}

				res[x+0] = BRedAdd(xpj[0]+qpjInv[v[0]], pj, bredParams)
				res[x+1] = BRedAdd(xpj[1]+qpjInv[v[1]], pj, bredParams)
				res[x+2] = BRedAdd(xpj[2]+qpjInv[v[2]], pj, bredParams)
				res[x+3] = BRedAdd(xpj[3]+qpjInv[v[3]], pj, bredParams)
				res[x+4] = BRedAdd(xpj[4]+qpjInv[v[4]], pj, bredParams)
				res[x+5] = BRedAdd(xpj[5]+qpjInv[v[5]], pj, bredParams)
				res[x+6] = BRedAdd(xpj[6]+qpjInv[v[6]], pj, bredParams)
				res[x+7] = BRedAdd(xpj[7]+qpjInv[v[7]], pj, bredParams)
			}
		}
	}
}

// DecomposeAndSplit decomposes takes a polynomial p(x) in basis Q, reduces it modulo qi, and returns
// the result in basis QP separately.
func (decomposer *Decomposer) DecomposeAndSplit(level, crtDecompLevel uint64, p0, p1Q, p1P *Poly) {

	alphai := decomposer.xalpha[crtDecompLevel]

	p0idxst := crtDecompLevel * decomposer.alpha
	p0idxed := p0idxst + alphai

	// First we check if the vector can simply by coping and rearanging elements (the case where no reconstruction is needed)
	if (p0idxed > level+1 && (level+1)%decomposer.nPprimes == 1) || alphai == 1 {

		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x = x + 8 {

			tmp := p0.Coeffs[p0idxst]

			for j := uint64(0); j < level+1; j++ {

				p1Q.Coeffs[j][x+0] = tmp[x+0]
				p1Q.Coeffs[j][x+1] = tmp[x+1]
				p1Q.Coeffs[j][x+2] = tmp[x+2]
				p1Q.Coeffs[j][x+3] = tmp[x+3]
				p1Q.Coeffs[j][x+4] = tmp[x+4]
				p1Q.Coeffs[j][x+5] = tmp[x+5]
				p1Q.Coeffs[j][x+6] = tmp[x+6]
				p1Q.Coeffs[j][x+7] = tmp[x+7]
			}

			for j := uint64(0); j < decomposer.nPprimes; j++ {

				p1P.Coeffs[j][x+0] = tmp[x+0]
				p1P.Coeffs[j][x+1] = tmp[x+1]
				p1P.Coeffs[j][x+2] = tmp[x+2]
				p1P.Coeffs[j][x+3] = tmp[x+3]
				p1P.Coeffs[j][x+4] = tmp[x+4]
				p1P.Coeffs[j][x+5] = tmp[x+5]
				p1P.Coeffs[j][x+6] = tmp[x+6]
				p1P.Coeffs[j][x+7] = tmp[x+7]
			}
		}

		// Else we apply a fast exact base conversion for the reconstruction
	} else {

		var index uint64
		if level >= alphai+crtDecompLevel*decomposer.alpha {
			//fmt.Println("A")
			index = decomposer.xalpha[crtDecompLevel] - 2
		} else {
			//fmt.Println("B")
			index = (level - 1) % decomposer.alpha
		}

		params := decomposer.modUpParams[crtDecompLevel][index]

		v := make([]uint64, 8, 8)
		vi := make([]float64, 8, 8)
		xpj := make([]uint64, 8, 8)

		y0 := make([]uint64, index+2, index+2)
		y1 := make([]uint64, index+2, index+2)
		y2 := make([]uint64, index+2, index+2)
		y3 := make([]uint64, index+2, index+2)
		y4 := make([]uint64, index+2, index+2)
		y5 := make([]uint64, index+2, index+2)
		y6 := make([]uint64, index+2, index+2)
		y7 := make([]uint64, index+2, index+2)

		var qibMont, qi, pj, mredParams uint64
		var qif float64

		//We loop over each coefficient and apply the basis extension
		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x = x + 8 {

			vi[0], vi[1], vi[2], vi[3], vi[4], vi[5], vi[6], vi[7] = 0, 0, 0, 0, 0, 0, 0, 0

			// Coefficients to be decomposed
			for i, j := uint64(0), p0idxst; i < index+2; i, j = i+1, j+1 {

				qibMont = params.qibMont[i]
				qi = params.Q[i]
				mredParams = params.mredParamsQ[i]
				qif = float64(qi)

				px := p0.Coeffs[j]
				py := p1Q.Coeffs[j]

				// For the coefficients to be decomposed, we can simplly copy them
				py[x+0] = px[x+0]
				py[x+1] = px[x+1]
				py[x+2] = px[x+2]
				py[x+3] = px[x+3]
				py[x+4] = px[x+4]
				py[x+5] = px[x+5]
				py[x+6] = px[x+6]
				py[x+7] = px[x+7]

				y0[i] = MRed(px[x+0], qibMont, qi, mredParams)
				y1[i] = MRed(px[x+1], qibMont, qi, mredParams)
				y2[i] = MRed(px[x+2], qibMont, qi, mredParams)
				y3[i] = MRed(px[x+3], qibMont, qi, mredParams)
				y4[i] = MRed(px[x+4], qibMont, qi, mredParams)
				y5[i] = MRed(px[x+5], qibMont, qi, mredParams)
				y6[i] = MRed(px[x+6], qibMont, qi, mredParams)
				y7[i] = MRed(px[x+7], qibMont, qi, mredParams)

				// Computation of the correction term v * Q%pi
				vi[0] += float64(y0[i]) / qif
				vi[1] += float64(y1[i]) / qif
				vi[2] += float64(y2[i]) / qif
				vi[3] += float64(y3[i]) / qif
				vi[4] += float64(y4[i]) / qif
				vi[5] += float64(y5[i]) / qif
				vi[6] += float64(y6[i]) / qif
				vi[7] += float64(y7[i]) / qif
			}

			// Index of the correction term
			v[0] = uint64(vi[0])
			v[1] = uint64(vi[1])
			v[2] = uint64(vi[2])
			v[3] = uint64(vi[3])
			v[4] = uint64(vi[4])
			v[5] = uint64(vi[5])
			v[6] = uint64(vi[6])
			v[7] = uint64(vi[7])

			// Coefficients of index smaller than the ones to be decomposer
			for j := uint64(0); j < p0idxst; j++ {

				xpj[0], xpj[1], xpj[2], xpj[3], xpj[4], xpj[5], xpj[6], xpj[7] = 0, 0, 0, 0, 0, 0, 0, 0

				pj = params.P[j]
				mredParams = params.mredParamsP[j]
				bredParams := params.bredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]
				res := p1Q.Coeffs[j]

				for i := uint64(0); i < index+2; i++ {

					xpj[0] += MRed(y0[i], qispjMont[i], pj, mredParams)
					xpj[1] += MRed(y1[i], qispjMont[i], pj, mredParams)
					xpj[2] += MRed(y2[i], qispjMont[i], pj, mredParams)
					xpj[3] += MRed(y3[i], qispjMont[i], pj, mredParams)
					xpj[4] += MRed(y4[i], qispjMont[i], pj, mredParams)
					xpj[5] += MRed(y5[i], qispjMont[i], pj, mredParams)
					xpj[6] += MRed(y6[i], qispjMont[i], pj, mredParams)
					xpj[7] += MRed(y7[i], qispjMont[i], pj, mredParams)

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj[0] = BRedAdd(xpj[0], pj, bredParams)
						xpj[1] = BRedAdd(xpj[1], pj, bredParams)
						xpj[2] = BRedAdd(xpj[2], pj, bredParams)
						xpj[3] = BRedAdd(xpj[3], pj, bredParams)
						xpj[4] = BRedAdd(xpj[4], pj, bredParams)
						xpj[5] = BRedAdd(xpj[5], pj, bredParams)
						xpj[6] = BRedAdd(xpj[6], pj, bredParams)
						xpj[7] = BRedAdd(xpj[7], pj, bredParams)
					}
				}

				res[x+0] = BRedAdd(xpj[0]+qpjInv[v[0]], pj, bredParams)
				res[x+1] = BRedAdd(xpj[1]+qpjInv[v[1]], pj, bredParams)
				res[x+2] = BRedAdd(xpj[2]+qpjInv[v[2]], pj, bredParams)
				res[x+3] = BRedAdd(xpj[3]+qpjInv[v[3]], pj, bredParams)
				res[x+4] = BRedAdd(xpj[4]+qpjInv[v[4]], pj, bredParams)
				res[x+5] = BRedAdd(xpj[5]+qpjInv[v[5]], pj, bredParams)
				res[x+6] = BRedAdd(xpj[6]+qpjInv[v[6]], pj, bredParams)
				res[x+7] = BRedAdd(xpj[7]+qpjInv[v[7]], pj, bredParams)

			}

			// Coefficients of index greater than the ones to be decomposer
			for j := decomposer.alpha * crtDecompLevel; j < level+1; j = j + 1 {

				xpj[0], xpj[1], xpj[2], xpj[3], xpj[4], xpj[5], xpj[6], xpj[7] = 0, 0, 0, 0, 0, 0, 0, 0

				pj = params.P[j]
				mredParams = params.mredParamsP[j]
				bredParams := params.bredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]
				res := p1Q.Coeffs[j]

				for i := uint64(0); i < index+2; i++ {

					xpj[0] += MRed(y0[i], qispjMont[i], pj, mredParams)
					xpj[1] += MRed(y1[i], qispjMont[i], pj, mredParams)
					xpj[2] += MRed(y2[i], qispjMont[i], pj, mredParams)
					xpj[3] += MRed(y3[i], qispjMont[i], pj, mredParams)
					xpj[4] += MRed(y4[i], qispjMont[i], pj, mredParams)
					xpj[5] += MRed(y5[i], qispjMont[i], pj, mredParams)
					xpj[6] += MRed(y6[i], qispjMont[i], pj, mredParams)
					xpj[7] += MRed(y7[i], qispjMont[i], pj, mredParams)

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj[0] = BRedAdd(xpj[0], pj, bredParams)
						xpj[1] = BRedAdd(xpj[1], pj, bredParams)
						xpj[2] = BRedAdd(xpj[2], pj, bredParams)
						xpj[3] = BRedAdd(xpj[3], pj, bredParams)
						xpj[4] = BRedAdd(xpj[4], pj, bredParams)
						xpj[5] = BRedAdd(xpj[5], pj, bredParams)
						xpj[6] = BRedAdd(xpj[6], pj, bredParams)
						xpj[7] = BRedAdd(xpj[7], pj, bredParams)
					}
				}

				res[x+0] = BRedAdd(xpj[0]+qpjInv[v[0]], pj, bredParams)
				res[x+1] = BRedAdd(xpj[1]+qpjInv[v[1]], pj, bredParams)
				res[x+2] = BRedAdd(xpj[2]+qpjInv[v[2]], pj, bredParams)
				res[x+3] = BRedAdd(xpj[3]+qpjInv[v[3]], pj, bredParams)
				res[x+4] = BRedAdd(xpj[4]+qpjInv[v[4]], pj, bredParams)
				res[x+5] = BRedAdd(xpj[5]+qpjInv[v[5]], pj, bredParams)
				res[x+6] = BRedAdd(xpj[6]+qpjInv[v[6]], pj, bredParams)
				res[x+7] = BRedAdd(xpj[7]+qpjInv[v[7]], pj, bredParams)

			}

			// Coefficients of the special primes
			for j, u := uint64(0), decomposer.nQprimes; j < decomposer.nPprimes; j, u = j+1, u+1 {

				xpj[0], xpj[1], xpj[2], xpj[3], xpj[4], xpj[5], xpj[6], xpj[7] = 0, 0, 0, 0, 0, 0, 0, 0

				pj = params.P[u]
				mredParams = params.mredParamsP[u]
				bredParams := params.bredParamsP[u]
				qpjInv := params.qpjInv[u]
				qispjMont := params.qispjMont[u]
				res := p1P.Coeffs[j]

				for i := uint64(0); i < index+2; i++ {

					xpj[0] += MRed(y0[i], qispjMont[i], pj, mredParams)
					xpj[1] += MRed(y1[i], qispjMont[i], pj, mredParams)
					xpj[2] += MRed(y2[i], qispjMont[i], pj, mredParams)
					xpj[3] += MRed(y3[i], qispjMont[i], pj, mredParams)
					xpj[4] += MRed(y4[i], qispjMont[i], pj, mredParams)
					xpj[5] += MRed(y5[i], qispjMont[i], pj, mredParams)
					xpj[6] += MRed(y6[i], qispjMont[i], pj, mredParams)
					xpj[7] += MRed(y7[i], qispjMont[i], pj, mredParams)

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj[0] = BRedAdd(xpj[0], pj, bredParams)
						xpj[1] = BRedAdd(xpj[1], pj, bredParams)
						xpj[2] = BRedAdd(xpj[2], pj, bredParams)
						xpj[3] = BRedAdd(xpj[3], pj, bredParams)
						xpj[4] = BRedAdd(xpj[4], pj, bredParams)
						xpj[5] = BRedAdd(xpj[5], pj, bredParams)
						xpj[6] = BRedAdd(xpj[6], pj, bredParams)
						xpj[7] = BRedAdd(xpj[7], pj, bredParams)
					}
				}

				res[x+0] = BRedAdd(xpj[0]+qpjInv[v[0]], pj, bredParams)
				res[x+1] = BRedAdd(xpj[1]+qpjInv[v[1]], pj, bredParams)
				res[x+2] = BRedAdd(xpj[2]+qpjInv[v[2]], pj, bredParams)
				res[x+3] = BRedAdd(xpj[3]+qpjInv[v[3]], pj, bredParams)
				res[x+4] = BRedAdd(xpj[4]+qpjInv[v[4]], pj, bredParams)
				res[x+5] = BRedAdd(xpj[5]+qpjInv[v[5]], pj, bredParams)
				res[x+6] = BRedAdd(xpj[6]+qpjInv[v[6]], pj, bredParams)
				res[x+7] = BRedAdd(xpj[7]+qpjInv[v[7]], pj, bredParams)

			}
		}
	}
}
