package ring

import (
	"math"
	"math/big"
	"math/bits"
	"unsafe"
)

// FastBasisExtender stores the necessary parameters for RNS basis extension.
// The used algorithm is from https://eprint.iacr.org/2018/117.pdf.
type FastBasisExtender struct {
	ringQ             *Ring
	ringP             *Ring
	paramsQtoP        []modupParams
	paramsPtoQ        []modupParams
	modDownparamsPtoQ [][]uint64
	modDownparamsQtoP [][]uint64

	polypoolQ *Poly
	polypoolP *Poly
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

func genModDownParams(ringQ, ringP *Ring) (params [][]uint64) {

	params = make([][]uint64, len(ringP.Modulus))

	bredParams := ringQ.BredParams
	mredParams := ringQ.MredParams

	for j := range ringP.Modulus {
		params[j] = make([]uint64, len(ringQ.Modulus))
		for i, qi := range ringQ.Modulus {
			params[j][i] = ModExp(ringP.Modulus[j], int(qi-2), qi)
			params[j][i] = MForm(params[j][i], qi, bredParams[i])

			if j > 0 {
				params[j][i] = MRed(params[j][i], params[j-1][i], qi, mredParams[i])
			}
		}
	}

	return
}

// NewFastBasisExtender creates a new FastBasisExtender, enabling RNS basis extension from Q to P and P to Q.
func NewFastBasisExtender(ringQ, ringP *Ring) *FastBasisExtender {

	newParams := new(FastBasisExtender)

	newParams.ringQ = ringQ
	newParams.ringP = ringP

	newParams.paramsQtoP = make([]modupParams, len(ringQ.Modulus))
	for i := range ringQ.Modulus {
		newParams.paramsQtoP[i] = basisextenderparameters(ringQ.Modulus[:i+1], ringP.Modulus)
	}

	newParams.paramsPtoQ = make([]modupParams, len(ringP.Modulus))
	for i := range ringP.Modulus {
		newParams.paramsPtoQ[i] = basisextenderparameters(ringP.Modulus[:i+1], ringQ.Modulus)
	}

	newParams.modDownparamsPtoQ = genModDownParams(ringQ, ringP)
	newParams.modDownparamsQtoP = genModDownParams(ringP, ringQ)

	newParams.polypoolQ = ringQ.NewPoly()
	newParams.polypoolP = ringP.NewPoly()

	return newParams
}

func basisextenderparameters(moduliQ, moduliP []uint64) modupParams {

	Q := make([]uint64, len(moduliQ))
	bredParamsQ := make([][]uint64, len(Q))
	mredParamsQ := make([]uint64, len(Q))
	for i, qi := range moduliQ {
		Q[i] = qi
		bredParamsQ[i] = BRedParams(qi)
		mredParamsQ[i] = MRedParams(qi)
	}

	P := make([]uint64, len(moduliP))
	bredParamsP := make([][]uint64, len(P))
	mredParamsP := make([]uint64, len(P))
	for i, pj := range moduliP {
		P[i] = pj
		bredParamsP[i] = BRedParams(pj)
		mredParamsP[i] = MRedParams(pj)
	}

	tmp := new(big.Int)
	QiB := new(big.Int)
	QiStar := new(big.Int)
	QiBarre := new(big.Int)

	modulusbigint := NewUint(1)
	for _, qi := range Q {
		modulusbigint.Mul(modulusbigint, NewUint(qi))
	}

	qibMont := make([]uint64, len(Q))
	qispjMont := make([][]uint64, len(P))

	for i := range P {
		qispjMont[i] = make([]uint64, len(Q))
	}

	for i, qi := range Q {

		QiB.SetUint64(qi)
		QiStar.Quo(modulusbigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)

		// (Q/Qi)^-1) * r (mod Qi) (in Montgomery form)
		qibMont[i] = MForm(QiBarre.Uint64(), qi, bredParamsQ[i])

		for j, pj := range P {
			// (Q/qi * r) (mod Pj) (in Montgomery form)
			qispjMont[j][i] = MForm(tmp.Mod(QiStar, NewUint(pj)).Uint64(), pj, bredParamsP[j])
		}
	}

	qpjInv := make([][]uint64, len(P))
	for j, pj := range P {
		qpjInv[j] = make([]uint64, len(Q)+1)
		// Correction Term (v*Q) mod each Pj
		v := pj - tmp.Mod(modulusbigint, NewUint(pj)).Uint64()
		qpjInv[j][0] = 0
		for i := 1; i < len(Q)+1; i++ {
			qpjInv[j][i] = CRed(qpjInv[j][i-1]+v, pj)
		}
	}

	return modupParams{Q: Q, P: P, qibMont: qibMont, qispjMont: qispjMont, qpjInv: qpjInv, bredParamsQ: bredParamsQ, mredParamsQ: mredParamsQ, bredParamsP: bredParamsP, mredParamsP: mredParamsP}
}

// ShallowCopy creates a shallow copy of this basis extender in which the read-only data-structures are
// shared with the receiver.
func (basisextender *FastBasisExtender) ShallowCopy() *FastBasisExtender {
	if basisextender == nil {
		return nil
	}
	return &FastBasisExtender{
		ringQ:             basisextender.ringQ,
		ringP:             basisextender.ringP,
		paramsQtoP:        basisextender.paramsQtoP,
		paramsPtoQ:        basisextender.paramsPtoQ,
		modDownparamsQtoP: basisextender.modDownparamsQtoP,
		modDownparamsPtoQ: basisextender.modDownparamsPtoQ,

		polypoolQ: basisextender.ringQ.NewPoly(),
		polypoolP: basisextender.ringP.NewPoly(),
	}
}

// ModUpQtoP extends the RNS basis of a polynomial from Q to QP.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel},
// it extends its basis from {Q0,Q1....Qlevel} to {Q0,Q1....Qlevel,P0,P1...Pj}
func (basisextender *FastBasisExtender) ModUpQtoP(levelQ, levelP int, polQ, polP *Poly) {
	modUpExact(polQ.Coeffs[:levelQ+1], polP.Coeffs[:levelP+1], basisextender.paramsQtoP[levelQ])
}

// ModUpPtoQ extends the RNS basis of a polynomial from P to PQ.
// Given a polynomial with coefficients in basis {P0,P1....Plevel},
// it extends its basis from {P0,P1....Plevel} to {Q0,Q1...Qj}
func (basisextender *FastBasisExtender) ModUpPtoQ(levelP, levelQ int, polP, polQ *Poly) {
	modUpExact(polP.Coeffs[:levelP+1], polQ.Coeffs[:levelQ+1], basisextender.paramsPtoQ[levelP])
}

// ModDownQPtoQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel} and {P0,P1...Pj},
// it reduces its basis from {Q0,Q1....Qlevel} and {P0,P1...Pj} to {Q0,Q1....Qlevel}
// and does a rounded integer division of the result by P.
func (basisextender *FastBasisExtender) ModDownQPtoQ(levelQ, levelP int, p1Q, p1P, p2Q *Poly) {

	ringQ := basisextender.ringQ
	modDownParams := basisextender.modDownparamsPtoQ
	polypool := basisextender.polypoolQ

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	basisextender.ModUpPtoQ(levelP, levelQ, p1P, polypool)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := 0; i < levelQ+1; i++ {

		qi := ringQ.Modulus[i]
		twoqi := qi << 1
		p1tmp := p1Q.Coeffs[i]
		p2tmp := p2Q.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := qi - modDownParams[levelP][i]
		mredParams := ringQ.MredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := 0; j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(y[0]+twoqi-x[0], params, qi, mredParams)
			z[1] = MRed(y[1]+twoqi-x[1], params, qi, mredParams)
			z[2] = MRed(y[2]+twoqi-x[2], params, qi, mredParams)
			z[3] = MRed(y[3]+twoqi-x[3], params, qi, mredParams)
			z[4] = MRed(y[4]+twoqi-x[4], params, qi, mredParams)
			z[5] = MRed(y[5]+twoqi-x[5], params, qi, mredParams)
			z[6] = MRed(y[6]+twoqi-x[6], params, qi, mredParams)
			z[7] = MRed(y[7]+twoqi-x[7], params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownQPtoQNTT reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi} and {P0,P1...Pj},
// it reduces its basis from {Q0,Q1....Qi} and {P0,P1...Pj} to {Q0,Q1....Qi}
// and does a rounded integer division of the result by P.
// Inputs must be in the NTT domain.
func (basisextender *FastBasisExtender) ModDownQPtoQNTT(levelQ, levelP int, p1Q, p1P, p2Q *Poly) {

	ringQ := basisextender.ringQ
	ringP := basisextender.ringP
	modDownParams := basisextender.modDownparamsPtoQ
	polypool := basisextender.polypoolQ

	// First we get the P basis part of p1 out of the NTT domain
	ringP.InvNTTLazyLvl(levelP, p1P, p1P)

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	basisextender.ModUpPtoQ(levelP, levelQ, p1P, polypool)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := 0; i < levelQ+1; i++ {

		qi := ringQ.Modulus[i]
		twoqi := qi << 1
		p1tmp := p1Q.Coeffs[i]
		p2tmp := p2Q.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := qi - modDownParams[levelP][i]
		mredParams := ringQ.MredParams[i]
		bredParams := ringQ.BredParams[i]
		nttPsi := ringQ.NttPsi[i]

		// First we switch back the relevant polypool CRT array back to the NTT domain
		NTTLazy(p3tmp, p3tmp, ringQ.N, nttPsi, qi, mredParams, bredParams)

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := 0; j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(y[0]+twoqi-x[0], params, qi, mredParams)
			z[1] = MRed(y[1]+twoqi-x[1], params, qi, mredParams)
			z[2] = MRed(y[2]+twoqi-x[2], params, qi, mredParams)
			z[3] = MRed(y[3]+twoqi-x[3], params, qi, mredParams)
			z[4] = MRed(y[4]+twoqi-x[4], params, qi, mredParams)
			z[5] = MRed(y[5]+twoqi-x[5], params, qi, mredParams)
			z[6] = MRed(y[6]+twoqi-x[6], params, qi, mredParams)
			z[7] = MRed(y[7]+twoqi-x[7], params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownQPtoP reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....QlevelQ} and {P0,P1...PlevelP},
// it reduces its basis from {Q0,Q1....QlevelQ} and {P0,P1...PlevelP} to {P0,P1...PlevelP}
// and does a floored integer division of the result by Q.
func (basisextender *FastBasisExtender) ModDownQPtoP(levelQ, levelP int, p1Q, p1P, p2P *Poly) {

	ringP := basisextender.ringP
	modDownParams := basisextender.modDownparamsQtoP
	polypool := basisextender.polypoolP

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	basisextender.ModUpQtoP(levelQ, levelP, p1Q, polypool)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := 0; i < levelP+1; i++ {

		qi := ringP.Modulus[i]
		twoqi := qi << 1
		p1tmp := p1P.Coeffs[i]
		p2tmp := p2P.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := qi - modDownParams[levelP][i]
		mredParams := ringP.MredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := 0; j < ringP.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(y[0]+twoqi-x[0], params, qi, mredParams)
			z[1] = MRed(y[1]+twoqi-x[1], params, qi, mredParams)
			z[2] = MRed(y[2]+twoqi-x[2], params, qi, mredParams)
			z[3] = MRed(y[3]+twoqi-x[3], params, qi, mredParams)
			z[4] = MRed(y[4]+twoqi-x[4], params, qi, mredParams)
			z[5] = MRed(y[5]+twoqi-x[5], params, qi, mredParams)
			z[6] = MRed(y[6]+twoqi-x[6], params, qi, mredParams)
			z[7] = MRed(y[7]+twoqi-x[7], params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// Caution, returns the values in [0, 2q-1]
func modUpExact(p1, p2 [][]uint64, params modupParams) {

	var v [8]uint64
	var y0, y1, y2, y3, y4, y5, y6, y7 [32]uint64

	// We loop over each coefficient and apply the basis extension
	for x := 0; x < len(p1[0]); x = x + 8 {

		reconstructRNS(len(p1), x, p1, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, params.Q, params.mredParamsQ, params.qibMont)

		for j := 0; j < len(p2); j++ {

			pj := params.P[j]
			qInv := params.mredParamsP[j]
			qpjInv := params.qpjInv[j]
			qispjMont := params.qispjMont[j]

			res := (*[8]uint64)(unsafe.Pointer(&p2[j][x]))

			multSum(res, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, len(p1), pj, qInv, qpjInv, qispjMont)
		}
	}
}

// Decomposer is a structure that stores the parameters of the arbitrary decomposer.
// This decomposer takes a p(x)_Q (in basis Q) and returns p(x) mod qi in basis QP, where
// qi = prod(Q_i) for 0<=i<=L, where L is the number of factors in P.
type Decomposer struct {
	nQprimes    int
	nPprimes    int
	alpha       int
	beta        int
	xalpha      []int
	modUpParams [][]modupParams
	QInt        *big.Int
	PInt        *big.Int
}

// Xalpha returns a slice that contains all the values of #Qi/#Pi.
func (decomposer *Decomposer) Xalpha() (xalpha []int) {
	return decomposer.xalpha
}

// NewDecomposer creates a new Decomposer.
func NewDecomposer(Q, P []uint64) (decomposer *Decomposer) {
	decomposer = new(Decomposer)

	decomposer.nQprimes = len(Q)
	decomposer.nPprimes = len(P)

	decomposer.QInt = NewUint(1)
	for i := range Q {
		decomposer.QInt.Mul(decomposer.QInt, NewUint(Q[i]))
	}

	decomposer.PInt = NewUint(1)
	for i := range P {
		decomposer.PInt.Mul(decomposer.PInt, NewUint(P[i]))
	}

	decomposer.alpha = len(P)
	decomposer.beta = int(math.Ceil(float64(len(Q)) / float64(decomposer.alpha)))

	decomposer.xalpha = make([]int, decomposer.beta)
	for i := range decomposer.xalpha {
		decomposer.xalpha[i] = decomposer.alpha
	}

	if len(Q)%decomposer.alpha != 0 {
		decomposer.xalpha[decomposer.beta-1] = len(Q) % decomposer.alpha
	}

	decomposer.modUpParams = make([][]modupParams, decomposer.beta)

	// Create a basis extension for each possible combination of [Qi,Pj] according to xalpha
	for i := 0; i < decomposer.beta; i++ {

		decomposer.modUpParams[i] = make([]modupParams, decomposer.xalpha[i]-1)

		for j := 0; j < decomposer.xalpha[i]-1; j++ {

			Qi := make([]uint64, j+2)
			Pi := make([]uint64, len(Q)+len(P))

			for k := 0; k < j+2; k++ {
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

// DecomposeAndSplit decomposes a polynomial p(x) in basis Q, reduces it modulo qi, and returns
// the result in basis QP separately.
func (decomposer *Decomposer) DecomposeAndSplit(levelQ, levelP, crtDecompLevel int, p0Q, p1Q, p1P *Poly) {

	alphai := decomposer.xalpha[crtDecompLevel]

	p0Qidxst := crtDecompLevel * decomposer.alpha
	p0Qidxed := p0Qidxst + alphai

	// First we check if the vector can simply by coping and rearranging elements (the case where no reconstruction is needed)
	if (p0Qidxed > levelQ+1 && (levelQ+1)%(levelP+1) == 1) || alphai == 1 {

		for j := 0; j < levelQ+1; j++ {
			copy(p1Q.Coeffs[j], p0Q.Coeffs[p0Qidxst])
		}

		for j := 0; j < decomposer.nPprimes; j++ {
			copy(p1P.Coeffs[j], p0Q.Coeffs[p0Qidxst])
		}

		// Otherwise, we apply a fast exact base conversion for the reconstruction
	} else {

		var index int
		if levelQ >= alphai+crtDecompLevel*decomposer.alpha {
			index = decomposer.xalpha[crtDecompLevel] - 2
		} else {
			index = (levelQ - 1) % decomposer.alpha
		}

		params := decomposer.modUpParams[crtDecompLevel][index]

		var v [8]uint64
		var vi [8]float64
		var y0, y1, y2, y3, y4, y5, y6, y7 [32]uint64
		var qibMont, qi, pj, mredParams uint64
		var qif float64

		// We loop over each coefficient and apply the basis extension
		for x := 0; x < len(p0Q.Coeffs[0]); x = x + 8 {

			vi[0], vi[1], vi[2], vi[3], vi[4], vi[5], vi[6], vi[7] = 0, 0, 0, 0, 0, 0, 0, 0

			// Coefficients to be decomposed
			for i, j := 0, p0Qidxst; i < index+2; i, j = i+1, j+1 {

				qibMont = params.qibMont[i]
				qi = params.Q[i]
				mredParams = params.mredParamsQ[i]
				qif = float64(qi)

				px := (*[8]uint64)(unsafe.Pointer(&p0Q.Coeffs[j][x]))
				py := (*[8]uint64)(unsafe.Pointer(&p1Q.Coeffs[j][x]))

				// For the coefficients to be decomposed, we can simply copy them
				py[0], py[1], py[2], py[3], py[4], py[5], py[6], py[7] = px[0], px[1], px[2], px[3], px[4], px[5], px[6], px[7]

				y0[i] = MRed(px[0], qibMont, qi, mredParams)
				y1[i] = MRed(px[1], qibMont, qi, mredParams)
				y2[i] = MRed(px[2], qibMont, qi, mredParams)
				y3[i] = MRed(px[3], qibMont, qi, mredParams)
				y4[i] = MRed(px[4], qibMont, qi, mredParams)
				y5[i] = MRed(px[5], qibMont, qi, mredParams)
				y6[i] = MRed(px[6], qibMont, qi, mredParams)
				y7[i] = MRed(px[7], qibMont, qi, mredParams)

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

			// Coefficients of index smaller than the ones to be decomposed
			for j := 0; j < p0Qidxst; j++ {

				pj = params.P[j]
				qInv := params.mredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]

				res := (*[8]uint64)(unsafe.Pointer(&p1Q.Coeffs[j][x]))

				multSum(res, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, index+2, pj, qInv, qpjInv, qispjMont)
			}

			// Coefficients of index greater than the ones to be decomposed
			for j := decomposer.alpha * crtDecompLevel; j < levelQ+1; j = j + 1 {

				pj = params.P[j]
				qInv := params.mredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]

				res := (*[8]uint64)(unsafe.Pointer(&p1Q.Coeffs[j][x]))

				multSum(res, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, index+2, pj, qInv, qpjInv, qispjMont)
			}

			// Coefficients of the special primes Pi
			for j, u := 0, decomposer.nQprimes; j < levelP+1; j, u = j+1, u+1 {

				pj = params.P[u]
				qInv := params.mredParamsP[u]
				qpjInv := params.qpjInv[u]
				qispjMont := params.qispjMont[u]

				res := (*[8]uint64)(unsafe.Pointer(&p1P.Coeffs[j][x]))

				multSum(res, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, index+2, pj, qInv, qpjInv, qispjMont)
			}
		}
	}
}

func reconstructRNS(index, x int, p [][]uint64, v *[8]uint64, y0, y1, y2, y3, y4, y5, y6, y7 *[32]uint64, Q, QInv, QbMont []uint64) {

	var vi [8]float64
	var qi, qiInv, qibMont uint64
	var qif float64

	for i := 0; i < index; i++ {

		qibMont = QbMont[i]
		qi = Q[i]
		qiInv = QInv[i]
		qif = float64(qi)
		pTmp := (*[8]uint64)(unsafe.Pointer(&p[i][x]))

		y0[i] = MRed(pTmp[0], qibMont, qi, qiInv)
		y1[i] = MRed(pTmp[1], qibMont, qi, qiInv)
		y2[i] = MRed(pTmp[2], qibMont, qi, qiInv)
		y3[i] = MRed(pTmp[3], qibMont, qi, qiInv)
		y4[i] = MRed(pTmp[4], qibMont, qi, qiInv)
		y5[i] = MRed(pTmp[5], qibMont, qi, qiInv)
		y6[i] = MRed(pTmp[6], qibMont, qi, qiInv)
		y7[i] = MRed(pTmp[7], qibMont, qi, qiInv)

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

	v[0] = uint64(vi[0])
	v[1] = uint64(vi[1])
	v[2] = uint64(vi[2])
	v[3] = uint64(vi[3])
	v[4] = uint64(vi[4])
	v[5] = uint64(vi[5])
	v[6] = uint64(vi[6])
	v[7] = uint64(vi[7])
}

// Caution, returns the values in [0, 2q-1]
func multSum(res, v *[8]uint64, y0, y1, y2, y3, y4, y5, y6, y7 *[32]uint64, index int, pj, qInv uint64, qpjInv, qispjMont []uint64) {

	var rlo, rhi [8]uint64
	var mhi, mlo, c, hhi uint64

	// Accumulates the sum on uint128 and does a lazy montgomery reduction at the end
	for i := 0; i < index; i++ {

		mhi, mlo = bits.Mul64(y0[i], qispjMont[i])
		rlo[0], c = bits.Add64(rlo[0], mlo, 0)
		rhi[0] += mhi + c

		mhi, mlo = bits.Mul64(y1[i], qispjMont[i])
		rlo[1], c = bits.Add64(rlo[1], mlo, 0)
		rhi[1] += mhi + c

		mhi, mlo = bits.Mul64(y2[i], qispjMont[i])
		rlo[2], c = bits.Add64(rlo[2], mlo, 0)
		rhi[2] += mhi + c

		mhi, mlo = bits.Mul64(y3[i], qispjMont[i])
		rlo[3], c = bits.Add64(rlo[3], mlo, 0)
		rhi[3] += mhi + c

		mhi, mlo = bits.Mul64(y4[i], qispjMont[i])
		rlo[4], c = bits.Add64(rlo[4], mlo, 0)
		rhi[4] += mhi + c

		mhi, mlo = bits.Mul64(y5[i], qispjMont[i])
		rlo[5], c = bits.Add64(rlo[5], mlo, 0)
		rhi[5] += mhi + c

		mhi, mlo = bits.Mul64(y6[i], qispjMont[i])
		rlo[6], c = bits.Add64(rlo[6], mlo, 0)
		rhi[6] += mhi + c

		mhi, mlo = bits.Mul64(y7[i], qispjMont[i])
		rlo[7], c = bits.Add64(rlo[7], mlo, 0)
		rhi[7] += mhi + c
	}

	hhi, _ = bits.Mul64(rlo[0]*qInv, pj)
	res[0] = rhi[0] - hhi + pj + qpjInv[v[0]]

	hhi, _ = bits.Mul64(rlo[1]*qInv, pj)
	res[1] = rhi[1] - hhi + pj + qpjInv[v[1]]

	hhi, _ = bits.Mul64(rlo[2]*qInv, pj)
	res[2] = rhi[2] - hhi + pj + qpjInv[v[2]]

	hhi, _ = bits.Mul64(rlo[3]*qInv, pj)
	res[3] = rhi[3] - hhi + pj + qpjInv[v[3]]

	hhi, _ = bits.Mul64(rlo[4]*qInv, pj)
	res[4] = rhi[4] - hhi + pj + qpjInv[v[4]]

	hhi, _ = bits.Mul64(rlo[5]*qInv, pj)
	res[5] = rhi[5] - hhi + pj + qpjInv[v[5]]

	hhi, _ = bits.Mul64(rlo[6]*qInv, pj)
	res[6] = rhi[6] - hhi + pj + qpjInv[v[6]]

	hhi, _ = bits.Mul64(rlo[7]*qInv, pj)
	res[7] = rhi[7] - hhi + pj + qpjInv[v[7]]
}
