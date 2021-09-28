package ring

import (
	"math"
	"math/big"
)

//FastBasisExtender Algorithm from https://eprint.iacr.org/2018/117.pdf
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

// NewFastBasisExtender is a struct storing all the pre-computed parameters to apply modUp from a basis Q to P
// and modDown from a basis QP to Q.
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
	params.qispjMont = make([][]uint64, len(Q))
	for i, qi := range Q {

		QiB.SetUint64(qi)
		QiStar.Quo(modulusbigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)

		// (Q/Qi)^-1) * r (mod Qi) (in Montgomery form)
		params.qibMont[i] = MForm(QiBarre.Uint64(), qi, params.bredParamsQ[i])

		params.qispjMont[i] = make([]uint64, len(P))
		for j, pj := range P {
			// (Q/qi * r) (mod Pj) (in Montgomery form)
			params.qispjMont[i][j] = MForm(tmp.Mod(QiStar, NewUint(pj)).Uint64(), pj, params.bredParamsP[j])
		}
	}

	var v uint64

	params.qpjInv = make([][]uint64, len(P))
	for j, pj := range P {
		params.qpjInv[j] = make([]uint64, len(Q)+1)
		// Correction Term (v*Q) mod each Pj
		v = pj - tmp.Mod(modulusbigint, NewUint(pj)).Uint64()
		params.qpjInv[j][0] = 0
		for i := 1; i < len(Q)+1; i++ {
			params.qpjInv[j][i] = CRed(params.qpjInv[j][i-1]+v, pj)
		}
	}

	return params
}

// ModUpSplitQP extends the basis of a polynomial
// Given a polynomial with coefficients in basis {Q0,Q1....Qi}
// Extends its basis from {Q0,Q1....Qi} to {Q0,Q1....Qi,P0,P1...Pj}
func (basisextender *FastBasisExtender) ModUpSplitQP(level uint64, p1, p2 *Poly) {
	modUpExact(p1.Coeffs[:level+1], p2.Coeffs[:uint64(len(basisextender.paramsQP.P))], basisextender.paramsQP)
}

// ModUpSplitPQ extends the basis of a polynomial
// Given a polynomial with coefficients in basis {P0,P1....Pi}
// Extends its basis from {P0,P1....Pi} to {Q0,Q1...Qj}
func (basisextender *FastBasisExtender) ModUpSplitPQ(level uint64, p1, p2 *Poly) {
	modUpExact(p1.Coeffs[:level+1], p2.Coeffs[:uint64(len(basisextender.paramsPQ.P))], basisextender.paramsPQ)
}

// ModDownNTTPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi,P0,P1...Pj}
// Reduces its basis from {Q0,Q1....Qi,P0,P1...Pj} to {Q0,Q1....Qi}
// and does a runded integer division of the result by P.
// Inputs must be in the NTT domain.
func (basisextender *FastBasisExtender) ModDownNTTPQ(level uint64, p1, p2 *Poly) {

	contextQ := basisextender.contextQ
	contextP := basisextender.contextP
	modDownParams := basisextender.modDownParamsPQ
	polypool := basisextender.polypoolQ

	// First we get the P basis part of p1 out of the NTT domain
	for j := 0; j < len(contextP.Modulus); j++ {
		InvNTT(p1.Coeffs[len(contextQ.Modulus)+j], p1.Coeffs[len(contextQ.Modulus)+j], contextP.N, contextP.GetNttPsiInv()[j], contextP.GetNttNInv()[j], contextP.Modulus[j], contextP.GetMredParams()[j])
	}

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	modUpExact(p1.Coeffs[len(contextQ.Modulus):len(contextQ.Modulus)+len(contextP.Modulus)], polypool.Coeffs[:level+1], basisextender.paramsPQ)

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
		for j := uint64(0); j < contextQ.N; j++ {
			p2tmp[j] = MRed(p1tmp[j]+(qi-p3tmp[j]), params, qi, mredParams)
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
		for j := uint64(0); j < contextQ.N; j++ {
			p2tmp[j] = MRed(p1tmp[j]+(qi-p3tmp[j]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi,P0,P1...Pj}
// Reduces its basis from {Q0,Q1....Qi,P0,P1...Pj} to {Q0,Q1....Qi}
// and does a runded integer division of the result by P.
func (basisextender *FastBasisExtender) ModDownPQ(level uint64, p1, p2 *Poly) {

	context := basisextender.contextQ
	modDownParams := basisextender.modDownParamsPQ
	polypool := basisextender.polypoolQ

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	modUpExact(p1.Coeffs[level+1:level+1+uint64(len(basisextender.paramsQP.P))], polypool.Coeffs[:level+1], basisextender.paramsPQ)

	// Finaly, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := context.Modulus[i]
		p1tmp := p1.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := modDownParams[i]
		mredParams := context.mredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = MRed(p1tmp[j]+(qi-p3tmp[j]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownSplitedPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi} and {P0,P1...Pj}
// Reduces its basis from {Q0,Q1....Qi} and {P0,P1...Pj} to {Q0,Q1....Qi}
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
		for j := uint64(0); j < contextQ.N; j++ {
			p2tmp[j] = MRed(p1tmp[j]+(qi-p3tmp[j]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownSplitedQP reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi} and {P0,P1...Pj}
// Reduces its basis from {Q0,Q1....Qi} and {P0,P1...Pj} to {P0,P1...Pj}
// and does a runded integer division of the result by Q.
func (basisextender *FastBasisExtender) ModDownSplitedQP(levelQ, levelP uint64, p1Q, p1P, p2 *Poly) {

	contextP := basisextender.contextP
	//contextQ := basisextender.contextQ
	modDownParams := basisextender.modDownParamsQP
	polypool := basisextender.polypoolP

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)

	//bigint_coeffs := make([]*big.Int, contextP.N)
	//contextQ.PolyToBigint(p1Q, bigint_coeffs)
	//fmt.Println("Q", bigint_coeffs[0])

	basisextender.ModUpSplitQP(levelQ, p1Q, polypool)

	//contextP.PolyToBigint(polypool, bigint_coeffs)
	//fmt.Println("P", bigint_coeffs[0])

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

	var v uint64
	var vi float64
	var xpj uint64

	y := make([]uint64, len(p1))

	//We loop over each coefficient and apply the basis extension
	for x := uint64(0); x < uint64(len(p1[0])); x++ {

		vi = 0

		for i := 0; i < len(p1); i++ {

			y[i] = MRed(p1[i][x], params.qibMont[i], params.Q[i], params.mredParamsQ[i])

			// Computation of the correction term v * Q%pi
			vi += float64(y[i]) / float64(params.Q[i])

		}

		// Index of the correction term
		v = uint64(vi)

		for j := 0; j < len(p2); j++ {

			xpj = 0

			for i := 0; i < len(p1); i++ {
				xpj += MRed(y[i], params.qispjMont[i][j], params.P[j], params.mredParamsP[j])

				if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
					xpj = BRedAdd(xpj, params.P[j], params.bredParamsP[j])
				}
			}

			p2[j][x] = BRedAdd(xpj+params.qpjInv[j][v], params.P[j], params.bredParamsP[j])

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

	var v uint64
	var vi float64
	var xpj uint64

	//fmt.Println(p0idxed, level + 1, (level+1)%decomposer.nPprimes)

	// First we check if the vector can simply by coping and rearanging elements (the case where no reconstruction is needed)
	if (p0idxed > level+1 && (level+1)%decomposer.nPprimes == 1) || alphai == 1 {

		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x++ {

			for j := uint64(0); j < level+decomposer.nPprimes+1; j++ {

				p1.Coeffs[j][x] = p0.Coeffs[p0idxst][x]
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

		/*
			fmt.Println()
			fmt.Println("CRT DECOMP :", crtDecompLevel)
			fmt.Println("Min threshold :", alphai + crtDecompLevel*decomposer.alpha)
			fmt.Println("INDEX :", index)
			fmt.Println("Level :", level)
			fmt.Println("xalpha :", decomposer.xalpha)
			fmt.Println()
		*/

		params := decomposer.modUpParams[crtDecompLevel][index]

		y := make([]uint64, index+2)

		//We loop over each coefficient and apply the basis extension
		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x++ {

			vi = 0

			// Coefficients to be decomposed
			for i := range y {

				// For the coefficients to be decomposed, we can simplly copy them
				p1.Coeffs[i+int(p0idxst)][x] = p0.Coeffs[i+int(p0idxst)][x]

				y[i] = MRed(p0.Coeffs[i+int(p0idxst)][x], params.qibMont[i], params.Q[i], params.mredParamsQ[i])

				// Computation of the correction term v * Q%pi
				vi += float64(y[i]) / float64(params.Q[i])
			}

			// Index of the correction term
			v = uint64(vi)

			// Coefficients of index smaller than the ones to be decomposer
			for j := uint64(0); j < p0idxst; j++ {

				xpj = 0

				for i := range y {
					xpj += MRed(y[i], params.qispjMont[i][j], params.P[j], params.mredParamsP[j])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = BRedAdd(xpj, params.P[j], params.bredParamsP[j])
					}
				}

				p1.Coeffs[j][x] = BRedAdd(xpj+params.qpjInv[j][v], params.P[j], params.bredParamsP[j])

			}

			// Coefficients of index greater than the ones to be decomposer
			for j := decomposer.alpha * crtDecompLevel; j < level+1; j = j + 1 {

				xpj = 0

				for i := range y {
					xpj += MRed(y[i], params.qispjMont[i][j], params.P[j], params.mredParamsP[j])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = BRedAdd(xpj, params.P[j], params.bredParamsP[j])
					}
				}

				p1.Coeffs[j][x] = BRedAdd(xpj+params.qpjInv[j][v], params.P[j], params.bredParamsP[j])

			}

			// Coefficients of the special primes
			for u, j := decomposer.nQprimes, level+1; j < level+1+decomposer.nPprimes; u, j = u+1, j+1 {

				xpj = 0

				for i := range y {
					xpj += MRed(y[i], params.qispjMont[i][u], params.P[u], params.mredParamsP[u])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = BRedAdd(xpj, params.P[u], params.bredParamsP[u])
					}
				}

				p1.Coeffs[j][x] = BRedAdd(xpj+params.qpjInv[u][v], params.P[u], params.bredParamsP[u])
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

	var v uint64
	var vi float64
	var xpj uint64

	// First we check if the vector can simply by coping and rearanging elements (the case where no reconstruction is needed)
	if (p0idxed > level+1 && (level+1)%decomposer.nPprimes == 1) || alphai == 1 {

		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x++ {

			for j := uint64(0); j < level+1; j++ {

				p1Q.Coeffs[j][x] = p0.Coeffs[p0idxst][x]
			}

			for j := uint64(0); j < decomposer.nPprimes; j++ {

				p1P.Coeffs[j][x] = p0.Coeffs[p0idxst][x]
			}
		}

		// Else we apply a fast exact base conversion for the reconstruction
	} else {

		var index uint64
		if level >= alphai+crtDecompLevel*decomposer.alpha {
			index = decomposer.xalpha[crtDecompLevel] - 2
		} else {
			index = (level - 1) % decomposer.alpha
		}

		params := decomposer.modUpParams[crtDecompLevel][index]

		y := make([]uint64, index+2)

		//We loop over each coefficient and apply the basis extension
		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x++ {

			vi = 0

			// Coefficients to be decomposed
			for i := range y {

				// For the coefficients to be decomposed, we can simplly copy them
				p1Q.Coeffs[i+int(p0idxst)][x] = p0.Coeffs[i+int(p0idxst)][x]

				y[i] = MRed(p0.Coeffs[i+int(p0idxst)][x], params.qibMont[i], params.Q[i], params.mredParamsQ[i])

				// Computation of the correction term v * Q%pi
				vi += float64(y[i]) / float64(params.Q[i])
			}

			// Index of the correction term
			v = uint64(vi)

			// Coefficients of index smaller than the ones to be decomposer
			for j := uint64(0); j < p0idxst; j++ {

				xpj = 0

				for i := range y {
					xpj += MRed(y[i], params.qispjMont[i][j], params.P[j], params.mredParamsP[j])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = BRedAdd(xpj, params.P[j], params.bredParamsP[j])
					}
				}

				p1Q.Coeffs[j][x] = BRedAdd(xpj+params.qpjInv[j][v], params.P[j], params.bredParamsP[j])

			}

			// Coefficients of index greater than the ones to be decomposer
			for j := decomposer.alpha * crtDecompLevel; j < level+1; j = j + 1 {

				xpj = 0

				for i := range y {
					xpj += MRed(y[i], params.qispjMont[i][j], params.P[j], params.mredParamsP[j])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = BRedAdd(xpj, params.P[j], params.bredParamsP[j])
					}
				}

				p1Q.Coeffs[j][x] = BRedAdd(xpj+params.qpjInv[j][v], params.P[j], params.bredParamsP[j])

			}

			// Coefficients of the special primes
			for j, u := uint64(0), decomposer.nQprimes; j < decomposer.nPprimes; j, u = j+1, u+1 {

				xpj = 0

				for i := range y {
					xpj += MRed(y[i], params.qispjMont[i][u], params.P[u], params.mredParamsP[u])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = BRedAdd(xpj, params.P[u], params.bredParamsP[u])
					}
				}

				p1P.Coeffs[j][x] = BRedAdd(xpj+params.qpjInv[u][v], params.P[u], params.bredParamsP[u])
			}
		}
	}
}
