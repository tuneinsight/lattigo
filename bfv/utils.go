package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math"
)

type ArbitraryDecomposer struct {
	nQprimes    uint64
	nPprimes    uint64
	alpha       uint64
	beta        uint64
	xalpha      []uint64
	modUpParams [][]*modupParams
	Q_bigint    *ring.Int
	P_bigint    *ring.Int
}

func NewArbitraryDecomposer(Q, P []uint64) (decomposer *ArbitraryDecomposer) {
	decomposer = new(ArbitraryDecomposer)

	decomposer.nQprimes = uint64(len(Q))
	decomposer.nPprimes = uint64(len(P))

	decomposer.Q_bigint = ring.NewUint(1)
	for i := range Q {
		decomposer.Q_bigint.Mul(decomposer.Q_bigint, ring.NewUint(Q[i]))
	}

	decomposer.P_bigint = ring.NewUint(1)
	for i := range P {
		decomposer.P_bigint.Mul(decomposer.P_bigint, ring.NewUint(P[i]))
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

//Algorithm from https://eprint.iacr.org/2018/117.pdf
type FastBasisExtender struct {
	paramsQP *modupParams
	paramsPQ *modupParams
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

func NewFastBasisExtender(Q, P []uint64) *FastBasisExtender {

	newParams := new(FastBasisExtender)

	newParams.paramsQP = basisextenderparameters(Q, P)
	newParams.paramsPQ = basisextenderparameters(P, Q)

	return newParams
}

func basisextenderparameters(Q, P []uint64) (params *modupParams) {

	params = new(modupParams)

	params.Q = make([]uint64, len(Q))
	params.bredParamsQ = make([][]uint64, len(Q))
	params.mredParamsQ = make([]uint64, len(Q))
	for i, qi := range Q {
		params.Q[i] = Q[i]
		params.bredParamsQ[i] = ring.BRedParams(qi)
		params.mredParamsQ[i] = ring.MRedParams(qi)
	}

	params.P = make([]uint64, len(P))
	params.bredParamsP = make([][]uint64, len(P))
	params.mredParamsP = make([]uint64, len(P))

	for i, pj := range P {
		params.P[i] = P[i]
		params.bredParamsP[i] = ring.BRedParams(pj)
		params.mredParamsP[i] = ring.MRedParams(pj)
	}

	tmp := new(ring.Int)
	QiB := new(ring.Int)
	QiStar := new(ring.Int)
	QiBarre := new(ring.Int)

	modulusbigint := ring.NewUint(1)
	for _, qi := range Q {
		modulusbigint.Mul(modulusbigint, ring.NewUint(qi))
	}

	params.qibMont = make([]uint64, len(Q))
	params.qispjMont = make([][]uint64, len(Q))
	for i, qi := range Q {

		QiB.SetUint(qi)
		QiStar.Div(modulusbigint, QiB)
		QiBarre.Inv(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)

		// (Q/Qi)^-1) * r (mod Qi) (in Montgomery form)
		params.qibMont[i] = ring.MForm(QiBarre.Uint64(), qi, params.bredParamsQ[i])

		params.qispjMont[i] = make([]uint64, len(P))
		for j, pj := range P {
			// (Q/qi * r) (mod Pj) (in Montgomery form)
			params.qispjMont[i][j] = ring.MForm(tmp.Mod(QiStar, ring.NewUint(pj)).Uint64(), pj, params.bredParamsP[j])
		}
	}

	var v uint64

	params.qpjInv = make([][]uint64, len(P))
	for j, pj := range P {
		params.qpjInv[j] = make([]uint64, len(Q)+1)
		// Correction Term (v*Q) mod each Pj
		v = pj - tmp.Mod(modulusbigint, ring.NewUint(pj)).Uint64()
		params.qpjInv[j][0] = 0
		for i := 1; i < len(Q)+1; i++ {
			params.qpjInv[j][i] = ring.CRed(params.qpjInv[j][i-1]+v, pj)
		}
	}

	return params
}

// Extends the basis of a ring
// Given a ring with coefficients in basis {Q0,Q1....Qi}
// Extends its basis from {Q0,Q1....Qi} to {Q0,Q1....Qi,P0,P1...Pj}
func (basisextender *FastBasisExtender) ModUp(level uint64, p1, p2 *ring.Poly) {

	if p1 != p2 {
		for i := range p1.Coeffs {
			for j := uint64(0); j < level; j++ {
				p2.Coeffs[i][j] = p1.Coeffs[i][j]
			}
		}
	}

	modUpExact(p1.Coeffs[:level+1], p2.Coeffs[level+1:level+1+uint64(len(basisextender.paramsQP.P))], basisextender.paramsQP)
}

func (basisextender *FastBasisExtender) ModDown(context *ring.Context, rescalParamsKeys []uint64, level uint64, p1, p2, polypool *ring.Poly) {

	// First we get the P basis part of p1 out of the NTT domain
	for j := uint64(0); j < uint64(len(basisextender.paramsQP.P)); j++ {
		ring.InvNTT(p1.Coeffs[level+1+j], p1.Coeffs[level+1+j], context.N, context.GetNttPsiInv()[level+1+j], context.GetNttNInv()[level+1+j], context.Modulus[level+1+j], context.GetMredParams()[level+1+j])
	}

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	modUpExact(p1.Coeffs[level+1:level+1+uint64(len(basisextender.paramsQP.P))], polypool.Coeffs[:level+1], basisextender.paramsPQ)

	// Finaly, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		// First we switch back the relevant polypool CRT array back to the NTT domain
		ring.NTT(polypool.Coeffs[i], polypool.Coeffs[i], context.N, context.GetNttPsi()[i], context.Modulus[i], context.GetMredParams()[i], context.GetBredParams()[i])

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = ring.MRed(p1.Coeffs[i][j]+(context.Modulus[i]-polypool.Coeffs[i][j]), rescalParamsKeys[i], context.Modulus[i], context.GetMredParams()[i])
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

func (decomposer *ArbitraryDecomposer) Decompose(level, crtDecompLevel uint64, p0, p1 *ring.Poly) {

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

				y[i] = ring.MRed(p0.Coeffs[i+int(p0idxst)][x], params.qibMont[i], params.Q[i], params.mredParamsQ[i])

				// Computation of the correction term v * Q%pi
				vi += float64(y[i]) / float64(params.Q[i])
			}

			// Index of the correction term
			v = uint64(vi)

			// Coefficients of index smaller than the ones to be decomposer
			for j := uint64(0); j < p0idxst; j++ {

				xpj = 0

				for i := range y {
					xpj += ring.MRed(y[i], params.qispjMont[i][j], params.P[j], params.mredParamsP[j])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = ring.BRedAdd(xpj, params.P[j], params.bredParamsP[j])
					}
				}

				p1.Coeffs[j][x] = ring.BRedAdd(xpj+params.qpjInv[j][v], params.P[j], params.bredParamsP[j])

			}

			// Coefficients of index greater than the ones to be decomposer
			for j := decomposer.alpha * crtDecompLevel; j < level+1; j = j + 1 {

				xpj = 0

				for i := range y {
					xpj += ring.MRed(y[i], params.qispjMont[i][j], params.P[j], params.mredParamsP[j])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = ring.BRedAdd(xpj, params.P[j], params.bredParamsP[j])
					}
				}

				p1.Coeffs[j][x] = ring.BRedAdd(xpj+params.qpjInv[j][v], params.P[j], params.bredParamsP[j])

			}

			// Coefficients of the special primes
			for u, j := decomposer.nQprimes, level+1; j < level+1+decomposer.nPprimes; u, j = u+1, j+1 {

				xpj = 0

				for i := range y {
					xpj += ring.MRed(y[i], params.qispjMont[i][u], params.P[u], params.mredParamsP[u])

					if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
						xpj = ring.BRedAdd(xpj, params.P[u], params.bredParamsP[u])
					}
				}

				p1.Coeffs[j][x] = ring.BRedAdd(xpj+params.qpjInv[u][v], params.P[u], params.bredParamsP[u])
			}
		}
	}
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

			y[i] = ring.MRed(p1[i][x], params.qibMont[i], params.Q[i], params.mredParamsQ[i])

			// Computation of the correction term v * Q%pi
			vi += float64(y[i]) / float64(params.Q[i])

		}

		// Index of the correction term
		v = uint64(vi)

		for j := 0; j < len(p2); j++ {

			xpj = 0

			for i := 0; i < len(p1); i++ {
				xpj += ring.MRed(y[i], params.qispjMont[i][j], params.P[j], params.mredParamsP[j])

				if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
					xpj = ring.BRedAdd(xpj, params.P[j], params.bredParamsP[j])
				}
			}

			p2[j][x] = ring.BRedAdd(xpj+params.qpjInv[j][v], params.P[j], params.bredParamsP[j])

		}
	}
}

// equalslice compares two slices of uint64 values, and return true if they are equal, else false.
func EqualSlice(a, b []uint64) bool {

	if len(a) != len(b) {
		return false
	}

	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// min returns the minimum value of the input slice of uint64 values.
func min(a, b uint64) (r uint64) {
	if a <= b {
		return a
	}
	return b
}

// max returns the maximum value of the input slice of uint64 values.
func max(a, b uint64) (r uint64) {
	if a >= b {
		return a
	}
	return b
}

// bitReverse64 returns the bit-reverse value of the input value, within a context of 2^bitLen.
func bitReverse64(index, bitLen uint64) uint64 {
	indexReverse := uint64(0)
	for i := uint64(0); i < bitLen; i++ {
		if (index>>i)&1 != 0 {
			indexReverse |= 1 << (bitLen - 1 - i)
		}
	}
	return indexReverse
}

// hammingWeight64 returns the hammingweight if the input value.
func hammingWeight64(x uint64) uint64 {
	x -= (x >> 1) & 0x5555555555555555
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
	return ((x * 0x0101010101010101) & 0xffffffffffffffff) >> 56
}
