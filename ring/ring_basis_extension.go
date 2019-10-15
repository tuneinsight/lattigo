package ring

//==========================================
//===== CRT BASIS EXTENSION PARAMETERS =====
//==========================================

// BasisExtender is the structure keeping all the pre-computed constants required to apply
// a basis extension of a polynomial in basis Q to a polynomial in basis Q + P.
// Algorithm from https://eprint.iacr.org/2018/117.pdf
type BasisExtender struct {
	contextQ *Context
	contextP *Context

	//Parameters for basis extension

	// (Q/Qi)^-1) * r (mod each Qi) (in Montgomery form)
	qibMont []uint64
	// Q/qi (mod each Pj) (in Montgomery form)
	qispjMont [][]uint64
	// Q*v (mod each Pj) for v in [1,...,k] where k is the number of Pj modulies
	qpjInv [][]uint64
	// Qi as a float128 variable
	qiFloat128 []Float128
}

// NewBasisExtender creates a new BasisExtender, that will be used to extend a polynomial in basis Q to a polynomial in basis Q + P.
func NewBasisExtender(contextQ, contextP *Context) (newParams *BasisExtender) {

	newParams = new(BasisExtender)

	var PjB Int
	var QiB Int
	var QiStar Int
	var QiBarre Int

	newParams.contextQ = contextQ
	newParams.contextP = contextP

	tmp := NewUint(1)

	// Q, P in Bigint

	newParams.qibMont = make([]uint64, len(contextQ.Modulus))
	newParams.qispjMont = make([][]uint64, len(contextQ.Modulus))
	newParams.qpjInv = make([][]uint64, len(contextP.Modulus))
	newParams.qiFloat128 = make([]Float128, len(contextQ.Modulus))

	for i, qi := range contextQ.Modulus {

		QiB.SetUint(qi)
		QiStar.Div(contextQ.ModulusBigint, &QiB)
		QiBarre.Inv(&QiStar, &QiB)
		QiBarre.Mod(&QiBarre, &QiB)

		// (Q/Qi)^-1) * r (mod Qi) (in Montgomery form)
		newParams.qibMont[i] = MForm(QiBarre.Uint64(), qi, contextQ.bredParams[i])

		// (Q/qi * r) (mod Pj) (in Montgomery form)
		newParams.qispjMont[i] = make([]uint64, len(contextP.Modulus))
		for j, pj := range contextP.Modulus {
			PjB.SetUint(pj)
			newParams.qispjMont[i][j] = MForm(tmp.Mod(&QiStar, &PjB).Uint64(), pj, contextP.bredParams[j])
		}

		newParams.qiFloat128[i] = Float128SetUint64(qi)
	}

	// Correction Term (v*Q) mod each Pj
	var v uint64
	for j, pj := range contextP.Modulus {
		// [Q]_{pi}
		PjB.SetUint(pj)
		v = pj - PjB.Mod(contextQ.ModulusBigint, &PjB).Uint64()
		newParams.qpjInv[j] = make([]uint64, len(contextQ.Modulus)+1)
		newParams.qpjInv[j][0] = 0
		for i := 1; i < len(contextQ.Modulus)+1; i++ {
			newParams.qpjInv[j][i] = CRed(newParams.qpjInv[j][i-1]+v, pj)
		}
	}

	return
}

// ExtendBasis extends the basis of a polynomial from Q to Q + P.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi},
// extends its basis from {Q0,Q1....Qi} to {Q0,Q1....Qi,P0,P1...Pj}.
func (Parameters *BasisExtender) ExtendBasis(p1, p2 *Poly) {

	var v uint64
	var vi, yiFloat128 Float128
	var xpj uint64

	y := make([]uint64, len(Parameters.contextQ.Modulus))

	// If the receiver is equal to p1, then extend the number of modulies of p1
	if p1 == p2 {
		coeffsNewBase := make([][]uint64, len(Parameters.contextP.Modulus))
		for i := 0; i < len(Parameters.contextP.Modulus); i++ {
			coeffsNewBase[i] = make([]uint64, Parameters.contextQ.N)
		}
		p1.Coeffs = append(p1.Coeffs, coeffsNewBase...)
	} else { // Else copies the qi coefficients of p1 on p2
		for i := range Parameters.contextQ.Modulus {
			for j := uint64(0); j < Parameters.contextQ.N; j++ {
				p2.Coeffs[i][j] = p1.Coeffs[i][j]
			}
		}
	}

	//We loop over each coefficient and apply the basis extension
	for x := uint64(0); x < Parameters.contextQ.N; x++ {

		vi[0], vi[1] = 0, 0

		for i, qi := range Parameters.contextQ.Modulus {

			y[i] = MRed(p1.Coeffs[i][x], Parameters.qibMont[i], qi, Parameters.contextQ.mredParams[i])

			yiFloat128[0] = float64(y[i] >> 12)
			yiFloat128[1] = float64(y[i]&0xfff) / float64(4096)

			yiFloat128 = Float128Div(yiFloat128, Parameters.qiFloat128[i])

			vi = Float128Add(vi, yiFloat128)
		}

		// Index of the correction term
		v = uint64(vi[0])

		//For each Pi we sum over the Qi
		for j, pj := range Parameters.contextP.Modulus {
			xpj = 0

			for i := range Parameters.contextQ.Modulus {
				xpj += MRed(y[i], Parameters.qispjMont[i][j], pj, Parameters.contextP.mredParams[j])

				if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
					xpj = BRedAdd(xpj, pj, Parameters.contextP.bredParams[j])
				}
			}

			p2.Coeffs[j+len(Parameters.contextQ.Modulus)][x] = BRedAdd(xpj+Parameters.qpjInv[j][v], pj, Parameters.contextP.bredParams[j])
		}
	}
}

// ExtendBasisApproximate approximates the basis extension (returns the value + some multiple of Q (equal to the correction term v) in the basis QP).
// The algorithm is identical to the ExtendBasis, except it doesn't make use of the correction term v that removes the additional multiple of Q
// introduced during the basis extension in the new basis P.
func (Parameters *BasisExtender) ExtendBasisApproximate(p1, p2 *Poly) {

	var xpj uint64

	y := make([]uint64, len(Parameters.contextQ.Modulus))

	// If the receiver is equal to p1, then extend the number of modulies of p1
	if p1 == p2 {
		coeffsNewBase := make([][]uint64, len(Parameters.contextP.Modulus))
		for i := 0; i < len(Parameters.contextP.Modulus); i++ {
			coeffsNewBase[i] = make([]uint64, Parameters.contextQ.N)
		}
		p1.Coeffs = append(p1.Coeffs, coeffsNewBase...)
	} else { // Else copies the qi coefficients of p1 on p2
		for i := range Parameters.contextQ.Modulus {
			for j := uint64(0); j < Parameters.contextQ.N; j++ {
				p2.Coeffs[i][j] = p1.Coeffs[i][j]
			}
		}
	}

	//We loop over each array of coeffs
	for x := uint64(0); x < Parameters.contextQ.N; x++ {

		for i, qi := range Parameters.contextQ.Modulus {
			y[i] = MRed(p1.Coeffs[i][x], Parameters.qibMont[i], qi, Parameters.contextQ.mredParams[i])
		}

		//For each Pi we sum over the Qi
		for j, pj := range Parameters.contextP.Modulus {
			xpj = 0

			for i := range Parameters.contextQ.Modulus {

				xpj += MRed(y[i], Parameters.qispjMont[i][j], pj, Parameters.contextP.mredParams[j])

				if i&7 == 6 { //Only every 7 addition, since we add one more 60 bit integer after the loop
					xpj = BRedAdd(xpj, pj, Parameters.contextP.bredParams[j])
				}
			}

			p2.Coeffs[j+len(Parameters.contextQ.Modulus)][x] = BRedAdd(xpj, pj, Parameters.contextP.bredParams[j])
		}
	}
}
