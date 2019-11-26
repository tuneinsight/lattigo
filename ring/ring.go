package ring

import (
	"github.com/ldsec/lattigo/utils"
	"math/big"
	"math/bits"
)

// Add adds p1 to p2 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Add(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = CRed(p1tmp[j]+p2tmp[j], qi)
		}
	}
}

// AddLvl adds p1 to p2 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) AddLvl(level uint64, p1, p2, p3 *Poly) {
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = CRed(p1tmp[j]+p2tmp[j], qi)
		}
	}
}

// AddNoMod adds p1 to p2 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*Qi -1].
func (context *Context) AddNoMod(p1, p2, p3 *Poly) {
	for i := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = p1tmp[j] + p2tmp[j]
		}
	}
}

// AddNoModLvl adds p1 to p2 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*Qi -1].
func (context *Context) AddNoModLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = p1tmp[j] + p2tmp[j]
		}
	}
}

// Sub subtracts p2 to p1 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Sub(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = CRed((p1tmp[j]+qi)-p2tmp[j], qi)
		}
	}
}

// SubLvl subtracts p2 to p1 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) SubLvl(level uint64, p1, p2, p3 *Poly) {
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = CRed((p1tmp[j]+qi)-p2tmp[j], qi)
		}
	}
}

// SubNoMod subtracts p2 to p1 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*Qi -1].
func (context *Context) SubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = (p1tmp[j] + qi) - p2tmp[j]
		}
	}
}

// SubNoModLvl subtracts p2 to p1 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*Qi -1].
func (context *Context) SubNoModLvl(level uint64, p1, p2, p3 *Poly) {
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = (p1tmp[j] + qi) - p2tmp[j]
		}
	}
}

// Neg sets all coefficients of p1 to there additive inverse, returning the result on p2.
func (context *Context) Neg(p1, p2 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = qi - p1tmp[j]
		}
	}
}

// NegLvl sets all coefficients of p1 to there additive inverse, returning the result on p2.
func (context *Context) NegLvl(level uint64, p1, p2 *Poly) {
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = qi - p1tmp[j]
		}
	}
}

// Reduce applies a modular reduction over the coefficients of p1 returning the result on p2.
func (context *Context) Reduce(p1, p2 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := context.bredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = BRedAdd(p1tmp[j], qi, bredParams)
		}
	}
}

// ReduceLvl applies a modular reduction over the coefficients of p1 returning the result on p2.
func (context *Context) ReduceLvl(level uint64, p1, p2 *Poly) {
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := context.bredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = BRedAdd(p1tmp[j], qi, bredParams)
		}
	}
}

// Mod applies a modular reduction by m over the coefficients of p1, returning the result on p2.
func (context *Context) Mod(p1 *Poly, m uint64, p2 *Poly) {
	params := BRedParams(m)
	for i := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = BRedAdd(p1tmp[j], m, params)
		}
	}
}

// AND applies a logical AND of m to the coefficients of p1,  returning the result on p2.
func (context *Context) AND(p1 *Poly, m uint64, p2 *Poly) {
	for i := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = p1tmp[j] & m
		}
	}
}

// OR applies a logical OR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) OR(p1 *Poly, m uint64, p2 *Poly) {
	for i := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = p1tmp[j] | m
		}
	}
}

// XOR applies a logical XOR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) XOR(p1 *Poly, m uint64, p2 *Poly) {
	for i := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = p1tmp[j] ^ m
		}
	}
}

// MulCoeffs multiplies p1 by p2 coefficient wise with a Barrett modular reduction, returning the result on p3.
func (context *Context) MulCoeffs(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := context.bredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = BRed(p1tmp[j], p2tmp[j], qi, bredParams)
		}
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient wise with a Barret modular reduction, adding the result to p3 with modular reduction.
func (context *Context) MulCoeffsAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := context.bredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = CRed(p3tmp[j]+BRed(p1tmp[j], p2tmp[j], qi, bredParams), qi)
		}
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient wise with a Barrett modular reduction, adding the result to p3 without modular reduction.
func (context *Context) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := context.bredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] += BRed(p1tmp[j], p2tmp[j], qi, bredParams)
		}
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, returning the result on p3.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomery(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = MRed(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

func (context *Context) MulCoeffsMontgomeryLvl(level uint64, p1, p2, p3 *Poly) {
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = MRed(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, adding the result to p3.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = CRed(p3tmp[j]+MRed(p1tmp[j], p2tmp[j], qi, mredParams), qi)
		}
	}
}

func (context *Context) MulCoeffsMontgomeryAndAddLvl(level uint64, p1, p2, p3 *Poly) {
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = CRed(p3tmp[j]+MRed(p1tmp[j], p2tmp[j], qi, mredParams), qi)
		}
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, adding the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] += MRed(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

func (context *Context) MulCoeffsMontgomeryAndAddNoModLvl(level uint64, p1, p2, p3 *Poly) {
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] += MRed(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, subtractsing the result to p3 with modular reduction.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = CRed(p3tmp[j]+(qi-MRed(p1tmp[j], p2tmp[j], qi, mredParams)), qi)
		}
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, subtractsing the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = p3tmp[j] + (qi - MRed(p1tmp[j], p2tmp[j], qi, mredParams))
		}
	}
}

// MulCoeffsConstant multiplies p1 by p2 coefficient wise with a constant time Barrett modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*Qi -1].
func (context *Context) MulCoeffsConstant(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := context.bredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = BRedConstant(p1tmp[j], p2tmp[j], qi, bredParams)
		}
	}
}

// MulCoeffsMontgomeryConstant multiplies p1 by p2 coefficient wise with a constant time Montgomery modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*Qi -1].
func (context *Context) MulCoeffsMontgomeryConstant(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p3tmp[j] = MRedConstant(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulPoly multiplies p1 by p2 and returns the result on p3.
func (context *Context) MulPoly(p1, p2, p3 *Poly) {

	a := context.NewPoly()
	b := context.NewPoly()

	context.NTT(p1, a)
	context.NTT(p2, b)
	context.MulCoeffs(a, b, p3)
	context.InvNTT(p3, p3)
}

// MulPolyMontgomery multiplies p1 by p2 and returns the result on p3.
// Expect wither p1 or p2 to be in Montgomery form for correctness.
func (context *Context) MulPolyMontgomery(p1, p2, p3 *Poly) {

	a := context.NewPoly()
	b := context.NewPoly()

	context.NTT(p1, a)
	context.NTT(p2, b)
	context.MulCoeffsMontgomery(a, b, p3)
	context.InvNTT(p3, p3)
}

// MulPolyNaive multiplies p1 by p2 with a naive convolution, returning the result on p3.
func (context *Context) MulPolyNaive(p1, p2, p3 *Poly) {

	p1Copy := p1.CopyNew()
	p2Copy := p2.CopyNew()

	context.MForm(p1Copy, p1Copy)

	context.AND(p3, 0, p3)

	for x, qi := range context.Modulus {

		p1tmp, p2tmp, p3tmp := p1Copy.Coeffs[x], p2Copy.Coeffs[x], p3.Coeffs[x]

		mredParams := context.mredParams[x]

		for i := uint64(0); i < context.N; i++ {

			for j := uint64(0); j < i; j++ {
				p3tmp[j] = CRed(p3tmp[j]+(qi-MRed(p1tmp[i], p2tmp[context.N-i+j], qi, mredParams)), qi)
			}

			for j := uint64(i); j < context.N; j++ {
				p3tmp[j] = CRed(p3tmp[j]+MRed(p1tmp[i], p2tmp[j-i], qi, mredParams), qi)
			}
		}
	}
}

// MulPolyNaiveMontgomery multiplies p1 by p2 with a naive convolution, returning the result on p3.
// Much faster than MulPolyNaive.
func (context *Context) MulPolyNaiveMontgomery(p1, p2, p3 *Poly) {

	p1Copy := p1.CopyNew()
	p2Copy := p2.CopyNew()

	context.AND(p3, 0, p3)

	for x, qi := range context.Modulus {

		p1tmp, p2tmp, p3tmp := p1Copy.Coeffs[x], p2Copy.Coeffs[x], p3.Coeffs[x]

		mredParams := context.mredParams[x]

		for i := uint64(0); i < context.N; i++ {

			for j := uint64(0); j < i; j++ {
				p3tmp[j] = CRed(p3tmp[j]+(qi-MRed(p1tmp[i], p2tmp[context.N-i+j], qi, mredParams)), qi)
			}

			for j := uint64(i); j < context.N; j++ {
				p3tmp[j] = CRed(p3tmp[j]+MRed(p1tmp[i], p2tmp[j-i], qi, mredParams), qi)
			}
		}
	}
}

// Exp raises p1 to p1^e, returning the result on p2.
// TODO : implement Montgomery ladder
func (context *Context) Exp(p1 *Poly, e uint64, p2 *Poly) {

	context.NTT(p1, p1)

	tmp := context.NewPoly()
	context.Add(tmp, p1, tmp)

	for i := range context.Modulus {
		p2tmp := p2.Coeffs[i]
		for x := uint64(0); x < context.N; x++ {
			p2tmp[x] = 1
		}
	}

	for i := e; i > 0; i >>= 1 {
		if (i & 1) == 1 {
			context.MulCoeffs(p2, tmp, p2)
		}
		context.MulCoeffs(tmp, p1, tmp)
	}

	context.InvNTT(p2, p2)
	context.InvNTT(p1, p2)
}

// AddScalar adds to each coefficient of p1 a scalar and applies a modular reduction, returing the result on p2.
func (context *Context) AddScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = CRed(p1tmp[j]+scalar, Qi)
		}
	}
}

// AddScalarBigint adds to each coefficient of p1 a big.Int scalar and applies a modular reduction, returing the result on p2.
func (context *Context) AddScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	tmp := new(big.Int)
	var scalarQi uint64
	for i, Qi := range context.Modulus {
		scalarQi = tmp.Mod(scalar, NewUint(Qi)).Uint64()
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = CRed(p1tmp[j]+scalarQi, Qi)
		}
	}
}

// SubScalar subtractss to each coefficient of p1 a scalar and applies a modular reduction, returing the result on p2.
func (context *Context) SubScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = CRed(p1tmp[j]+(Qi-scalar), Qi)
		}
	}
}

// SubScalarBigint subtractss to each coefficient of p1 a big.Int scalar and applies a modular reduction, returing the result on p2.
func (context *Context) SubScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	tmp := new(big.Int)
	var scalarQi uint64
	for i, Qi := range context.Modulus {
		scalarQi = tmp.Mod(scalar, NewUint(Qi)).Uint64()
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = CRed(p1tmp[j]+(Qi-scalarQi), Qi)
		}
	}
}

// MulScalar multiplies each coefficient of p1 by a scalar and applies a modular reduction, returning the result on p2.
func (context *Context) MulScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	var scalarMont uint64
	for i, Qi := range context.Modulus {
		scalarMont = MForm(BRedAdd(scalar, Qi, context.bredParams[i]), Qi, context.bredParams[i])
		mredParams := context.mredParams[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], scalarMont, Qi, mredParams)
		}
	}
}

// MulScalarLvl multiplies each coefficient of p1 by a scalar and applies a modular reduction, returning the result on p2.
func (context *Context) MulScalarLvl(level uint64, p1 *Poly, scalar uint64, p2 *Poly) {
	var Qi, scalarMont uint64
	for i := uint64(0); i < level+1; i++ {
		Qi = context.Modulus[i]
		scalarMont = MForm(BRedAdd(scalar, Qi, context.bredParams[i]), Qi, context.bredParams[i])
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], scalarMont, Qi, mredParams)
		}
	}
}

// MulScalarBigint multiplies each coefficientsof p1 by a big.Int scalar and applies a modular reduction, returning the result on p2.
// To be used when the scalar is bigger than 64 bits.
func (context *Context) MulScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	scalarQi := new(big.Int)
	var scalarMont uint64
	for i, Qi := range context.Modulus {
		scalarQi.Mod(scalar, NewUint(Qi))
		scalarMont = MForm(BRedAdd(scalarQi.Uint64(), Qi, context.bredParams[i]), Qi, context.bredParams[i])
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], scalarMont, Qi, mredParams)
		}
	}
}

// MulScalarBigintLvl multiplies each coefficientsof p1 by a big.Int scalar and applies a modular reduction, returning the result on p2.
// To be used when the scalar is bigger than 64 bits.
func (context *Context) MulScalarBigintLvl(level uint64, p1 *Poly, scalar *big.Int, p2 *Poly) {

	scalarQi := new(big.Int)
	var Qi, scalarMont uint64

	for i := uint64(0); i < level+1; i++ {
		Qi = context.Modulus[i]
		scalarQi.Mod(scalar, NewUint(Qi))
		scalarMont = MForm(BRedAdd(scalarQi.Uint64(), Qi, context.bredParams[i]), Qi, context.bredParams[i])
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], scalarMont, Qi, mredParams)
		}
	}
}

// Shift circulary shifts the coefficients of the polynomial p1 by n to the left and returns the result on the receiver polynomial.
func (context *Context) Shift(p1 *Poly, n uint64, p2 *Poly) {
	mask := uint64((1 << context.N) - 1)
	for i := range context.Modulus {
		p2.Coeffs[i] = append(p1.Coeffs[i][(n&mask):], p1.Coeffs[i][:(n&mask)]...)
	}
}

// MForm setss p1 in conventional form to its Montgomeryform, returning the result on p2.
func (context *Context) MForm(p1, p2 *Poly) {

	for i, qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := context.bredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = MForm(p1tmp[j], qi, bredParams)
		}
	}
}

// MFormLvl setss p1 in conventional form to its Montgomeryform, returning the result on p2.
func (context *Context) MFormLvl(level uint64, p1, p2 *Poly) {

	var qi uint64
	var bredParams []uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		bredParams = context.bredParams[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = MForm(p1tmp[j], qi, bredParams)
		}
	}
}

// InvMForm setss p1 in Montgomeryform to its conventional form, returning the result on p2.
func (context *Context) InvMForm(p1, p2 *Poly) {

	for i, qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = InvMForm(p1tmp[j], qi, mredParams)
		}
	}
}

// PermuteNTTIndex computes the index table for PermuteNTT.
func PermuteNTTIndex(gen, N uint64) (index []uint64) {

	var mask, logN, tmp1, tmp2 uint64

	logN = uint64(bits.Len64(N) - 1)

	mask = (N << 1) - 1

	index = make([]uint64, N)

	for i := uint64(0); i < N; i++ {
		tmp1 = 2*utils.BitReverse64(i, logN) + 1

		tmp2 = ((gen * tmp1 & mask) - 1) >> 1

		index[i] = utils.BitReverse64(tmp2, logN)
	}

	return
}

// PermuteNTT applies the galois transform on a polynomial in the NTT domain.
// It maps the coefficients x^i to x^(gen*i)
// Careful, not inplace!
func PermuteNTT(polIn *Poly, gen uint64, polOut *Poly) {

	var N, tmp uint64

	N = uint64(len(polIn.Coeffs[0]))

	index := PermuteNTTIndex(gen, N)

	for j := uint64(0); j < N; j++ {

		tmp = index[j]

		for i := 0; i < len(polIn.Coeffs); i++ {

			polOut.Coeffs[i][j] = polIn.Coeffs[i][tmp]
		}
	}
}

// PermuteNTT applies the galois transform on a polynomial in the NTT domain.
// It maps the coefficients x^i to x^(gen*i) using the PermuteNTTIndex table.
// Careful, not inplace!
func PermuteNTTWithIndex(polIn *Poly, index []uint64, polOut *Poly) {

	var tmp uint64

	for j := uint64(0); j < uint64(len(polIn.Coeffs[0])); j++ {

		tmp = index[j]

		for i := 0; i < len(polIn.Coeffs); i++ {
			polOut.Coeffs[i][j] = polIn.Coeffs[i][tmp]
		}
	}
}

// Permute applies the galois transform on a polynonial outside of the NTT domain.
// It maps the coefficients x^i to x^(gen*i)
// Careful, not inplace!
func (context *Context) Permute(polIn *Poly, gen uint64, polOut *Poly) {

	var mask, index, indexRaw, logN, tmp uint64

	mask = context.N - 1

	logN = uint64(bits.Len64(mask))

	for i := uint64(0); i < context.N; i++ {

		indexRaw = i * gen

		index = indexRaw & mask

		tmp = (indexRaw >> logN) & 1

		for j, qi := range context.Modulus {

			polOut.Coeffs[j][index] = polIn.Coeffs[j][i]*(tmp^1) | (qi-polIn.Coeffs[j][i])*tmp
		}
	}
}

// MulByPow2New multiplies the input polynomial by 2^pow2 and returns the result on a new polynomial.
func (context *Context) MulByPow2New(p1 *Poly, pow2 uint64) (p2 *Poly) {
	p2 = context.NewPoly()
	context.MulByPow2(p1, pow2, p2)
	return
}

// MulByPow2 multiplies the input polynomial by 2^pow2 and returns the result on the receiver polynomial.
func (context *Context) MulByPow2(p1 *Poly, pow2 uint64, p2 *Poly) {
	context.MForm(p1, p2)
	for i, Qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = PowerOf2(p1tmp[j], pow2, Qi, mredParams)
		}
	}
}

// MulByPow2Lvl multiplies the input polynomial by 2^pow2 and returns the result on the receiver polynomial.
func (context *Context) MulByPow2Lvl(level uint64, p1 *Poly, pow2 uint64, p2 *Poly) {
	context.MFormLvl(level, p1, p2)
	var qi uint64
	var mredParams uint64
	for i := uint64(0); i < level+1; i++ {
		qi = context.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams = context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = PowerOf2(p1tmp[j], pow2, qi, mredParams)
		}
	}
}

// MultByMonomialNew multiplies the input polynomial by x^monomialDeg and returns the result on a new polynomial.
func (context *Context) MultByMonomialNew(p1 *Poly, monomialDeg uint64) (p2 *Poly) {
	p2 = context.NewPoly()
	context.MultByMonomial(p1, monomialDeg, p2)
	return
}

// MultByMonomial multiplies the input polynomial by x^monomialDeg and returns the result on the receiver polynomial.
func (context *Context) MultByMonomial(p1 *Poly, monomialDeg uint64, p2 *Poly) {

	var shift uint64

	shift = monomialDeg % (context.N << 1)

	if shift == 0 {

		for i := range context.Modulus {
			p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
			for j := uint64(0); j < context.N; j++ {

				p2tmp[j] = p1tmp[j]

			}
		}

	} else {

		tmpx := context.NewPoly()

		if shift < context.N {

			for i := range context.Modulus {
				p1tmp, tmpxT := p1.Coeffs[i], tmpx.Coeffs[i]
				for j := uint64(0); j < context.N; j++ {

					tmpxT[j] = p1tmp[j]

				}
			}

		} else {

			for i, qi := range context.Modulus {
				p1tmp, tmpxT := p1.Coeffs[i], tmpx.Coeffs[i]
				for j := uint64(0); j < context.N; j++ {
					tmpxT[j] = qi - p1tmp[j]

				}
			}
		}

		shift %= context.N

		for i, qi := range context.Modulus {
			p2tmp, tmpxT := p2.Coeffs[i], tmpx.Coeffs[i]
			for j := uint64(0); j < shift; j++ {
				p2tmp[j] = qi - tmpxT[context.N-shift+j]
			}
		}

		for i := range context.Modulus {
			p2tmp, tmpxT := p2.Coeffs[i], tmpx.Coeffs[i]
			for j := shift; j < context.N; j++ {
				p2tmp[j] = tmpxT[j-shift]

			}
		}
	}
}

// MulByVectorMontgomery multiplies p1 by a vector of uint64 coefficients and returns the result on p2.
func (context *Context) MulByVectorMontgomery(p1 *Poly, vector []uint64, p2 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], vector[j], qi, mredParams)
		}
	}
}

// MulByVectorMontgomeryAndAddNoMod multiplies p1 by a vector of uint64 coefficients and adds the result on p2 without modular reduction.
func (context *Context) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint64, p2 *Poly) {
	for i, qi := range context.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := context.mredParams[i]
		for j := uint64(0); j < context.N; j++ {
			p2tmp[j] += MRed(p1tmp[j], vector[j], qi, mredParams)
		}
	}
}

// BitReverse applies a bit reverse permutation on the coefficients of the input polynomial and returns the result on the receiver polynomial.
// Can safely be used for inplace permutation.
func (context *Context) BitReverse(p1, p2 *Poly) {
	bitLenOfN := uint64(bits.Len64(context.N) - 1)

	if p1 != p2 {
		for i := range context.Modulus {
			p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
			for j := uint64(0); j < context.N; j++ {
				p2tmp[utils.BitReverse64(j, bitLenOfN)] = p1tmp[j]
			}
		}
	} else { // In place in case p1 = p2
		for x := range context.Modulus {
			p2tmp := p2.Coeffs[x]
			for i := uint64(0); i < context.N; i++ {
				j := utils.BitReverse64(i, bitLenOfN)
				if i < j {
					p2tmp[i], p2tmp[j] = p2tmp[j], p2tmp[i]
				}
			}
		}
	}
}

// Rotate applies a Galoi Automorphism on p1 in NTT form,
// rotating the coefficients to the right by n, returning the result on p2.
// Requires the data to permuted in bitreversal order before applying NTT.
func (context *Context) Rotate(p1 *Poly, n uint64, p2 *Poly) {

	var root, gal uint64

	n &= (1 << context.N) - 1

	for i, qi := range context.Modulus {

		mredParams := context.mredParams[i]

		root = MRed(context.psiMont[i], context.psiMont[i], qi, mredParams)

		root = modexpMontgomery(root, n, qi, mredParams, context.bredParams[i])

		gal = MForm(1, qi, context.bredParams[i])

		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]

		for j := uint64(1); j < context.N; j++ {

			gal = MRed(gal, root, qi, mredParams)

			p2tmp[j] = MRed(p1tmp[j], gal, qi, mredParams)

		}
	}
}
