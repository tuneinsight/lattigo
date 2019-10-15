package ring

import (
	"math/bits"
)

// Add adds p1 to p2 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Add(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed(p1.Coeffs[i][j]+p2.Coeffs[i][j], qi)
		}
	}
}

// AddNoMod adds p1 to p2 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*Qi -1].
func (context *Context) AddNoMod(p1, p2, p3 *Poly) {
	for i := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = p1.Coeffs[i][j] + p2.Coeffs[i][j]
		}
	}
}

// Sub subtracts p2 to p1 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Sub(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed((p1.Coeffs[i][j]+qi)-p2.Coeffs[i][j], qi)
		}
	}
}

// SubNoMod subtracts p2 to p1 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*Qi -1].
func (context *Context) SubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = (p1.Coeffs[i][j] + qi) - p2.Coeffs[i][j]
		}
	}
}

// Neg sets all coefficients of p1 to their additive inverse, returning the result on p2.
func (context *Context) Neg(p1, p2 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = qi - p1.Coeffs[i][j]
		}
	}
}

// Reduce applies a modular reduction over the coefficients of p1 returning the result on p2.
func (context *Context) Reduce(p1, p2 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = BRedAdd(p1.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// Mod applies a modular reduction by m over the coefficients of p1, returning the result on p2.
func (context *Context) Mod(p1 *Poly, m uint64, p2 *Poly) {
	params := BRedParams(m)
	for i := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = BRedAdd(p1.Coeffs[i][j], m, params)
		}
	}
}

// AND applies a logical AND of m to the coefficients of p1,  returning the result on p2.
func (context *Context) AND(p1 *Poly, m uint64, p2 *Poly) {
	for i := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = p1.Coeffs[i][j] & m
		}
	}
}

// OR applies a logical OR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) OR(p1 *Poly, m uint64, p2 *Poly) {
	for i := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = p1.Coeffs[i][j] | m
		}
	}
}

// XOR applies a logical XOR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) XOR(p1 *Poly, m uint64, p2 *Poly) {
	for i := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = p1.Coeffs[i][j] ^ m
		}
	}
}

// MulCoeffs multiplies p1 by p2 coefficient wise with a Barrett modular reduction, returning the result on p3.
func (context *Context) MulCoeffs(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = BRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient wise with a Barret modular reduction, adding the result to p3 with modular reduction.
func (context *Context) MulCoeffsAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed(p3.Coeffs[i][j]+BRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.bredParams[i]), qi)
		}
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient wise with a Barrett modular reduction, adding the result to p3 without modular reduction.
func (context *Context) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] += BRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, returning the result on p3.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomery(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i])
		}
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, adding the result to p3.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed(p3.Coeffs[i][j]+MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i]), qi)
		}
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, adding the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] += MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i])
		}
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, subtracting the result to p3 with modular reduction.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed(p3.Coeffs[i][j]+(qi-MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i])), qi)
		}
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient wise with a Montgomery modular reduction, subtracting the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in Montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = p3.Coeffs[i][j] + (qi - MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i]))
		}
	}
}

// MulCoeffsConstant multiplies p1 by p2 coefficient wise with a constant time Barrett modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*Qi -1].
func (context *Context) MulCoeffsConstant(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = BRedConstant(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// MulCoeffsConstantMontgomery multiplies p1 by p2 coefficient wise with a constant time Montgomery modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*Qi -1].
func (context *Context) MulCoeffsConstantMontgomery(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p3.Coeffs[i][j] = MRedConstant(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i])
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
// Expects either p1 or p2 to be in Montgomery form for correctness.
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

		for i := uint64(0); i < context.N; i++ {

			for j := uint64(0); j < i; j++ {
				p3.Coeffs[x][j] = CRed(p3.Coeffs[x][j]+(qi-MRed(p1Copy.Coeffs[x][i], p2Copy.Coeffs[x][context.N-i+j], qi, context.mredParams[x])), qi)
			}

			for j := uint64(i); j < context.N; j++ {
				p3.Coeffs[x][j] = CRed(p3.Coeffs[x][j]+MRed(p1Copy.Coeffs[x][i], p2Copy.Coeffs[x][j-i], qi, context.mredParams[x]), qi)
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

		for i := uint64(0); i < context.N; i++ {

			for j := uint64(0); j < i; j++ {
				p3.Coeffs[x][j] = CRed(p3.Coeffs[x][j]+(qi-MRed(p1Copy.Coeffs[x][i], p2Copy.Coeffs[x][context.N-i+j], qi, context.mredParams[x])), qi)
			}

			for j := uint64(i); j < context.N; j++ {
				p3.Coeffs[x][j] = CRed(p3.Coeffs[x][j]+MRed(p1Copy.Coeffs[x][i], p2Copy.Coeffs[x][j-i], qi, context.mredParams[x]), qi)
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
		for x := uint64(0); x < context.N; x++ {
			p2.Coeffs[i][x] = 1
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

// AddScalar adds to each coefficient of p1 a scalar and applies a modular reduction, returning the result on p2.
func (context *Context) AddScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = CRed(p1.Coeffs[i][j]+scalar, Qi)
		}
	}
}

// SubScalar subtracts to each coefficient of p1 a scalar and applies a modular reduction, returning the result on p2.
func (context *Context) SubScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = CRed(p1.Coeffs[i][j]+(Qi-scalar), Qi)
		}
	}
}

// MulScalar multiplies each coefficient of p1 by a scalar and applies a modular reduction, returning the result on p2.
func (context *Context) MulScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	var scalarMont uint64
	for i, Qi := range context.Modulus {
		scalarMont = MForm(BRedAdd(scalar, Qi, context.bredParams[i]), Qi, context.bredParams[i])
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = MRed(p1.Coeffs[i][j], scalarMont, Qi, context.mredParams[i])
		}
	}
}

// MulScalarBigint multiplies each coefficient of p1 by an Int scalar and applies a modular reduction, returning the result on p2.
// To be used when the scalar is bigger than 64 bits.
func (context *Context) MulScalarBigint(p1 *Poly, scalar *Int, p2 *Poly) {

	var QiB Int
	var coeff Int

	for i, Qi := range context.Modulus {
		QiB.SetInt(int64(Qi))
		for j := uint64(0); j < context.N; j++ {
			coeff.SetUint(p1.Coeffs[i][j])
			coeff.Mul(&coeff, scalar)
			coeff.Mod(&coeff, &QiB)
			p2.Coeffs[i][j] = coeff.Uint64()
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

// MForm sets p1 in conventional form to its Montgomery form, returning the result on p2.
func (context *Context) MForm(p1, p2 *Poly) {

	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = MForm(p1.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// InvMForm sets p1 in Montgomery form to its conventional form, returning the result on p2.
func (context *Context) InvMForm(p1, p2 *Poly) {

	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = InvMForm(p1.Coeffs[i][j], qi, context.mredParams[i])
		}
	}
}

// PermuteNTT applies the galois transform on a polynomial in the NTT domain.
// It maps the coefficients x^i to x^(gen*i)
// Careful, not inplace!
func PermuteNTT(polIn *Poly, gen uint64, polOut *Poly) {

	var N, mask, logN, tmp, index uint64

	N = uint64(len(polIn.Coeffs[0]))

	logN = uint64(bits.Len64(N) - 1)

	mask = (N << 1) - 1

	for j := uint64(0); j < N; j++ {

		index = 2*bitReverse64(j, logN) + 1

		tmp = ((gen * index & mask) - 1) >> 1

		index = bitReverse64(tmp, logN)

		for i := 0; i < len(polIn.Coeffs); i++ {

			polOut.Coeffs[i][j] = polIn.Coeffs[i][index]
		}
	}
}

// Permute applies the galois transform on a polynonial outside of the NTT domain.
// It maps the coefficients x^i to x^(gen*i).
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
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = PowerOf2(p2.Coeffs[i][j], pow2, Qi, context.mredParams[i])
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
			for j := uint64(0); j < context.N; j++ {

				p2.Coeffs[i][j] = p1.Coeffs[i][j]

			}
		}

	} else {

		tmpx := context.NewPoly()

		if shift < context.N {

			for i := range context.Modulus {
				for j := uint64(0); j < context.N; j++ {

					tmpx.Coeffs[i][j] = p1.Coeffs[i][j]

				}
			}

		} else {

			for i, qi := range context.Modulus {
				for j := uint64(0); j < context.N; j++ {

					tmpx.Coeffs[i][j] = qi - p1.Coeffs[i][j]

				}
			}
		}

		shift %= context.N

		for i, qi := range context.Modulus {
			for j := uint64(0); j < shift; j++ {

				p2.Coeffs[i][j] = qi - tmpx.Coeffs[i][context.N-shift+j]

			}
		}

		for i := range context.Modulus {
			for j := shift; j < context.N; j++ {

				p2.Coeffs[i][j] = tmpx.Coeffs[i][j-shift]

			}
		}
	}
}

// MulByVectorMontgomery multiplies p1 by a vector of uint64 coefficients and returns the result on p2.
func (context *Context) MulByVectorMontgomery(p1 *Poly, vector []uint64, p2 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] = MRed(p1.Coeffs[i][j], vector[j], qi, context.mredParams[i])
		}
	}
}

// MulByVectorMontgomeryAndAddNoMod multiplies p1 by a vector of uint64 coefficients and adds the result on p2 without modular reduction.
func (context *Context) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint64, p2 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint64(0); j < context.N; j++ {
			p2.Coeffs[i][j] += MRed(p1.Coeffs[i][j], vector[j], qi, context.mredParams[i])
		}
	}
}

// BitReverse applies a bit reverse permutation on the coefficients of the input polynomial and returns the result on the receiver polynomial.
// Can safely be used for inplace permutation.
func (context *Context) BitReverse(p1, p2 *Poly) {
	bitLenOfN := uint64(bits.Len64(context.N) - 1)

	if p1 != p2 {
		for i := range context.Modulus {
			for j := uint64(0); j < context.N; j++ {
				p2.Coeffs[i][bitReverse64(j, bitLenOfN)] = p1.Coeffs[i][j]
			}
		}
	} else { // In place in case p1 = p2
		for x := range context.Modulus {
			for i := uint64(0); i < context.N; i++ {
				j := bitReverse64(i, bitLenOfN)
				if i < j {
					p2.Coeffs[x][i], p2.Coeffs[x][j] = p2.Coeffs[x][j], p2.Coeffs[x][i]
				}
			}
		}
	}
}

// Rotate applies a Galois Automorphism on p1 in NTT form,
// rotating the coefficients to the right by n, returning the result on p2.
// Requires the data to permuted in bitreversal order before applying NTT.
func (context *Context) Rotate(p1 *Poly, n uint64, p2 *Poly) {

	var root, gal uint64

	n &= (1 << context.N) - 1

	for i, qi := range context.Modulus {

		root = MRed(context.psiMont[i], context.psiMont[i], qi, context.mredParams[i])

		root = modexpMontgomery(root, n, qi, context.mredParams[i], context.bredParams[i])

		gal = MForm(1, qi, context.bredParams[i])

		for j := uint64(1); j < context.N; j++ {

			gal = MRed(gal, root, qi, context.mredParams[i])

			p2.Coeffs[i][j] = MRed(p1.Coeffs[i][j], gal, qi, context.mredParams[i])

		}
	}
}
