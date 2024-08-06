package ring

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Add evaluates p3 = p1 + p2 coefficient-wise in the ring.
func (r Ring) Add(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.Add(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// AddLazy evaluates p3 = p1 + p2 coefficient-wise in the ring, with p3 in [0, 2*modulus-1].
func (r Ring) AddLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.AddLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// Sub evaluates p3 = p1 - p2 coefficient-wise in the ring.
func (r Ring) Sub(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.Sub(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// SubLazy evaluates p3 = p1 - p2 coefficient-wise in the ring, with p3 in [0, 2*modulus-1].
func (r Ring) SubLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.SubLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// Neg evaluates p2 = -p1 coefficient-wise in the ring.
func (r Ring) Neg(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.Neg(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// Reduce evaluates p2 = p1 coefficient-wise mod modulus in the ring.
func (r Ring) Reduce(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.Reduce(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// ReduceLazy evaluates p2 = p1 coefficient-wise mod modulus in the ring, with p2 in [0, 2*modulus-1].
func (r Ring) ReduceLazy(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.ReduceLazy(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// MulCoeffsBarrett evaluates p3 = p1 * p2 coefficient-wise in the ring, with Barrett reduction.
func (r Ring) MulCoeffsBarrett(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsBarrett(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsBarrettLazy evaluates p3 = p1 * p2 coefficient-wise in the ring, with Barrett reduction, with p3 in [0, 2*modulus-1].
func (r Ring) MulCoeffsBarrettLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsBarrettLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsBarrettThenAdd evaluates p3 = p3 + p1 * p2 coefficient-wise in the ring, with Barrett reduction.
func (r Ring) MulCoeffsBarrettThenAdd(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsBarrettThenAdd(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsBarrettThenAddLazy evaluates p3 = p1 * p2 coefficient-wise in the ring, with Barrett reduction, with p3 in [0, 2*modulus-1].
func (r Ring) MulCoeffsBarrettThenAddLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsBarrettThenAddLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomery evaluates p3 = p1 * p2 coefficient-wise in the ring, with Montgomery reduction.
func (r Ring) MulCoeffsMontgomery(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomery(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomeryLazy evaluates p3 = p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
func (r Ring) MulCoeffsMontgomeryLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomeryLazyThenNeg evaluates p3 = -p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
func (r Ring) MulCoeffsMontgomeryLazyThenNeg(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryLazyThenNeg(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomeryThenAdd evaluates p3 = p3 + p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
func (r Ring) MulCoeffsMontgomeryThenAdd(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryThenAdd(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomeryThenAddLazy evaluates p3 = p3 + p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
func (r Ring) MulCoeffsMontgomeryThenAddLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryThenAddLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomeryLazyThenAddLazy evaluates p3 = p3 + p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 3*modulus-2].
func (r Ring) MulCoeffsMontgomeryLazyThenAddLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryLazyThenAddLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomeryThenSub evaluates p3 = p3 - p1 * p2 coefficient-wise in the ring, with Montgomery reduction.
func (r Ring) MulCoeffsMontgomeryThenSub(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryThenSub(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomeryThenSubLazy evaluates p3 = p3 - p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 2*modulus-1].
func (r Ring) MulCoeffsMontgomeryThenSubLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryThenSubLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// MulCoeffsMontgomeryLazyThenSubLazy evaluates p3 = p3 - p1 * p2 coefficient-wise in the ring, with Montgomery reduction, with p3 in [0, 3*modulus-2].
func (r Ring) MulCoeffsMontgomeryLazyThenSubLazy(p1, p2, p3 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryLazyThenSubLazy(p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i])
	}
}

// AddScalar evaluates p2 = p1 + scalar coefficient-wise in the ring.
func (r Ring) AddScalar(p1 Poly, scalar uint64, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.AddScalar(p1.Coeffs[i], scalar, p2.Coeffs[i])
	}
}

// AddScalarBigint evaluates p2 = p1 + scalar coefficient-wise in the ring.
func (r Ring) AddScalarBigint(p1 Poly, scalar *big.Int, p2 Poly) {
	tmp := new(big.Int)
	for i, s := range r.SubRings[:r.level+1] {
		s.AddScalar(p1.Coeffs[i], tmp.Mod(scalar, bignum.NewInt(s.Modulus)).Uint64(), p2.Coeffs[i])
	}
}

// AddDoubleRNSScalar evaluates p2 = p1[:N/2] + scalar0 || p1[N/2] + scalar1 coefficient-wise in the ring,
// with the scalar values expressed in the CRT decomposition at a given level.
func (r Ring) AddDoubleRNSScalar(p1 Poly, scalar0, scalar1 RNSScalar, p2 Poly) {
	NHalf := r.N() >> 1
	for i, s := range r.SubRings[:r.level+1] {
		s.AddScalar(p1.Coeffs[i][:NHalf], scalar0[i], p2.Coeffs[i][:NHalf])
		s.AddScalar(p1.Coeffs[i][NHalf:], scalar1[i], p2.Coeffs[i][NHalf:])
	}
}

// SubDoubleRNSScalar evaluates p2 = p1[:N/2] - scalar0 || p1[N/2] - scalar1 coefficient-wise in the ring,
// with the scalar values expressed in the CRT decomposition at a given level.
func (r Ring) SubDoubleRNSScalar(p1 Poly, scalar0, scalar1 RNSScalar, p2 Poly) {
	NHalf := r.N() >> 1
	for i, s := range r.SubRings[:r.level+1] {
		s.SubScalar(p1.Coeffs[i][:NHalf], scalar0[i], p2.Coeffs[i][:NHalf])
		s.SubScalar(p1.Coeffs[i][NHalf:], scalar1[i], p2.Coeffs[i][NHalf:])
	}
}

// SubScalar evaluates p2 = p1 - scalar coefficient-wise in the ring.
func (r Ring) SubScalar(p1 Poly, scalar uint64, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.SubScalar(p1.Coeffs[i], scalar, p2.Coeffs[i])
	}
}

// SubScalarBigint evaluates p2 = p1 - scalar coefficient-wise in the ring.
func (r Ring) SubScalarBigint(p1 Poly, scalar *big.Int, p2 Poly) {
	tmp := new(big.Int)
	for i, s := range r.SubRings[:r.level+1] {
		s.SubScalar(p1.Coeffs[i], tmp.Mod(scalar, bignum.NewInt(s.Modulus)).Uint64(), p2.Coeffs[i])
	}
}

// MulScalar evaluates p2 = p1 * scalar coefficient-wise in the ring.
func (r Ring) MulScalar(p1 Poly, scalar uint64, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulScalarMontgomery(p1.Coeffs[i], MForm(scalar, s.Modulus, s.BRedConstant), p2.Coeffs[i])
	}
}

// MulScalarThenAdd evaluates p2 = p2 + p1 * scalar coefficient-wise in the ring.
func (r Ring) MulScalarThenAdd(p1 Poly, scalar uint64, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulScalarMontgomeryThenAdd(p1.Coeffs[i], MForm(scalar, s.Modulus, s.BRedConstant), p2.Coeffs[i])
	}
}

// MulRNSScalarMontgomery evaluates p2 = p1 * scalar coefficient-wise in the ring, with a scalar value expressed in the CRT decomposition at a given level.
// It assumes the scalar decomposition to be in Montgomery form.
func (r Ring) MulRNSScalarMontgomery(p1 Poly, scalar RNSScalar, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulScalarMontgomery(p1.Coeffs[i], scalar[i], p2.Coeffs[i])
	}
}

// MulScalarThenSub evaluates p2 = p2 - p1 * scalar coefficient-wise in the ring.
func (r Ring) MulScalarThenSub(p1 Poly, scalar uint64, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		scalarNeg := MForm(s.Modulus-BRedAdd(scalar, s.Modulus, s.BRedConstant), s.Modulus, s.BRedConstant)
		s.MulScalarMontgomeryThenAdd(p1.Coeffs[i], scalarNeg, p2.Coeffs[i])
	}
}

// MulScalarBigint evaluates p2 = p1 * scalar coefficient-wise in the ring.
func (r Ring) MulScalarBigint(p1 Poly, scalar *big.Int, p2 Poly) {
	scalarQi := new(big.Int)
	for i, s := range r.SubRings[:r.level+1] {
		scalarQi.Mod(scalar, bignum.NewInt(s.Modulus))
		s.MulScalarMontgomery(p1.Coeffs[i], MForm(scalarQi.Uint64(), s.Modulus, s.BRedConstant), p2.Coeffs[i])
	}
}

// MulScalarBigintThenAdd evaluates p2 = p1 * scalar coefficient-wise in the ring.
func (r Ring) MulScalarBigintThenAdd(p1 Poly, scalar *big.Int, p2 Poly) {
	scalarQi := new(big.Int)
	for i, s := range r.SubRings[:r.level+1] {
		scalarQi.Mod(scalar, bignum.NewInt(s.Modulus))
		s.MulScalarMontgomeryThenAdd(p1.Coeffs[i], MForm(scalarQi.Uint64(), s.Modulus, s.BRedConstant), p2.Coeffs[i])
	}
}

// MulDoubleRNSScalar evaluates p2 = p1[:N/2] * scalar0 || p1[N/2] * scalar1 coefficient-wise in the ring,
// with the scalar values expressed in the CRT decomposition at a given level.
func (r Ring) MulDoubleRNSScalar(p1 Poly, scalar0, scalar1 RNSScalar, p2 Poly) {
	NHalf := r.N() >> 1
	for i, s := range r.SubRings[:r.level+1] {
		s.MulScalarMontgomery(p1.Coeffs[i][:NHalf], MForm(scalar0[i], s.Modulus, s.BRedConstant), p2.Coeffs[i][:NHalf])
		s.MulScalarMontgomery(p1.Coeffs[i][NHalf:], MForm(scalar1[i], s.Modulus, s.BRedConstant), p2.Coeffs[i][NHalf:])
	}
}

// MulDoubleRNSScalarThenAdd evaluates p2 = p2 + p1[:N/2] * scalar0 || p1[N/2] * scalar1 coefficient-wise in the ring,
// with the scalar values expressed in the CRT decomposition at a given level.
func (r Ring) MulDoubleRNSScalarThenAdd(p1 Poly, scalar0, scalar1 RNSScalar, p2 Poly) {
	NHalf := r.N() >> 1
	for i, s := range r.SubRings[:r.level+1] {
		s.MulScalarMontgomeryThenAdd(p1.Coeffs[i][:NHalf], MForm(scalar0[i], s.Modulus, s.BRedConstant), p2.Coeffs[i][:NHalf])
		s.MulScalarMontgomeryThenAdd(p1.Coeffs[i][NHalf:], MForm(scalar1[i], s.Modulus, s.BRedConstant), p2.Coeffs[i][NHalf:])
	}
}

// EvalPolyScalar evaluate p2 = p1(scalar) coefficient-wise in the ring.
func (r Ring) EvalPolyScalar(p1 []Poly, scalar uint64, p2 Poly) {
	p2.Copy(p1[len(p1)-1])
	for i := len(p1) - 1; i > 0; i-- {
		r.MulScalar(p2, scalar, p2)
		r.Add(p2, p1[i-1], p2)
	}
}

// Shift evaluates p2 = p2<<<k coefficient-wise in the ring.
func (r Ring) Shift(p1 Poly, k int, p2 Poly) {
	for i := range p1.Coeffs {
		utils.RotateSliceAllocFree(p1.Coeffs[i], k, p2.Coeffs[i])
	}
}

// MForm evaluates p2 = p1 * (2^64)^-1 coefficient-wise in the ring.
func (r Ring) MForm(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MForm(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// MFormLazy evaluates p2 = p1 * (2^64)^-1 coefficient-wise in the ring with p2 in [0, 2*modulus-1].
func (r Ring) MFormLazy(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MFormLazy(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// IMForm evaluates p2 = p1 * 2^64 coefficient-wise in the ring.
func (r Ring) IMForm(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.IMForm(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// MultByMonomial evaluates p2 = p1 * X^k coefficient-wise in the ring.
func (r Ring) MultByMonomial(p1 Poly, k int, p2 Poly) {

	N := r.N()

	shift := (k + (N << 1)) % (N << 1)

	if shift == 0 {

		for i := range r.SubRings[:r.level+1] {
			p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
			for j := 0; j < N; j++ {
				p2tmp[j] = p1tmp[j]
			}
		}

	} else {

		tmpx := r.NewPoly()

		if shift < N {

			for i := range r.SubRings[:r.level+1] {
				p1tmp, tmpxT := p1.Coeffs[i], tmpx.Coeffs[i]
				for j := 0; j < N; j++ {
					tmpxT[j] = p1tmp[j]
				}
			}

		} else {

			for i, s := range r.SubRings[:r.level+1] {
				qi := s.Modulus
				p1tmp, tmpxT := p1.Coeffs[i], tmpx.Coeffs[i]
				for j := 0; j < N; j++ {
					tmpxT[j] = qi - p1tmp[j]
				}
			}
		}

		shift %= N

		for i, s := range r.SubRings[:r.level+1] {
			qi := s.Modulus
			p2tmp, tmpxT := p2.Coeffs[i], tmpx.Coeffs[i]
			for j := 0; j < shift; j++ {
				p2tmp[j] = qi - tmpxT[N-shift+j]
			}
		}

		for i := range r.SubRings[:r.level+1] {
			p2tmp, tmpxT := p2.Coeffs[i], tmpx.Coeffs[i]
			for j := shift; j < N; j++ {
				p2tmp[j] = tmpxT[j-shift]

			}
		}
	}
}

// MulByVectorMontgomery evaluates p2 = p1 * vector coefficient-wise in the ring.
func (r Ring) MulByVectorMontgomery(p1 Poly, vector []uint64, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomery(p1.Coeffs[i], vector, p2.Coeffs[i])
	}
}

// MulByVectorMontgomeryThenAddLazy evaluates p2 = p2 + p1 * vector coefficient-wise in the ring.
func (r Ring) MulByVectorMontgomeryThenAddLazy(p1 Poly, vector []uint64, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.MulCoeffsMontgomeryThenAddLazy(p1.Coeffs[i], vector, p2.Coeffs[i])
	}
}

// MapSmallDimensionToLargerDimensionNTT maps Y = X^{N/n} -> X directly in the NTT domain
func MapSmallDimensionToLargerDimensionNTT(polSmall, polLarge Poly) {
	gap := len(polLarge.Coeffs[0]) / len(polSmall.Coeffs[0])
	for j := range polSmall.Coeffs {
		tmp0 := polSmall.Coeffs[j]
		tmp1 := polLarge.Coeffs[j]
		for i := range polSmall.Coeffs[0] {
			coeff := tmp0[i]
			for w := 0; w < gap; w++ {
				tmp1[i*gap+w] = coeff
			}
		}
	}
}
