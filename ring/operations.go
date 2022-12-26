package ring

import (
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// Add adds p1 to p2 coefficient-wise and writes the result on p3.
func (r *Ring) Add(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		AddVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus)
	}
}

// AddNoMod adds p1 to p2 coefficient-wise without modular reduction and writes the result on p3.
func (r *Ring) AddNoMod(p1, p2, p3 *Poly) {
	N := r.N()
	for i := range r.Tables[:r.level+1] {
		AddVecNoMod(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N])
	}
}

// Sub subtracts p2 to p1 coefficient-wise and writes the result on p3.
func (r *Ring) Sub(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		SubVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus)
	}
}

// SubNoMod subtracts p2 to p1 coefficient-wise without modular reduction and writes the result on p3.
func (r *Ring) SubNoMod(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		SubVecNomod(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus)
	}
}

// Neg negates the coefficients of p1 and and writes the result on p2.
func (r *Ring) Neg(p1, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		NegVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], table.Modulus)
	}
}

// Reduce applies a modular reduction on the coefficients of p1 and writes the result on p2.
func (r *Ring) Reduce(p1, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		ReduceVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], table.Modulus, table.BRedParams)
	}
}

// ReduceConstant applies a modular reduction on the coefficients of p1 and writes the result on p2.
// Returns values in [0, 2q-1]
func (r *Ring) ReduceConstant(p1, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		ReduceConstantVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], table.Modulus, table.BRedParams)
	}
}

// Mod applies a modular reduction by m on the coefficients of p1 and writes the result on p2.
func (r *Ring) Mod(p1 *Poly, m uint64, p2 *Poly) {
	bredParams := BRedParams(m)
	N := r.N()
	for i := range r.Tables[:r.level+1] {
		ModVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], m, bredParams)
	}
}

// MulCoeffs multiplies p1 by p2 coefficient-wise with a
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffs(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.BRedParams)
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient-wise with
// a Barret modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsAndAdd(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsAndAddVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.BRedParams)
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient-wise with a Barrett
// modular reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsAndAddNoModVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.BRedParams)
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsMontgomery(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsMontgomeryConstant multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsMontgomeryConstant(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryConstantVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsMontgomeryConstantAndNeg multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and writes the negative of the result on p3.
func (r *Ring) MulCoeffsMontgomeryConstantAndNeg(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryConstantAndNeg(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryAndAddVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient-wise with a Montgomery modular
// reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryAndAddNoModVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsMontgomeryConstantAndAddNoMod multiplies p1 by p2 coefficient-wise with a constant-time Montgomery
// modular reduction and adds the result to p3 without modular reduction.
// Return values in [0, 3q-1]
func (r *Ring) MulCoeffsMontgomeryConstantAndAddNoMod(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryConstantAndAddNoModVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r *Ring) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryAndSubVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and subtracts the result from p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryAndSubNoMod(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsMontgomeryConstantAndSubNoMod multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and subtracts the result from p3 without modular reduction.
// Return values in [0, 3q-1]
func (r *Ring) MulCoeffsMontgomeryConstantAndSubNoMod(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryConstantAndSubNoMod(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulCoeffsConstant multiplies p1 by p2 coefficient-wise with a constant-time
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsConstant(p1, p2, p3 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsConstantVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], p3.Coeffs[i][:N], table.Modulus, table.BRedParams)
	}
}

// AddScalar adds a scalar to each coefficient of p1 and writes the result on p2.
func (r *Ring) AddScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		AddScalarVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], scalar, table.Modulus)
	}
}

// AddScalarBigint adds a big.Int scalar to each coefficient of p1 and writes the result on p2.
func (r *Ring) AddScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	N := r.N()
	tmp := new(big.Int)
	for i, table := range r.Tables[:r.level+1] {
		AddScalarVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], tmp.Mod(scalar, NewUint(table.Modulus)).Uint64(), table.Modulus)
	}
}

// SubScalar subtracts a scalar from each coefficient of p1 and writes the result on p2.
func (r *Ring) SubScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		SubScalarVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], scalar, table.Modulus)
	}
}

// SubScalarBigint subtracts a big.Int scalar from each coefficient of p1 and writes the result on p2.
func (r *Ring) SubScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	N := r.N()
	tmp := new(big.Int)
	for i, table := range r.Tables[:r.level+1] {
		SubScalarVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], tmp.Mod(scalar, NewUint(table.Modulus)).Uint64(), table.Modulus)
	}
}

// MulScalar multiplies each coefficient of p1 by a scalar and writes the result on p2.
func (r *Ring) MulScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulScalarMontgomeryVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], MForm(BRedAdd(scalar, table.Modulus, table.BRedParams), table.Modulus, table.BRedParams), table.Modulus, table.MRedParams)
	}
}

// MulScalarAndAdd multiplies each coefficient of p1 by a scalar and adds the result on p2.
func (r *Ring) MulScalarAndAdd(p1 *Poly, scalar uint64, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		table := table
		modulus := table.Modulus
		bredParams := table.BRedParams
		mredParams := table.MRedParams
		MulScalarMontgomeryAndAddVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], MForm(BRedAdd(scalar, modulus, bredParams), modulus, bredParams), modulus, mredParams)
	}
}

// MulRNSScalarMontgomery multiplies p with a scalar value expressed in the CRT decomposition at a given level.
// It assumes the scalar decomposition to be in Montgomery form.
func (r *Ring) MulRNSScalarMontgomery(p *Poly, scalar RNSScalar, pOut *Poly) {
	for i, table := range r.Tables[:r.level+1] {
		Qi := table.Modulus
		scalar := scalar[i]
		p1tmp, p2tmp := p.Coeffs[i], pOut.Coeffs[i]
		mredParams := table.MRedParams
		MulScalarMontgomeryVec(p1tmp, p2tmp, scalar, Qi, mredParams)
	}
}

// MulScalarAndSub multiplies each coefficient of p1 by a scalar and subtracts the result on p2.
func (r *Ring) MulScalarAndSub(p1 *Poly, scalar uint64, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulScalarMontgomeryAndAddVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], MForm(table.Modulus-BRedAdd(scalar, table.Modulus, table.BRedParams), table.Modulus, table.BRedParams), table.Modulus, table.MRedParams)
	}
}

// MulScalarBigint multiplies each coefficient of p1 by a big.Int scalar
// and writes the result on p2.
func (r *Ring) MulScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	N := r.N()
	scalarQi := new(big.Int)
	for i, table := range r.Tables[:r.level+1] {
		scalarQi.Mod(scalar, NewUint(table.Modulus))
		MulScalarMontgomeryVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], MForm(BRedAdd(scalarQi.Uint64(), table.Modulus, table.BRedParams), table.Modulus, table.BRedParams), table.Modulus, table.MRedParams)
	}
}

// EvalPolyScalar evaluate the polynomial pol at pk and writes the result in p3
func (r *Ring) EvalPolyScalar(pol []*Poly, scalar uint64, pOut *Poly) {
	pOut.Copy(pol[len(pol)-1])
	for i := len(pol) - 1; i > 0; i-- {
		r.MulScalar(pOut, scalar, pOut)
		r.Add(pOut, pol[i-1], pOut)
	}
}

// Shift circularly shifts the coefficients of the polynomial p1 by k positions to the left and writes the result on p2.
func (r *Ring) Shift(p1 *Poly, k int, p2 *Poly) {
	for i := range p1.Coeffs {
		utils.RotateUint64SliceAllocFree(p1.Coeffs[i], k, p2.Coeffs[i])
	}
}

// MForm switches p1 to the Montgomery domain and writes the result on p2.
func (r *Ring) MForm(p1, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MFormVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], table.Modulus, table.BRedParams)
	}
}

// MFormConstant switches p1 to the Montgomery domain and writes the result on p2.
// Result is in the range [0, 2q-1]
func (r *Ring) MFormConstant(p1, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MFormConstantVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], table.Modulus, table.BRedParams)
	}
}

// InvMForm switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r *Ring) InvMForm(p1, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		InvMFormVec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulByPow2 multiplies p1 by 2^pow2 and writes the result on p2.
func (r *Ring) MulByPow2(p1 *Poly, pow2 int, p2 *Poly) {
	N := r.N()
	r.MForm(p1, p2)
	for i, table := range r.Tables[:r.level+1] {
		MulByPow2Vec(p1.Coeffs[i][:N], p2.Coeffs[i][:N], pow2, table.Modulus, table.MRedParams)
	}
}

// MultByMonomial multiplies p1 by x^monomialDeg and writes the result on p2.
func (r *Ring) MultByMonomial(p1 *Poly, monomialDeg int, p2 *Poly) {

	N := r.N()

	shift := monomialDeg % (N << 1)

	if shift == 0 {

		for i := range r.Tables {
			p1tmp, p2tmp := p1.Coeffs[i][:N], p2.Coeffs[i][:N]
			for j := 0; j < N; j++ {
				p2tmp[j] = p1tmp[j]
			}
		}

	} else {

		tmpx := r.NewPoly()

		if shift < N {

			for i := range r.Tables {
				p1tmp, tmpxT := p1.Coeffs[i][:N], tmpx.Coeffs[i]
				for j := 0; j < N; j++ {
					tmpxT[j] = p1tmp[j]
				}
			}

		} else {

			for i, table := range r.Tables {
				qi := table.Modulus
				p1tmp, tmpxT := p1.Coeffs[i][:N], tmpx.Coeffs[i]
				for j := 0; j < N; j++ {
					tmpxT[j] = qi - p1tmp[j]
				}
			}
		}

		shift %= N

		for i, table := range r.Tables {
			qi := table.Modulus
			p2tmp, tmpxT := p2.Coeffs[i][:N], tmpx.Coeffs[i]
			for j := 0; j < shift; j++ {
				p2tmp[j] = qi - tmpxT[N-shift+j]
			}
		}

		for i := range r.Tables {
			p2tmp, tmpxT := p2.Coeffs[i][:N], tmpx.Coeffs[i]
			for j := shift; j < N; j++ {
				p2tmp[j] = tmpxT[j-shift]

			}
		}
	}
}

// MulByVectorMontgomery multiplies p1 by a vector of uint64 coefficients and writes the result on p2.
func (r *Ring) MulByVectorMontgomery(p1 *Poly, vector []uint64, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryVec(p1.Coeffs[i][:N], vector, p2.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MulByVectorMontgomeryAndAddNoMod multiplies p1 by a vector of uint64 coefficients and adds the result on p2 without modular reduction.
func (r *Ring) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint64, p2 *Poly) {
	N := r.N()
	for i, table := range r.Tables[:r.level+1] {
		MulCoeffsMontgomeryAndAddNoModVec(p1.Coeffs[i][:N], vector, p2.Coeffs[i][:N], table.Modulus, table.MRedParams)
	}
}

// MapSmallDimensionToLargerDimensionNTT maps Y = X^{N/n} -> X directly in the NTT domain
func MapSmallDimensionToLargerDimensionNTT(polSmall, polLarge *Poly) {
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

// BitReverse applies a bit reverse permutation on the coefficients of p1 and writes the result on p2.
// In can safely be used for in-place permutation.
func (r *Ring) BitReverse(p1, p2 *Poly) {

	N := r.N()

	bitLenOfN := uint64(bits.Len64(uint64(N)) - 1)

	if p1 != p2 {
		for i := range r.Tables {
			p1tmp, p2tmp := p1.Coeffs[i][:N], p2.Coeffs[i][:N]
			for j := 0; j < N; j++ {
				p2tmp[utils.BitReverse64(uint64(j), bitLenOfN)] = p1tmp[j]
			}
		}
	} else { // In place in case p1 = p2
		for x := range r.Tables {
			p2tmp := p2.Coeffs[x]
			for i := 0; i < N; i++ {
				j := utils.BitReverse64(uint64(i), bitLenOfN)
				if i < int(j) {
					p2tmp[i], p2tmp[j] = p2tmp[j], p2tmp[i]
				}
			}
		}
	}
}

// Log2OfInnerSum returns the bit-size of the sum of all the coefficients (in absolute value) of a Poly.
func (r *Ring) Log2OfInnerSum(poly *Poly) (logSum int) {
	sumRNS := make([]uint64, r.level+1)
	var sum uint64
	for i, table := range r.Tables[:r.level+1] {

		qi := table.Modulus
		qiHalf := qi >> 1
		coeffs := poly.Coeffs[i]
		sum = 0

		for j := 0; j < r.N(); j++ {

			v := coeffs[j]

			if v >= qiHalf {
				sum = CRed(sum+qi-v, qi)
			} else {
				sum = CRed(sum+v, qi)
			}
		}

		sumRNS[i] = sum
	}

	var smallNorm = true
	for i := 1; i < r.level+1; i++ {
		smallNorm = smallNorm && (sumRNS[0] == sumRNS[i])
	}

	if !smallNorm {
		var crtReconstruction *big.Int

		sumBigInt := NewUint(0)
		QiB := new(big.Int)
		tmp := new(big.Int)
		modulusBigint := r.ModulusAtLevel[r.level]

		for i, table := range r.Tables[:r.level+1] {
			QiB.SetUint64(table.Modulus)
			crtReconstruction = new(big.Int).Quo(modulusBigint, QiB)
			tmp.ModInverse(crtReconstruction, QiB)
			tmp.Mod(tmp, QiB)
			crtReconstruction.Mul(crtReconstruction, tmp)
			sumBigInt.Add(sumBigInt, tmp.Mul(NewUint(sumRNS[i]), crtReconstruction))
		}

		sumBigInt.Mod(sumBigInt, modulusBigint)

		logSum = sumBigInt.BitLen()
	} else {
		logSum = bits.Len64(sumRNS[0])
	}

	return
}
