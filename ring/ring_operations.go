package ring

import (
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/utils"
)

func (r *Ring) minLevelTernary(p1, p2, p3 *Poly) int {
	return utils.MinInt(utils.MinInt(len(r.Modulus)-1, p1.Level()), utils.MinInt(p2.Level(), p3.Level()))
}

func (r *Ring) minLevelBinary(p1, p2 *Poly) int {
	return utils.MinInt(utils.MinInt(len(r.Modulus)-1, p1.Level()), p2.Level())
}

// Add adds p1 to p2 coefficient-wise and writes the result on p3.
func (r *Ring) Add(p1, p2, p3 *Poly) {
	r.AddLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// AddLvl adds p1 to p2 coefficient-wise for the moduli from
// q_0 up to q_level and writes the result on p3.
func (r *Ring) AddLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		AddVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i])
	}
}

// AddNoMod adds p1 to p2 coefficient-wise without
// modular reduction and writes the result on p3.
func (r *Ring) AddNoMod(p1, p2, p3 *Poly) {
	r.AddNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// AddNoModLvl adds p1 to p2 coefficient-wise without modular reduction
// for the moduli from q_0 up to q_level and writes the result on p3.
func (r *Ring) AddNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		AddVecNoMod(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N])
	}
}

// Sub subtracts p2 to p1 coefficient-wise and writes the result on p3.
func (r *Ring) Sub(p1, p2, p3 *Poly) {
	r.SubLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// SubLvl subtracts p2 to p1 coefficient-wise and writes the result on p3.
func (r *Ring) SubLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		SubVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i])
	}
}

// SubNoMod subtracts p2 to p1 coefficient-wise without
// modular reduction and returns the result on p3.
func (r *Ring) SubNoMod(p1, p2, p3 *Poly) {
	r.SubNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// SubNoModLvl subtracts p2 to p1 coefficient-wise without modular reduction
// for the moduli from q_0 up to q_level and writes the result on p3.
func (r *Ring) SubNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		SubVecNomod(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i])
	}
}

// Neg sets all coefficients of p1 to their additive inverse and writes the result on p2.
func (r *Ring) Neg(p1, p2 *Poly) {
	r.NegLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// NegLvl sets the coefficients of p1 to their additive inverse for
// the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) NegLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		NegVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], r.Modulus[i])
	}
}

// Reduce applies a modular reduction on the coefficients of p1 and writes the result on p2.
func (r *Ring) Reduce(p1, p2 *Poly) {
	r.ReduceLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// ReduceLvl applies a modular reduction on the coefficients of p1
// for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) ReduceLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		ReduceVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], r.Modulus[i], r.BredParams[i])
	}
}

// ReduceConstant applies a modular reduction on the coefficients of p1 and writes the result on p2.
// Return values in [0, 2q-1]
func (r *Ring) ReduceConstant(p1, p2 *Poly) {
	r.ReduceConstantLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// ReduceConstantLvl applies a modular reduction on the coefficients of p1
// for the moduli from q_0 up to q_level and writes the result on p2.
// Return values in [0, 2q-1]
func (r *Ring) ReduceConstantLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		ReduceConstantVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], r.Modulus[i], r.BredParams[i])
	}
}

// Mod applies a modular reduction by m on the coefficients of p1 and writes the result on p2.
func (r *Ring) Mod(p1 *Poly, m uint64, p2 *Poly) {
	r.ModLvl(r.minLevelBinary(p1, p2), p1, m, p2)
}

// ModLvl applies a modular reduction by m on the coefficients of p1 and writes the result on p2.
func (r *Ring) ModLvl(level int, p1 *Poly, m uint64, p2 *Poly) {
	bredParams := BRedParams(m)
	for i := 0; i < level+1; i++ {
		ModVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], m, bredParams)
	}
}

// MulCoeffs multiplies p1 by p2 coefficient-wise, performs a
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffs(p1, p2, p3 *Poly) {
	r.MulCoeffsLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsLvl multiplies p1 by p2 coefficient-wise, performs a
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.BredParams[i])
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient-wise with
// a Barret modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsAndAdd(p1, p2, p3 *Poly) {
	r.MulCoeffsAndAddLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsAndAddLvl multiplies p1 by p2 coefficient-wise with
// a Barret modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsAndAddLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsAndAddVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.BredParams[i])
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient-wise with a Barrett
// modular reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {
	r.MulCoeffsAndAddNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsAndAddNoModLvl multiplies p1 by p2 coefficient-wise with a Barrett
// modular reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsAndAddNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsAndAddNoModVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.BredParams[i])
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and returns the result on p3.
func (r *Ring) MulCoeffsMontgomery(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction for the moduli from q_0 up to q_level and returns the result on p3.
func (r *Ring) MulCoeffsMontgomeryLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsMontgomeryConstant multiplies p1 by p2 coefficient-wise with a
// constant-time Montgomery modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsMontgomeryConstant(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryConstantLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryConstantLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction for the moduli from q_0 up to q_level and returns the result on p3.
func (r *Ring) MulCoeffsMontgomeryConstantLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryConstantVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsMontgomeryConstantAndNegLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction for the moduli from q_0 up to q_level and returns the negative result on p3.
func (r *Ring) MulCoeffsMontgomeryConstantAndNegLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryConstantAndNeg(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryAndAddLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryAndAddLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction for the moduli from q_0 up to q_level and adds the result to p3.
func (r *Ring) MulCoeffsMontgomeryAndAddLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryAndAddVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryAndAddNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryAndAddNoModLvl multiplies p1 by p2 coefficient-wise with a Montgomery modular
// reduction for the moduli from q_0 up to q_level and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndAddNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryAndAddNoModVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsMontgomeryConstantAndAddNoMod multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3 without modular reduction.
// Return values in [0, 3q-1]
func (r *Ring) MulCoeffsMontgomeryConstantAndAddNoMod(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryConstantAndAddNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryConstantAndAddNoModLvl multiplies p1 by p2 coefficient-wise with a constant-time Montgomery
// modular reduction for the moduli from q_0 up to q_level and adds the result to p3 without modular reduction.
// Return values in [0, 3q-1]
func (r *Ring) MulCoeffsMontgomeryConstantAndAddNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryConstantAndAddNoModVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r *Ring) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryAndSubLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryAndSubLvl multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r *Ring) MulCoeffsMontgomeryAndSubLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryAndSubVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and subtracts the result from p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryAndSubNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryAndSubNoModLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and subtracts the result from p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndSubNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryAndSubNoMod(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsMontgomeryConstantAndSubNoModLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and subtracts the result from p3 without modular reduction.
// Return values in [0, 3q-1]
func (r *Ring) MulCoeffsMontgomeryConstantAndSubNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryConstantAndSubNoMod(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulCoeffsConstant multiplies p1 by p2 coefficient-wise with a constant-time
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsConstant(p1, p2, p3 *Poly) {
	r.MulCoeffsConstantLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsConstantLvl multiplies p1 by p2 coefficient-wise with a constant-time
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsConstantLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsConstantVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], p3.Coeffs[i][:r.N], r.Modulus[i], r.BredParams[i])
	}
}

// AddScalar adds a scalar to each coefficient of p1 and writes the result on p2.
func (r *Ring) AddScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	r.AddScalarLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// AddScalarLvl adds a scalar to each coefficient of p1 and writes the result on p2.
func (r *Ring) AddScalarLvl(level int, p1 *Poly, scalar uint64, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		AddScalarVec(p1.Coeffs[i][:r.N], p1.Coeffs[i][:r.N], scalar, r.Modulus[i])
	}
}

// AddScalarBigint adds a big.Int scalar to each coefficient of p1 and writes the result on p2.
func (r *Ring) AddScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	r.AddScalarBigintLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// AddScalarBigintLvl adds a big.Int scalar to each coefficient of p1 and writes the result on p2.
func (r *Ring) AddScalarBigintLvl(level int, p1 *Poly, scalar *big.Int, p2 *Poly) {
	tmp := new(big.Int)
	for i := 0; i < level+1; i++ {
		AddScalarVec(p1.Coeffs[i][:r.N], p1.Coeffs[i][:r.N], tmp.Mod(scalar, NewUint(r.Modulus[i])).Uint64(), r.Modulus[i])
	}
}

// SubScalar subtracts a scalar from each coefficient of p1 and writes the result on p2.
func (r *Ring) SubScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	r.SubScalarLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// SubScalarLvl subtracts a scalar from each coefficient of p1 and writes the result on p2.
func (r *Ring) SubScalarLvl(level int, p1 *Poly, scalar uint64, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		SubScalarVec(p1.Coeffs[i][:r.N], p1.Coeffs[i][:r.N], scalar, r.Modulus[i])
	}
}

// SubScalarBigint subtracts a big.Int scalar from each coefficient of p1 and writes the result on p2.
func (r *Ring) SubScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	r.SubScalarBigintLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// SubScalarBigintLvl subtracts a big.Int scalar from each coefficient of p1 and writes the result on p2.
func (r *Ring) SubScalarBigintLvl(level int, p1 *Poly, scalar *big.Int, p2 *Poly) {
	tmp := new(big.Int)
	for i := 0; i < level+1; i++ {
		SubScalarVec(p1.Coeffs[i][:r.N], p1.Coeffs[i][:r.N], tmp.Mod(scalar, NewUint(r.Modulus[i])).Uint64(), r.Modulus[i])
	}
}

// MulScalar multiplies each coefficient of p1 by a scalar and writes the result on p2.
func (r *Ring) MulScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	r.MulScalarLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// MulScalarLvl multiplies each coefficient of p1 by a scalar for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) MulScalarLvl(level int, p1 *Poly, scalar uint64, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		MulScalarMontgomeryVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], MForm(BRedAdd(scalar, r.Modulus[i], r.BredParams[i]), r.Modulus[i], r.BredParams[i]), r.Modulus[i], r.MredParams[i])
	}
}

// MulScalarAndAdd multiplies each coefficient of p1 by a scalar and adds the result on p2.
func (r *Ring) MulScalarAndAdd(p1 *Poly, scalar uint64, p2 *Poly) {
	r.MulScalarAndAddLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// MulScalarAndAddLvl multiplies each coefficient of p1 by a scalar for the moduli from q_0 up to q_level and adds the result on p2.
func (r *Ring) MulScalarAndAddLvl(level int, p1 *Poly, scalar uint64, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		MulScalarMontgomeryAndAddVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], MForm(BRedAdd(scalar, r.Modulus[i], r.BredParams[i]), r.Modulus[i], r.BredParams[i]), r.Modulus[i], r.MredParams[i])
	}
}

// MulScalarBigint multiplies each coefficient of p1 by a big.Int scalar and writes the result on p2.
func (r *Ring) MulScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	r.MulScalarBigintLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// MulScalarBigintLvl multiplies each coefficient of p1 by a big.Int scalar
//for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) MulScalarBigintLvl(level int, p1 *Poly, scalar *big.Int, p2 *Poly) {
	scalarQi := new(big.Int)
	for i := 0; i < level+1; i++ {
		scalarQi.Mod(scalar, NewUint(r.Modulus[i]))
		MulScalarMontgomeryVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], MForm(BRedAdd(scalarQi.Uint64(), r.Modulus[i], r.BredParams[i]), r.Modulus[i], r.BredParams[i]), r.Modulus[i], r.MredParams[i])
	}
}

// Shift circulary shifts the coefficients of the polynomial p1 by n positions to the left and writes the result on p2.
func (r *Ring) Shift(p1 *Poly, n int, p2 *Poly) {
	mask := (1 << r.N) - 1
	for i := range r.Modulus {
		p2.Coeffs[i] = append(p1.Coeffs[i][(n&mask):], p1.Coeffs[i][:(n&mask)]...)
	}
}

// MForm switches p1 to the Montgomery domain and writes the result on p2.
func (r *Ring) MForm(p1, p2 *Poly) {
	r.MFormLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// MFormLvl switches p1 to the Montgomery domain for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) MFormLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		MFormVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], r.Modulus[i], r.BredParams[i])
	}
}

// MFormConstantLvl switches p1 to the Montgomery domain for the moduli from q_0 up to q_level and writes the result on p2.
// Result is in the range [0, 2q-1]
func (r *Ring) MFormConstantLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		MFormConstantVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], r.Modulus[i], r.BredParams[i])
	}
}

// InvMForm switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r *Ring) InvMForm(p1, p2 *Poly) {
	r.InvMFormLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// InvMFormLvl switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r *Ring) InvMFormLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		InvMFormVec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulByPow2New multiplies p1 by 2^pow2 and returns the result in a new polynomial p2.
func (r *Ring) MulByPow2New(p1 *Poly, pow2 int) (p2 *Poly) {
	p2 = r.NewPoly()
	r.MulByPow2(p1, pow2, p2)
	return
}

// MulByPow2 multiplies p1 by 2^pow2 and writes the result on p2.
func (r *Ring) MulByPow2(p1 *Poly, pow2 int, p2 *Poly) {
	r.MulByPow2Lvl(r.minLevelBinary(p1, p2), p1, pow2, p2)
}

// MulByPow2Lvl multiplies p1 by 2^pow2 for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) MulByPow2Lvl(level int, p1 *Poly, pow2 int, p2 *Poly) {
	r.MFormLvl(level, p1, p2)
	for i := 0; i < level+1; i++ {
		MulByPow2Vec(p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N], pow2, r.Modulus[i], r.MredParams[i])
	}
}

// MultByMonomialNew multiplies p1 by x^monomialDeg and writes the result on a new polynomial p2.
func (r *Ring) MultByMonomialNew(p1 *Poly, monomialDeg int) (p2 *Poly) {
	p2 = r.NewPoly()
	r.MultByMonomial(p1, monomialDeg, p2)
	return
}

// MultByMonomial multiplies p1 by x^monomialDeg and writes the result on p2.
func (r *Ring) MultByMonomial(p1 *Poly, monomialDeg int, p2 *Poly) {

	shift := monomialDeg % (r.N << 1)

	if shift == 0 {

		for i := range r.Modulus {
			p1tmp, p2tmp := p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N]
			for j := 0; j < r.N; j++ {
				p2tmp[j] = p1tmp[j]
			}
		}

	} else {

		tmpx := r.NewPoly()

		if shift < r.N {

			for i := range r.Modulus {
				p1tmp, tmpxT := p1.Coeffs[i][:r.N], tmpx.Coeffs[i]
				for j := 0; j < r.N; j++ {
					tmpxT[j] = p1tmp[j]
				}
			}

		} else {

			for i, qi := range r.Modulus {
				p1tmp, tmpxT := p1.Coeffs[i][:r.N], tmpx.Coeffs[i]
				for j := 0; j < r.N; j++ {
					tmpxT[j] = qi - p1tmp[j]
				}
			}
		}

		shift %= r.N

		for i, qi := range r.Modulus {
			p2tmp, tmpxT := p2.Coeffs[i][:r.N], tmpx.Coeffs[i]
			for j := 0; j < shift; j++ {
				p2tmp[j] = qi - tmpxT[r.N-shift+j]
			}
		}

		for i := range r.Modulus {
			p2tmp, tmpxT := p2.Coeffs[i][:r.N], tmpx.Coeffs[i]
			for j := shift; j < r.N; j++ {
				p2tmp[j] = tmpxT[j-shift]

			}
		}
	}
}

// MulByVectorMontgomery multiplies p1 by a vector of uint64 coefficients and writes the result on p2.
func (r *Ring) MulByVectorMontgomery(p1 *Poly, vector []uint64, p2 *Poly) {
	r.MulByVectorMontgomeryLvl(r.minLevelBinary(p1, p2), p1, vector, p2)
}

// MulByVectorMontgomeryLvl multiplies p1 by a vector of uint64 coefficients and writes the result on p2.
func (r *Ring) MulByVectorMontgomeryLvl(level int, p1 *Poly, vector []uint64, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryVec(p1.Coeffs[i][:r.N], vector, p2.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
	}
}

// MulByVectorMontgomeryAndAddNoMod multiplies p1 by a vector of uint64 coefficients and adds the result on p2 without modular reduction.
func (r *Ring) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint64, p2 *Poly) {
	r.MulByVectorMontgomeryAndAddNoModLvl(r.minLevelBinary(p1, p2), p1, vector, p2)
}

// MulByVectorMontgomeryAndAddNoModLvl multiplies p1 by a vector of uint64 coefficients and adds the result on p2 without modular reduction.
func (r *Ring) MulByVectorMontgomeryAndAddNoModLvl(level int, p1 *Poly, vector []uint64, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		MulCoeffsMontgomeryAndAddNoModVec(p1.Coeffs[i][:r.N], vector, p2.Coeffs[i][:r.N], r.Modulus[i], r.MredParams[i])
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
	bitLenOfN := uint64(bits.Len64(uint64(r.N)) - 1)

	if p1 != p2 {
		for i := range r.Modulus {
			p1tmp, p2tmp := p1.Coeffs[i][:r.N], p2.Coeffs[i][:r.N]
			for j := 0; j < r.N; j++ {
				p2tmp[utils.BitReverse64(uint64(j), bitLenOfN)] = p1tmp[j]
			}
		}
	} else { // In place in case p1 = p2
		for x := range r.Modulus {
			p2tmp := p2.Coeffs[x]
			for i := 0; i < r.N; i++ {
				j := utils.BitReverse64(uint64(i), bitLenOfN)
				if i < int(j) {
					p2tmp[i], p2tmp[j] = p2tmp[j], p2tmp[i]
				}
			}
		}
	}
}

// Rotate applies a Galois automorphism on p1 in NTT form,
// rotating the coefficients to the right by n positions, and writes the result on p2.
// It requires the data to be permuted in bit-reversal order before applying the NTT.
func (r *Ring) Rotate(p1 *Poly, n int, p2 *Poly) {

	var root, gal uint64

	n &= (1 << r.N) - 1

	for i, qi := range r.Modulus {

		mredParams := r.MredParams[i]

		root = MRed(r.PsiMont[i], r.PsiMont[i], qi, mredParams)

		root = ModexpMontgomery(root, n, qi, mredParams, r.BredParams[i])

		gal = MForm(1, qi, r.BredParams[i])

		p1tmp, p2tmp := p1.Coeffs[i][:r.N], p1.Coeffs[i][:r.N]

		for j := 1; j < r.N; j++ {

			gal = MRed(gal, root, qi, mredParams)

			p2tmp[j] = MRed(p1tmp[j], gal, qi, mredParams)

		}
	}
}
