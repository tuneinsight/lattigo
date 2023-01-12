package ring

// Add evaluates p3 = p1 + p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) Add(p1, p2, p3 []uint64) {
	addvec(p1, p2, p3, s.Modulus)
}

// AddLazy evaluates p3 = p1 + p2.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddLazy(p1, p2, p3 []uint64) {
	addlazyvec(p1, p2, p3)
}

// Sub evaluates p3 = p1 - p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) Sub(p1, p2, p3 []uint64) {
	subvec(p1, p2, p3, s.Modulus)
}

// SubLazy evaluates p3 = p1 - p2.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) SubLazy(p1, p2, p3 []uint64) {
	sublazyvec(p1, p2, p3, s.Modulus)
}

// Neg evaluates p2 = -p1 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) Neg(p1, p2 []uint64) {
	negvec(p1, p2, s.Modulus)
}

// Reduce evaluates p2 = p1 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) Reduce(p1, p2 []uint64) {
	reducevec(p1, p2, s.Modulus, s.BRedConstant)
}

// ReduceLazy evaluates p2 = p1 (mod modulus) with p2 in range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) ReduceLazy(p1, p2 []uint64) {
	reducelazyvec(p1, p2, s.Modulus, s.BRedConstant)
}

// MulCoeffsLazy evaluates p3 = p1*p2.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsLazy(p1, p2, p3 []uint64) {
	mulcoeffslazyvec(p1, p2, p3)
}

// MulCoeffsLazyThenAddLazy evaluates p3 = p3 + p1*p2.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsLazyThenAddLazy(p1, p2, p3 []uint64) {
	mulcoeffslazythenaddlazyvec(p1, p2, p3)
}

// MulCoeffsBarrett evaluates p3 = p1*p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsBarrett(p1, p2, p3 []uint64) {
	mulcoeffsbarrettvec(p1, p2, p3, s.Modulus, s.BRedConstant)
}

// MulCoeffsBarrettLazy evaluates p3 = p1*p2 (mod modulus) with p3 in [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsBarrettLazy(p1, p2, p3 []uint64) {
	mulcoeffsbarrettlazyvec(p1, p2, p3, s.Modulus, s.BRedConstant)
}

// MulCoeffsBarrettThenAdd evaluates p3 = p3 + (p1*p2) (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsBarrettThenAdd(p1, p2, p3 []uint64) {
	mulcoeffsthenaddvec(p1, p2, p3, s.Modulus, s.BRedConstant)
}

// MulCoeffsBarrettThenAddLazy evaluates p3 = p3 + p1*p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsBarrettThenAddLazy(p1, p2, p3 []uint64) {
	mulcoeffsbarrettthenaddlazyvec(p1, p2, p3, s.Modulus, s.BRedConstant)
}

// MulCoeffsMontgomery evaluates p3 = p1*p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomery(p1, p2, p3 []uint64) {
	mulcoeffsmontgomeryvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// MulCoeffsMontgomeryLazy evaluates p3 = p1*p2 (mod modulus) with p3 in range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryLazy(p1, p2, p3 []uint64) {
	mulcoeffsmontgomerylazyvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// MulCoeffsMontgomeryThenAdd evaluates p3 = p3 + (p1*p2) (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryThenAdd(p1, p2, p3 []uint64) {
	mulcoeffsmontgomerythenaddvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// MulCoeffsMontgomeryThenAddLazy evaluates p3 = p3 + (p1*p2 (mod modulus)).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryThenAddLazy(p1, p2, p3 []uint64) {
	mulcoeffsmontgomerythenaddlazyvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// MulCoeffsMontgomeryLazyThenAddLazy evaluates p3 = p3 + p1*p2 (mod modulus) with p3 in range [0, 3modulus-2].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryLazyThenAddLazy(p1, p2, p3 []uint64) {
	mulcoeffsmontgomerylazythenaddlazyvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// MulCoeffsMontgomeryThenSub evaluates p3 = p3 - p1*p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryThenSub(p1, p2, p3 []uint64) {
	mulcoeffsmontgomerythensubvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// MulCoeffsMontgomeryThenSubLazy evaluates p3 = p3 - p1*p2 (mod modulus) with p3 in range [0, 2*modulus-2].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryThenSubLazy(p1, p2, p3 []uint64) {
	mulcoeffsmontgomerythensublazyvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// MulCoeffsMontgomeryLazyThenSubLazy evaluates p3 = p3 - p1*p2 (mod modulus) with p3 in range [0, 3*modulus-2].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryLazyThenSubLazy(p1, p2, p3 []uint64) {
	mulcoeffsmontgomerylazythensublazyvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// MulCoeffsMontgomeryLazyThenNeg evaluates p3 = - p1*p2 (mod modulus) with p3 in range [0, 2*modulus-2].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryLazyThenNeg(p1, p2, p3 []uint64) {
	mulcoeffsmontgomerylazythenNegvec(p1, p2, p3, s.Modulus, s.MRedConstant)
}

// AddLazyThenMulScalarMontgomery evaluates p3 = (p1+p2)*scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddLazyThenMulScalarMontgomery(p1, p2 []uint64, scalarMont uint64, p3 []uint64) {
	addlazythenmulscalarmontgomeryvec(p1, p2, scalarMont, p3, s.Modulus, s.MRedConstant)
}

// AddScalarLazyThenMulScalarMontgomery evaluates p3 = (scalarMont0+p2)*scalarMont1 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddScalarLazyThenMulScalarMontgomery(p1 []uint64, scalar0, scalarMont1 uint64, p2 []uint64) {
	addscalarlazythenmulscalarmontgomeryvec(p1, scalar0, scalarMont1, p2, s.Modulus, s.MRedConstant)
}

// AddScalar evaluates p2 = p1 + scalar (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddScalar(p1 []uint64, scalar uint64, p2 []uint64) {
	addscalarvec(p1, scalar, p2, s.Modulus)
}

// AddScalarLazy evaluates p2 = p1 + scalar.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddScalarLazy(p1 []uint64, scalar uint64, p2 []uint64) {
	addscalarlazyvec(p1, scalar, p2)
}

// AddScalarLazyThenNegTwoModulusLazy evaluates p2 = 2*modulus - p1 + scalar.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddScalarLazyThenNegTwoModulusLazy(p1 []uint64, scalar uint64, p2 []uint64) {
	addscalarlazythenNegTwoModuluslazyvec(p1, scalar, p2, s.Modulus)
}

// SubScalar evaluates p2 = p1 - scalar (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) SubScalar(p1 []uint64, scalar uint64, p2 []uint64) {
	subscalarvec(p1, scalar, p2, s.Modulus)
}

// MulScalarMontgomery evaluates p2 = p1*scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulScalarMontgomery(p1 []uint64, scalarMont uint64, p2 []uint64) {
	mulscalarmontgomeryvec(p1, scalarMont, p2, s.Modulus, s.MRedConstant)
}

// MulScalarMontgomeryLazy evaluates p2 = p1*scalarMont (mod modulus) with p2 in range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulScalarMontgomeryLazy(p1 []uint64, scalarMont uint64, p2 []uint64) {
	mulscalarmontgomerylazyvec(p1, scalarMont, p2, s.Modulus, s.MRedConstant)
}

// MulScalarMontgomeryThenAdd evaluates p2 = p2 + p1*scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulScalarMontgomeryThenAdd(p1 []uint64, scalarMont uint64, p2 []uint64) {
	mulscalarmontgomerythenaddvec(p1, scalarMont, p2, s.Modulus, s.MRedConstant)
}

// MulScalarMontgomeryThenAddScalar evaluates p2 = scalar + p1*scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulScalarMontgomeryThenAddScalar(p1 []uint64, scalar0, scalarMont1 uint64, p2 []uint64) {
	mulscalarmontgomerythenaddscalarvec(p1, scalar0, scalarMont1, p2, s.Modulus, s.MRedConstant)
}

// SubThenMulScalarMontgomeryTwoModulus evaluates p3 = (p1 + twomodulus - p2) * scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) SubThenMulScalarMontgomeryTwoModulus(p1, p2 []uint64, scalarMont uint64, p3 []uint64) {
	subthenmulscalarmontgomeryTwoModulusvec(p1, p2, scalarMont, p3, s.Modulus, s.MRedConstant)
}

// NTT evaluates p2 = NTT(p1).
func (s *SubRing) NTT(p1, p2 []uint64) {
	s.ntt.Forward(p1, p2)
}

// NTTLazy evaluates p2 = NTT(p1) with p2 in [0, 2*modulus-1].
func (s *SubRing) NTTLazy(p1, p2 []uint64) {
	s.ntt.ForwardLazy(p1, p2)
}

// INTT evaluates p2 = INTT(p1).
func (s *SubRing) INTT(p1, p2 []uint64) {
	s.ntt.Backward(p1, p2)
}

// INTTLazy evaluates p2 = INTT(p1) with p2 in [0, 2*modulus-1].
func (s *SubRing) INTTLazy(p1, p2 []uint64) {
	s.ntt.BackwardLazy(p1, p2)
}

// MForm evaluates p2 = p1 * 2^64 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MForm(p1, p2 []uint64) {
	mformvec(p1, p2, s.Modulus, s.BRedConstant)
}

// MFormLazy evaluates p2 = p1 * 2^64 (mod modulus) with p2 in the range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MFormLazy(p1, p2 []uint64) {
	mformlazyvec(p1, p2, s.Modulus, s.BRedConstant)
}

// IMForm evaluates p2 = p1 * (2^64)^-1 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) IMForm(p1, p2 []uint64) {
	imformvec(p1, p2, s.Modulus, s.MRedConstant)
}
