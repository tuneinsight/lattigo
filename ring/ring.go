// Package ring implements RNS-accelerated modular arithmetic operations for polynomials, including:
// RNS basis extension; RNS rescaling; number theoretic transform (NTT); uniform, Gaussian and ternary sampling.
package ring

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

const (
	// GaloisGen is an integer of order N/2 modulo M that spans Z_M with the integer -1.
	// The j-th ring automorphism takes the root zeta to zeta^(5j).
	GaloisGen uint64 = 5

	// MinimumRingDegreeForLoopUnrolledOperations is the minimum ring degree required to
	// safely perform loop-unrolled operations
	MinimumRingDegreeForLoopUnrolledOperations = 8
)

// Type is the type of ring used by the cryptographic scheme
type Type int

// RingStandard and RingConjugateInvariant are two types of Rings.
const (
	Standard           = Type(0) // Z[X]/(X^N + 1) (Default)
	ConjugateInvariant = Type(1) // Z[X+X^-1]/(X^2N + 1)
)

// String returns the string representation of the ring Type
func (rt Type) String() string {
	switch rt {
	case Standard:
		return "Standard"
	case ConjugateInvariant:
		return "ConjugateInvariant"
	default:
		return "Invalid"
	}
}

// UnmarshalJSON reads a JSON byte slice into the receiver Type
func (rt *Type) UnmarshalJSON(b []byte) error {
	var s string
	if err := json.Unmarshal(b, &s); err != nil {
		return err
	}
	switch s {
	default:
		return fmt.Errorf("invalid ring type: %s", s)
	case "Standard":
		*rt = Standard
	case "ConjugateInvariant":
		*rt = ConjugateInvariant
	}

	return nil
}

// MarshalJSON marshals the receiver Type into a JSON []byte
func (rt Type) MarshalJSON() ([]byte, error) {
	return json.Marshal(rt.String())
}

// Ring is a structure that keeps all the variables required to operate on a polynomial represented in this ring.
type Ring struct {
	SubRings []*SubRing

	// Product of the Moduli for each level
	ModulusAtLevel []*big.Int

	// Rescaling parameters (RNS division)
	RescaleConstants [][]uint64

	level int
}

// ConjugateInvariantRing returns the conjugate invariant ring of the receiver ring.
// If `r.Type()==ConjugateInvariant`, then the method returns the receiver.
// if `r.Type()==Standard`, then the method returns a ring with ring degree N/2.
// The returned Ring is a shallow copy of the receiver.
func (r Ring) ConjugateInvariantRing() (*Ring, error) {

	var err error

	if r.Type() == ConjugateInvariant {
		return &r, nil
	}

	cr := r

	cr.SubRings = make([]*SubRing, len(r.SubRings))

	factors := make([][]uint64, len(r.SubRings))

	for i, s := range r.SubRings {

		if cr.SubRings[i], err = NewSubRingWithCustomNTT(s.N>>1, s.Modulus, NewNumberTheoreticTransformerConjugateInvariant, int(s.NthRoot)); err != nil {
			return nil, err
		}

		factors[i] = s.Factors // Allocates factor for faster generation
	}

	return &cr, cr.generateNTTConstants(nil, factors)
}

// StandardRing returns the standard ring of the receiver ring.
// If `r.Type()==Standard`, then the method returns the receiver.
// if `r.Type()==ConjugateInvariant`, then the method returns a ring with ring degree 2N.
// The returned Ring is a shallow copy of the receiver.
func (r Ring) StandardRing() (*Ring, error) {

	var err error

	if r.Type() == Standard {
		return &r, nil
	}

	sr := r

	sr.SubRings = make([]*SubRing, len(r.SubRings))

	factors := make([][]uint64, len(r.SubRings))

	for i, s := range r.SubRings {

		if sr.SubRings[i], err = NewSubRingWithCustomNTT(s.N<<1, s.Modulus, NewNumberTheoreticTransformerStandard, int(s.NthRoot)); err != nil {
			return nil, err
		}

		factors[i] = s.Factors // Allocates factor for faster generation
	}

	return &sr, sr.generateNTTConstants(nil, factors)
}

// N returns the ring degree.
func (r Ring) N() int {
	return r.SubRings[0].N
}

// LogN returns log2(ring degree).
func (r Ring) LogN() int {
	return bits.Len64(uint64(r.N() - 1))
}

// LogModuli returns the size of the extended modulus P in bits
func (r Ring) LogModuli() (logmod float64) {
	for _, qi := range r.ModuliChain() {
		logmod += math.Log2(float64(qi))
	}
	return
}

// NthRoot returns the multiplicative order of the primitive root.
func (r Ring) NthRoot() uint64 {
	return r.SubRings[0].NthRoot
}

// ModuliChainLength returns the number of primes in the RNS basis of the ring.
func (r Ring) ModuliChainLength() int {
	return len(r.SubRings)
}

// Level returns the level of the current ring.
func (r Ring) Level() int {
	return r.level
}

// AtLevel returns an instance of the target ring that operates at the target level.
// This instance is thread safe and can be use concurrently with the base ring.
func (r Ring) AtLevel(level int) *Ring {

	// Sanity check
	if level < 0 {
		panic("level cannot be negative")
	}

	// Sanity check
	if level > r.MaxLevel() {
		panic("level cannot be larger than max level")
	}

	return &Ring{
		SubRings:         r.SubRings,
		ModulusAtLevel:   r.ModulusAtLevel,
		RescaleConstants: r.RescaleConstants,
		level:            level,
	}
}

// MaxLevel returns the maximum level allowed by the ring (#NbModuli -1).
func (r Ring) MaxLevel() int {
	return r.ModuliChainLength() - 1
}

// ModuliChain returns the list of primes in the modulus chain.
func (r Ring) ModuliChain() (moduli []uint64) {
	moduli = make([]uint64, len(r.SubRings))
	for i := range r.SubRings {
		moduli[i] = r.SubRings[i].Modulus
	}

	return
}

// Modulus returns the modulus of the target ring at the currently
// set level in *big.Int.
func (r Ring) Modulus() *big.Int {
	return r.ModulusAtLevel[r.level]
}

// MRedConstants returns the concatenation of the Montgomery constants
// of the target ring.
func (r Ring) MRedConstants() (MRC []uint64) {
	MRC = make([]uint64, len(r.SubRings))
	for i := range r.SubRings {
		MRC[i] = r.SubRings[i].MRedConstant
	}

	return
}

// BRedConstants returns the concatenation of the Barrett constants
// of the target ring.
func (r Ring) BRedConstants() (BRC [][2]uint64) {
	BRC = make([][2]uint64, len(r.SubRings))
	for i := range r.SubRings {
		BRC[i] = r.SubRings[i].BRedConstant
	}

	return
}

// NewRing creates a new RNS Ring with degree N and coefficient moduli Moduli with Standard NTT. N must be a power of two larger than 8. Moduli should be
// a non-empty []uint64 with distinct prime elements. All moduli must also be equal to 1 modulo 2*N.
// An error is returned with a nil *Ring in the case of non NTT-enabling parameters.
func NewRing(N int, Moduli []uint64) (r *Ring, err error) {
	return NewRingWithCustomNTT(N, Moduli, NewNumberTheoreticTransformerStandard, 2*N)
}

// NewRingConjugateInvariant creates a new RNS Ring with degree N and coefficient moduli Moduli with Conjugate Invariant NTT. N must be a power of two larger than 8. Moduli should be
// a non-empty []uint64 with distinct prime elements. All moduli must also be equal to 1 modulo 4*N.
// An error is returned with a nil *Ring in the case of non NTT-enabling parameters.
func NewRingConjugateInvariant(N int, Moduli []uint64) (r *Ring, err error) {
	return NewRingWithCustomNTT(N, Moduli, NewNumberTheoreticTransformerConjugateInvariant, 4*N)
}

// NewRingFromType creates a new RNS Ring with degree N and coefficient moduli Moduli for which the type of NTT is determined by the ringType argument.
// If ringType==Standard, the ring is instantiated with standard NTT with the Nth root of unity 2*N. If ringType==ConjugateInvariant, the ring
// is instantiated with a ConjugateInvariant NTT with Nth root of unity 4*N. N must be a power of two larger than 8.
// Moduli should be a non-empty []uint64 with distinct prime elements. All moduli must also be equal to 1 modulo the root of unity.
// An error is returned with a nil *Ring in the case of non NTT-enabling parameters.
func NewRingFromType(N int, Moduli []uint64, ringType Type) (r *Ring, err error) {
	switch ringType {
	case Standard:
		return NewRingWithCustomNTT(N, Moduli, NewNumberTheoreticTransformerStandard, 2*N)
	case ConjugateInvariant:
		return NewRingWithCustomNTT(N, Moduli, NewNumberTheoreticTransformerConjugateInvariant, 4*N)
	default:
		return nil, fmt.Errorf("invalid ring type")
	}
}

// NewRingWithCustomNTT creates a new RNS Ring with degree N and coefficient moduli Moduli with user-defined NTT transform and primitive Nth root of unity.
// ModuliChain should be a non-empty []uint64 with distinct prime elements.
// All moduli must also be equal to 1 modulo the root of unity.
// N must be a power of two larger than 8. An error is returned with a nil *Ring in the case of non NTT-enabling parameters.
func NewRingWithCustomNTT(N int, ModuliChain []uint64, ntt func(*SubRing, int) NumberTheoreticTransformer, NthRoot int) (r *Ring, err error) {
	r = new(Ring)

	// Checks if N is a power of 2
	if N < MinimumRingDegreeForLoopUnrolledOperations || (N&(N-1)) != 0 && N != 0 {
		return nil, fmt.Errorf("invalid ring degree: must be a power of 2 greater than %d", MinimumRingDegreeForLoopUnrolledOperations)
	}

	if len(ModuliChain) == 0 {
		return nil, fmt.Errorf("invalid ModuliChain (must be a non-empty []uint64)")
	}

	if !utils.AllDistinct(ModuliChain) {
		return nil, fmt.Errorf("invalid ModuliChain (moduli are not distinct)")
	}

	// Computes bigQ for all levels
	r.ModulusAtLevel = make([]*big.Int, len(ModuliChain))
	r.ModulusAtLevel[0] = bignum.NewInt(ModuliChain[0])
	for i := 1; i < len(ModuliChain); i++ {
		r.ModulusAtLevel[i] = new(big.Int).Mul(r.ModulusAtLevel[i-1], bignum.NewInt(ModuliChain[i]))
	}

	r.SubRings = make([]*SubRing, len(ModuliChain))

	for i := range r.SubRings {
		if r.SubRings[i], err = NewSubRingWithCustomNTT(N, ModuliChain[i], ntt, NthRoot); err != nil {
			return nil, err
		}
	}

	r.RescaleConstants = rewRescaleConstants(r.SubRings)

	r.level = len(ModuliChain) - 1

	return r, r.generateNTTConstants(nil, nil)
}

// Type returns the Type of the first subring which might be either `Standard` or `ConjugateInvariant`.
func (r *Ring) Type() Type {
	return r.SubRings[0].Type()
}

func rewRescaleConstants(subRings []*SubRing) (rescaleConstants [][]uint64) {

	rescaleConstants = make([][]uint64, len(subRings)-1)

	for j := len(subRings) - 1; j > 0; j-- {

		qj := subRings[j].Modulus

		rescaleConstants[j-1] = make([]uint64, j)

		for i := 0; i < j; i++ {
			qi := subRings[i].Modulus
			rescaleConstants[j-1][i] = MForm(qi-ModExp(qj, qi-2, qi), qi, subRings[i].BRedConstant)
		}
	}

	return
}

// generateNTTConstants checks that N has been correctly initialized, and checks that each modulus is a prime congruent to 1 mod 2N (i.e. NTT-friendly).
// Then, it computes the variables required for the NTT. The purpose of ValidateParameters is to validate that the moduli allow the NTT, and to compute the
// NTT parameters.
func (r *Ring) generateNTTConstants(primitiveRoots []uint64, factors [][]uint64) (err error) {

	for i := range r.SubRings {

		if primitiveRoots != nil && factors != nil {
			r.SubRings[i].PrimitiveRoot = primitiveRoots[i]
			r.SubRings[i].Factors = factors[i]
		}

		if err = r.SubRings[i].generateNTTConstants(); err != nil {
			return
		}
	}

	return nil
}

// NewPoly creates a new polynomial with all coefficients set to 0.
func (r Ring) NewPoly() Poly {
	return NewPoly(r.N(), r.level)
}

// NewMonomialXi returns a polynomial X^{i}.
func (r Ring) NewMonomialXi(i int) (p Poly) {

	p = r.NewPoly()

	N := r.N()

	i &= (N << 1) - 1

	if i >= N {
		i -= N << 1
	}

	for k, s := range r.SubRings[:r.level+1] {

		if i < 0 {
			p.Coeffs[k][N+i] = s.Modulus - 1
		} else {
			p.Coeffs[k][i] = 1
		}
	}

	return
}

// SetCoefficientsBigint sets the coefficients of p1 from an array of Int variables.
func (r Ring) SetCoefficientsBigint(coeffs []*big.Int, p1 Poly) {

	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, table := range r.SubRings[:r.level+1] {

		QiBigint.SetUint64(table.Modulus)

		p1Coeffs := p1.Coeffs[i]

		for j, coeff := range coeffs {
			p1Coeffs[j] = coeffTmp.Mod(coeff, QiBigint).Uint64()
		}
	}
}

// PolyToString reconstructs p1 and returns the result in an array of string.
func (r Ring) PolyToString(p1 Poly) []string {

	coeffsBigint := make([]*big.Int, r.N())
	r.PolyToBigint(p1, 1, coeffsBigint)
	coeffsString := make([]string, len(coeffsBigint))

	for i := range coeffsBigint {
		coeffsString[i] = coeffsBigint[i].String()
	}

	return coeffsString
}

// PolyToBigint reconstructs p1 and returns the result in an array of Int.
// gap defines coefficients X^{i*gap} that will be reconstructed.
// For example, if gap = 1, then all coefficients are reconstructed, while
// if gap = 2 then only coefficients X^{2*i} are reconstructed.
func (r Ring) PolyToBigint(p1 Poly, gap int, coeffsBigint []*big.Int) {

	crtReconstruction := make([]*big.Int, r.level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := r.ModulusAtLevel[r.level]

	for i, table := range r.SubRings[:r.level+1] {
		QiB.SetUint64(table.Modulus)
		crtReconstruction[i] = new(big.Int).Quo(modulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	N := r.N()

	for i, j := 0, 0; j < N; i, j = i+1, j+gap {

		tmp.SetUint64(0)
		coeffsBigint[i] = new(big.Int)

		for k := 0; k < r.level+1; k++ {
			coeffsBigint[i].Add(coeffsBigint[i], tmp.Mul(bignum.NewInt(p1.Coeffs[k][j]), crtReconstruction[k]))
		}

		coeffsBigint[i].Mod(coeffsBigint[i], modulusBigint)
	}
}

// PolyToBigintCentered reconstructs p1 and returns the result in an array of Int.
// Coefficients are centered around Q/2
// gap defines coefficients X^{i*gap} that will be reconstructed.
// For example, if gap = 1, then all coefficients are reconstructed, while
// if gap = 2 then only coefficients X^{2*i} are reconstructed.
func (r Ring) PolyToBigintCentered(p1 Poly, gap int, coeffsBigint []*big.Int) {

	crtReconstruction := make([]*big.Int, r.level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := r.ModulusAtLevel[r.level]

	for i, table := range r.SubRings[:r.level+1] {
		QiB.SetUint64(table.Modulus)
		crtReconstruction[i] = new(big.Int).Quo(modulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	modulusBigintHalf := new(big.Int)
	modulusBigintHalf.Rsh(modulusBigint, 1)

	N := r.N()

	var sign int
	for i, j := 0, 0; j < N; i, j = i+1, j+gap {

		tmp.SetUint64(0)
		coeffsBigint[i].SetUint64(0)

		for k := 0; k < r.level+1; k++ {
			coeffsBigint[i].Add(coeffsBigint[i], tmp.Mul(bignum.NewInt(p1.Coeffs[k][j]), crtReconstruction[k]))
		}

		coeffsBigint[i].Mod(coeffsBigint[i], modulusBigint)

		// Centers the coefficients
		sign = coeffsBigint[i].Cmp(modulusBigintHalf)

		if sign == 1 || sign == 0 {
			coeffsBigint[i].Sub(coeffsBigint[i], modulusBigint)
		}
	}
}

// Equal checks if p1 = p2 in the given Ring.
func (r Ring) Equal(p1, p2 Poly) bool {

	for i := 0; i < r.level+1; i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	r.Reduce(p1, p1)
	r.Reduce(p2, p2)

	return p1.Equal(&p2)
}

// ringParametersLiteral is a struct to store the minimum information
// to uniquely identify a Ring and be able to reconstruct it efficiently.
// This struct's purpose is to facilitate the marshalling of Rings.
type ringParametersLiteral []subRingParametersLiteral

// parametersLiteral returns the RingParametersLiteral of the Ring.
func (r Ring) parametersLiteral() ringParametersLiteral {
	p := make([]subRingParametersLiteral, len(r.SubRings))

	for i, s := range r.SubRings {
		p[i] = s.parametersLiteral()
	}

	return ringParametersLiteral(p)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (r Ring) MarshalBinary() (data []byte, err error) {
	return r.MarshalJSON()
}

// UnmarshalBinary decodes a slice of bytes generated by MarshalBinary or MarshalJSON on the object.
func (r *Ring) UnmarshalBinary(data []byte) (err error) {
	return r.UnmarshalJSON(data)
}

// MarshalJSON encodes the object into a binary form on a newly allocated slice of bytes with the json codec.
func (r Ring) MarshalJSON() (data []byte, err error) {
	return json.Marshal(r.parametersLiteral())
}

// UnmarshalJSON decodes a slice of bytes generated by MarshalJSON or MarshalBinary on the object.
func (r *Ring) UnmarshalJSON(data []byte) (err error) {

	p := ringParametersLiteral{}

	if err = json.Unmarshal(data, &p); err != nil {
		return
	}

	var rr *Ring
	if rr, err = newRingFromparametersLiteral(p); err != nil {
		return
	}

	*r = *rr

	return
}

// newRingFromparametersLiteral creates a new Ring from the provided RingParametersLiteral.
func newRingFromparametersLiteral(p ringParametersLiteral) (r *Ring, err error) {

	r = new(Ring)

	r.SubRings = make([]*SubRing, len(p))

	r.level = len(p) - 1

	for i := range r.SubRings {

		if r.SubRings[i], err = newSubRingFromParametersLiteral(p[i]); err != nil {
			return
		}

		if i > 0 {
			if r.SubRings[i].N != r.SubRings[i-1].N || r.SubRings[i].NthRoot != r.SubRings[i-1].NthRoot {
				return nil, fmt.Errorf("invalid SubRings: all SubRings must have the same ring degree and NthRoot")
			}
		}
	}

	r.ModulusAtLevel = make([]*big.Int, len(r.SubRings))

	r.ModulusAtLevel[0] = new(big.Int).SetUint64(r.SubRings[0].Modulus)

	for i := 1; i < len(r.SubRings); i++ {
		r.ModulusAtLevel[i] = new(big.Int).Mul(r.ModulusAtLevel[i-1], new(big.Int).SetUint64(r.SubRings[i].Modulus))
	}

	r.RescaleConstants = rewRescaleConstants(r.SubRings)

	return
}

// Log2OfStandardDeviation returns base 2 logarithm of the standard deviation of the coefficients
// of the polynomial.
func (r Ring) Log2OfStandardDeviation(poly Poly) (std float64) {

	N := r.N()

	prec := uint(128)

	coeffs := make([]*big.Int, N)

	for i := 0; i < N; i++ {
		coeffs[i] = new(big.Int)
	}

	r.PolyToBigintCentered(poly, 1, coeffs)

	mean := bignum.NewFloat(0, prec)
	tmp := bignum.NewFloat(0, prec)

	for i := 0; i < N; i++ {
		mean.Add(mean, tmp.SetInt(coeffs[i]))
	}

	mean.Quo(mean, bignum.NewFloat(float64(N), prec))

	stdFloat := bignum.NewFloat(0, prec)

	for i := 0; i < N; i++ {
		tmp.SetInt(coeffs[i])
		tmp.Sub(tmp, mean)
		tmp.Mul(tmp, tmp)
		stdFloat.Add(stdFloat, tmp)
	}

	stdFloat.Quo(stdFloat, bignum.NewFloat(float64(N-1), prec))

	stdFloat.Sqrt(stdFloat)

	stdF64, _ := stdFloat.Float64()

	return math.Log2(stdF64)
}
