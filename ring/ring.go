// Package ring implements RNS-accelerated modular arithmetic operations for polynomials, including:
// RNS basis extension; RNS rescaling; number theoretic transform (NTT); uniform, Gaussian and ternary sampling.
package ring

import (
	//"bytes"
	//"encoding/gob"
	"encoding/json"
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// GaloisGen is an integer of order N/2 modulo M that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

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
	NumberTheoreticTransformer

	Tables []*Table

	// Product of the Moduli for each level
	ModulusAtLevel []*big.Int

	// Rescaling parameters (RNS division)
	RescaleParams [][]uint64

	level int
}

// ConjugateInvariantRing returns the conjugate invariant ring of the receiver ring.
// If `r.Type()==ConjugateInvariant`, then the method returns the receiver.
// if `r.Type()==Standard`, then the method returns a ring with ring degree N/2.
// The returned Ring is a shallow copy of the receiver.
func (r *Ring) ConjugateInvariantRing() (*Ring, error) {

	if r.Type() == ConjugateInvariant {
		return r, nil
	}

	cr := *r
	cr.Tables = make([]*Table, len(r.Tables))
	for i := range cr.Tables {
		cr.Tables[i] = NewTable(r.N()>>1, r.Tables[i].Modulus)
		cr.Tables[i].Factors = r.Tables[i].Factors
		if err := r.Tables[i].GenNTTParams(uint64(cr.N()) << 2); err != nil {
			return nil, err
		}
	}

	cr.NumberTheoreticTransformer = NumberTheoreticTransformerConjugateInvariant{}

	return &cr, nil
}

// StandardRing returns the standard ring of the receiver ring.
// If `r.Type()==Standard`, then the method returns the receiver.
// if `r.Type()==ConjugateInvariant`, then the method returns a ring with ring degree 2N.
// The returned Ring is a shallow copy of the receiver.
func (r *Ring) StandardRing() (*Ring, error) {
	if r.Type() == Standard {
		return r, nil
	}

	sr := *r
	sr.Tables = make([]*Table, len(r.Tables))
	for i := range sr.Tables {
		sr.Tables[i] = NewTable(r.N()<<1, r.Tables[i].Modulus)
		sr.Tables[i].Factors = r.Tables[i].Factors
		if err := r.Tables[i].GenNTTParams(uint64(sr.N()) << 1); err != nil {
			return nil, err
		}
	}

	sr.NumberTheoreticTransformer = NumberTheoreticTransformerStandard{}

	return &sr, nil
}

// N returns the ring degree.
func (r *Ring) N() int {
	return r.Tables[0].N
}

// NthRoot returns the multiplicative order of the primitive root.
func (r *Ring) NthRoot() uint64 {
	return r.Tables[0].NthRoot
}

// NbModuli returns the number of primes in the RNS basis of the ring.
func (r *Ring) NbModuli() int {
	return len(r.Tables)
}

// Level returns the level of the current ring.
func (r *Ring) Level() int {
	return r.level
}

// AtLevel returns a shallowcopy of the target ring that operates at the target level.
func (r *Ring) AtLevel(level int) *Ring {

	if level < 0 {
		panic("level cannot be negative")
	}

	if level > r.MaxLevel() {
		panic("level cannot be larger than max level")
	}

	return &Ring{
		NumberTheoreticTransformer: r.NumberTheoreticTransformer,
		Tables:                     r.Tables,
		ModulusAtLevel:             r.ModulusAtLevel,
		RescaleParams:              r.RescaleParams,
		level:                      level,
	}
}

// MaxLevel returns the maximum level allowed by the ring (#NbModuli -1).
func (r *Ring) MaxLevel() int {
	return r.NbModuli() - 1
}

// Moduli returns the list of primes in the modulus chain.
func (r *Ring) Moduli() (moduli []uint64) {
	moduli = make([]uint64, len(r.Tables))
	for i := range r.Tables {
		moduli[i] = r.Tables[i].Modulus
	}

	return
}

// Modulus returns the modulus of the target ring at the currently
// set level in *big.Int.
func (r *Ring) Modulus() *big.Int {
	return r.ModulusAtLevel[r.level]
}

// MRedParams returns the concatenation of the Montgomery parameters
// of the target ring.
func (r *Ring) MRedParams() (mredparams []uint64) {
	mredparams = make([]uint64, len(r.Tables))
	for i := range r.Tables {
		mredparams[i] = r.Tables[i].MRedParams
	}

	return
}

// MRedParams returns the concatenation of the Barrett parameters
// of the target ring.
func (r *Ring) BRedParams() (bredparams [][]uint64) {
	bredparams = make([][]uint64, len(r.Tables))
	for i := range r.Tables {
		bredparams[i] = r.Tables[i].BRedParams
	}

	return
}

// NewRing creates a new RNS Ring with degree N and coefficient moduli Moduli with Standard NTT. N must be a power of two larger than 8. Moduli should be
// a non-empty []uint64 with distinct prime elements. All moduli must also be equal to 1 modulo 2*N.
// An error is returned with a nil *Ring in the case of non NTT-enabling parameters.
func NewRing(N int, Moduli []uint64) (r *Ring, err error) {
	return NewRingWithCustomNTT(N, Moduli, NumberTheoreticTransformerStandard{}, 2*N)
}

// NewRingConjugateInvariant creates a new RNS Ring with degree N and coefficient moduli Moduli with Conjugate Invariant NTT. N must be a power of two larger than 8. Moduli should be
// a non-empty []uint64 with distinct prime elements. All moduli must also be equal to 1 modulo 4*N.
// An error is returned with a nil *Ring in the case of non NTT-enabling parameters.
func NewRingConjugateInvariant(N int, Moduli []uint64) (r *Ring, err error) {
	return NewRingWithCustomNTT(N, Moduli, NumberTheoreticTransformerConjugateInvariant{}, 4*N)
}

// NewRingFromType creates a new RNS Ring with degree N and coefficient moduli Moduli for which the type of NTT is determined by the ringType argument.
// If ringType==Standard, the ring is instantiated with standard NTT with the Nth root of unity 2*N. If ringType==ConjugateInvariant, the ring
// is instantiated with a ConjugateInvariant NTT with Nth root of unity 4*N. N must be a power of two larger than 8.
// Moduli should be a non-empty []uint64 with distinct prime elements. All moduli must also be equal to 1 modulo the root of unity.
// An error is returned with a nil *Ring in the case of non NTT-enabling parameters.
func NewRingFromType(N int, Moduli []uint64, ringType Type) (r *Ring, err error) {
	switch ringType {
	case Standard:
		return NewRingWithCustomNTT(N, Moduli, NumberTheoreticTransformerStandard{}, 2*N)
	case ConjugateInvariant:
		return NewRingWithCustomNTT(N, Moduli, NumberTheoreticTransformerConjugateInvariant{}, 4*N)
	default:
		return nil, fmt.Errorf("invalid ring type")
	}
}

// NewRingWithCustomNTT creates a new RNS Ring with degree N and coefficient moduli Moduli with user-defined NTT transform and primitive Nth root of unity.
// Moduli should be a non-empty []uint64 with distinct prime elements. All moduli must also be equal to 1 modulo the root of unity.
// N must be a power of two larger than 8. An error is returned with a nil *Ring in the case of non NTT-enabling parameters.
func NewRingWithCustomNTT(N int, Moduli []uint64, ntt NumberTheoreticTransformer, NthRoot int) (r *Ring, err error) {
	r = new(Ring)
	err = r.SetParameters(N, Moduli)
	if err != nil {
		return nil, err
	}

	r.NumberTheoreticTransformer = ntt

	err = r.GenNTTParams(uint64(NthRoot), nil, nil)
	if err != nil {
		return r, err
	}

	return r, nil
}

// Type returns the Type of the ring which might be either `Standard` or `ConjugateInvariant`.
func (r *Ring) Type() Type {
	switch r.NumberTheoreticTransformer.(type) {
	case NumberTheoreticTransformerStandard:
		return Standard
	case NumberTheoreticTransformerConjugateInvariant:
		return ConjugateInvariant
	default:
		panic("invalid NumberTheoreticTransformer type")
	}
}

// SetParameters initializes a *Ring by setting the required pre-computed values (except for the NTT-related values, which are set by the
// GenNTTParams function).
func (r *Ring) SetParameters(N int, Modulus []uint64) error {

	// Checks if N is a power of 2
	if (N < 16) || (N&(N-1)) != 0 && N != 0 {
		return errors.New("invalid ring degree (must be a power of 2 >= 8)")
	}

	if len(Modulus) == 0 {
		return errors.New("invalid modulus (must be a non-empty []uint64)")
	}

	if !utils.AllDistinct(Modulus) {
		return errors.New("invalid modulus (moduli are not distinct)")
	}

	// Computes bigQ for all levels
	r.ModulusAtLevel = make([]*big.Int, len(Modulus))
	r.ModulusAtLevel[0] = NewUint(Modulus[0])
	for i := 1; i < len(Modulus); i++ {
		r.ModulusAtLevel[i] = new(big.Int).Mul(r.ModulusAtLevel[i-1], NewUint(Modulus[i]))
	}

	r.Tables = make([]*Table, len(Modulus))

	for i := range r.Tables {
		r.Tables[i] = NewTable(N, Modulus[i])
	}

	r.RescaleParams = rewRescaleParams(r.Tables)

	r.level = len(Modulus) - 1

	return nil
}

func rewRescaleParams(tables []*Table) (rescaleParams [][]uint64) {

	rescaleParams = make([][]uint64, len(tables)-1)

	for j := len(tables) - 1; j > 0; j-- {

		qj := tables[j].Modulus

		rescaleParams[j-1] = make([]uint64, j)

		for i := 0; i < j; i++ {

			qi := tables[i].Modulus
			bredParams := tables[i].BRedParams

			rescaleParams[j-1][i] = MForm(qi-ModExp(qj, qi-2, qi), qi, bredParams)
		}
	}

	return
}

// GenNTTParams checks that N has been correctly initialized, and checks that each modulus is a prime congruent to 1 mod 2N (i.e. NTT-friendly).
// Then, it computes the variables required for the NTT. The purpose of ValidateParameters is to validate that the moduli allow the NTT, and to compute the
// NTT parameters.
func (r *Ring) GenNTTParams(NthRoot uint64, primitiveRoots []uint64, factors [][]uint64) (err error) {

	for i := range r.Tables {
		if primitiveRoots != nil && factors != nil {
			r.Tables[i].PrimitiveRoot = primitiveRoots[i]
			r.Tables[i].Factors = factors[i]
		}

		if err = r.Tables[i].GenNTTParams(NthRoot); err != nil {
			return
		}
	}

	return nil
}

// NewPoly creates a new polynomial with all coefficients set to 0.
func (r *Ring) NewPoly() *Poly {
	return NewPoly(r.N(), r.level)
}

// SetCoefficientsBigint sets the coefficients of p1 from an array of Int variables.
func (r *Ring) SetCoefficientsBigint(coeffs []*big.Int, p1 *Poly) {

	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, table := range r.Tables[:r.level+1] {

		QiBigint.SetUint64(table.Modulus)

		p1Coeffs := p1.Coeffs[i]

		for j, coeff := range coeffs {
			p1Coeffs[j] = coeffTmp.Mod(coeff, QiBigint).Uint64()
		}
	}
}

// PolyToString reconstructs p1 and returns the result in an array of string.
func (r *Ring) PolyToString(p1 *Poly) []string {

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
func (r *Ring) PolyToBigint(p1 *Poly, gap int, coeffsBigint []*big.Int) {

	crtReconstruction := make([]*big.Int, r.level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := r.ModulusAtLevel[r.level]

	for i, table := range r.Tables[:r.level+1] {
		QiB.SetUint64(table.Modulus)
		crtReconstruction[i] = new(big.Int).Quo(modulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	for i, j := 0, 0; j < r.N(); i, j = i+1, j+gap {

		tmp.SetUint64(0)
		coeffsBigint[i] = new(big.Int)

		for k := 0; k < r.level+1; k++ {
			coeffsBigint[i].Add(coeffsBigint[i], tmp.Mul(NewUint(p1.Coeffs[k][j]), crtReconstruction[k]))
		}

		coeffsBigint[i].Mod(coeffsBigint[i], modulusBigint)
	}
}

// PolyToBigintCentered reconstructs p1 and returns the result in an array of Int.
// Coefficients are centered around Q/2
// gap defines coefficients X^{i*gap} that will be reconstructed.
// For example, if gap = 1, then all coefficients are reconstructed, while
// if gap = 2 then only coefficients X^{2*i} are reconstructed.
func (r *Ring) PolyToBigintCentered(p1 *Poly, gap int, coeffsBigint []*big.Int) {

	crtReconstruction := make([]*big.Int, r.level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := r.ModulusAtLevel[r.level]

	for i, table := range r.Tables[:r.level+1] {
		QiB.SetUint64(table.Modulus)
		crtReconstruction[i] = new(big.Int).Quo(modulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	modulusBigintHalf := new(big.Int)
	modulusBigintHalf.Rsh(modulusBigint, 1)

	var sign int
	for i, j := 0, 0; j < r.N(); i, j = i+1, j+gap {

		tmp.SetUint64(0)
		coeffsBigint[i].SetUint64(0)

		for k := 0; k < r.level+1; k++ {
			coeffsBigint[i].Add(coeffsBigint[i], tmp.Mul(NewUint(p1.Coeffs[k][j]), crtReconstruction[k]))
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
func (r *Ring) Equal(p1, p2 *Poly) bool {

	for i := 0; i < r.level+1; i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	r.Reduce(p1, p1)
	r.Reduce(p2, p2)

	for i := 0; i < r.level+1; i++ {
		for j := 0; j < r.N(); j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}

// MarshalBinarySize returns the size in bytes of the target Ring.
func (r *Ring) MarshalBinarySize() (dataLen int) {
	dataLen++ // Type
	dataLen++ // #Tables
	dataLen++ // level
	for i := range r.Tables {
		dataLen += r.Tables[i].MarshalBinarySize()
	}

	return
}

// MarshalBinary encodes the target ring on a slice of bytes.
func (r *Ring) MarshalBinary() (data []byte, err error) {
	data = make([]byte, r.MarshalBinarySize())
	_, err = r.Encode(data)
	return
}

// UnmarshalBinary decodes a slice of bytes on the target ring.
func (r *Ring) UnmarshalBinary(data []byte) (err error) {
	var ptr int
	if ptr, err = r.Decode(data); err != nil {
		return
	}

	if ptr != len(data) {
		return fmt.Errorf("remaining unparsed data")
	}

	return
}

// Encode encodes the target Ring on a slice of bytes and returns
// the number of bytes written.
func (r *Ring) Encode(data []byte) (ptr int, err error) {
	data[ptr] = uint8(r.Type())
	ptr++
	data[ptr] = uint8(len(r.Tables))
	ptr++
	data[ptr] = uint8(r.level)
	ptr++

	var inc int
	for i := range r.Tables {
		if inc, err = r.Tables[i].Encode(data[ptr:]); err != nil {
			return
		}

		ptr += inc
	}

	return
}

// Decode decodes the input slice of bytes on the target Ring and
// returns the number of bytes read.
func (r *Ring) Decode(data []byte) (ptr int, err error) {
	ringType := Type(data[ptr])
	ptr++

	r.Tables = make([]*Table, data[ptr])
	ptr++

	r.level = int(data[ptr])
	ptr++

	var inc int
	for i := range r.Tables {

		r.Tables[i] = new(Table)

		if inc, err = r.Tables[i].Decode(data[ptr:]); err != nil {
			return
		}
		ptr += inc

		if i > 0 {
			if r.Tables[i].N != r.Tables[i-1].N || r.Tables[i].NthRoot != r.Tables[i-1].NthRoot {
				return ptr, fmt.Errorf("invalid tables: all tables must have the same ring degree and NthRoot")
			}
		}
	}

	switch ringType {
	case Standard:
		r.NumberTheoreticTransformer = NumberTheoreticTransformerStandard{}

		if int(r.NthRoot()) < r.N()<<1 {
			return ptr, fmt.Errorf("invalid ring type: NthRoot must be at least 2N but is %dN", int(r.NthRoot())/r.N())
		}

	case ConjugateInvariant:
		r.NumberTheoreticTransformer = NumberTheoreticTransformerConjugateInvariant{}

		if int(r.NthRoot()) < r.N()<<2 {
			return ptr, fmt.Errorf("invalid ring type: NthRoot must be at least 4N but is %dN", int(r.NthRoot())/r.N())
		}

	default:
		return ptr, fmt.Errorf("invalid ring type")
	}

	r.ModulusAtLevel = make([]*big.Int, len(r.Tables))

	r.ModulusAtLevel[0] = new(big.Int).SetUint64(r.Tables[0].Modulus)

	for i := 1; i < len(r.Tables); i++ {
		r.ModulusAtLevel[i] = new(big.Int).Mul(r.ModulusAtLevel[i-1], new(big.Int).SetUint64(r.Tables[i].Modulus))
	}

	r.RescaleParams = rewRescaleParams(r.Tables)

	return
}
