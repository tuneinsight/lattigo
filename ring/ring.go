// Package ring implements RNS-accelerated modular arithmetic operations for polynomials, including:
// RNS basis extension; RNS rescaling; number theoretic transform (NTT); uniform, Gaussian and ternary sampling.
package ring

import (
	"bytes"
	"encoding/gob"
	"encoding/json"
	"errors"
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/utils"
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
	NumberTheoreticTransformer

	// Polynomial nb.Coefficients
	N int

	// Moduli
	Modulus []uint64

	// 2^bit_length(Qi) - 1
	Mask []uint64

	// Indicates whether NTT can be used with the current ring.
	AllowsNTT bool

	// Product of the Moduli for each level
	ModulusAtLevel []*big.Int

	// Fast reduction parameters
	BredParams [][]uint64
	MredParams []uint64

	RescaleParams [][]uint64

	//NTT Parameters
	NthRoot    uint64
	PsiMont    []uint64 //2N-th primitive root in Montgomery form
	PsiInvMont []uint64 //2N-th inverse primitive root in Montgomery form

	NttPsi    [][]uint64 //powers of the inverse of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	NttPsiInv [][]uint64 //powers of the inverse of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	NttNInv   []uint64   //[N^-1] mod Qi in Montgomery form
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
	err = r.setParameters(N, Moduli)
	if err != nil {
		return nil, err
	}

	r.NumberTheoreticTransformer = ntt

	err = r.genNTTParams(uint64(NthRoot))
	if err != nil {
		return r, err
	}

	return r, nil
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
	cr.N = r.N >> 1
	cr.NumberTheoreticTransformer = NumberTheoreticTransformerConjugateInvariant{}
	return &cr, cr.genNTTParams(uint64(cr.N) << 2)
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
	sr.N = r.N << 1
	sr.NumberTheoreticTransformer = NumberTheoreticTransformerStandard{}
	return &sr, sr.genNTTParams(uint64(sr.N) << 1)
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

// setParameters initializes a *Ring by setting the required pre-computed values (except for the NTT-related values, which are set by the
// genNTTParams function).
func (r *Ring) setParameters(N int, Modulus []uint64) error {

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

	r.AllowsNTT = false

	r.N = N

	r.Modulus = make([]uint64, len(Modulus))
	r.Mask = make([]uint64, len(Modulus))

	for i, qi := range Modulus {
		r.Modulus[i] = qi
		r.Mask[i] = (1 << uint64(bits.Len64(qi))) - 1
	}

	// Computes bigQ for all levels
	r.ModulusAtLevel = make([]*big.Int, len(r.Modulus))
	r.ModulusAtLevel[0] = NewUint(r.Modulus[0])
	for i := 1; i < len(r.Modulus); i++ {
		r.ModulusAtLevel[i] = new(big.Int).Mul(r.ModulusAtLevel[i-1], NewUint(r.Modulus[i]))
	}

	// Computes the fast reduction parameters
	r.BredParams = make([][]uint64, len(r.Modulus))
	r.MredParams = make([]uint64, len(r.Modulus))

	for i, qi := range r.Modulus {

		// Computes the fast modular reduction parameters for the Ring
		r.BredParams[i] = BRedParams(qi)

		// If qi is not a power of 2, we can compute the MRedParams (otherwise, it
		// would return an error as there is no valid Montgomery form mod a power of 2)
		if (qi&(qi-1)) != 0 && qi != 0 {
			r.MredParams[i] = MRedParams(qi)
		}
	}

	return nil
}

// genNTTParams checks that N has been correctly initialized, and checks that each modulus is a prime congruent to 1 mod 2N (i.e. NTT-friendly).
// Then, it computes the variables required for the NTT. The purpose of ValidateParameters is to validate that the moduli allow the NTT, and to compute the
// NTT parameters.
func (r *Ring) genNTTParams(NthRoot uint64) error {

	if r.N == 0 || r.Modulus == nil {
		return errors.New("invalid r parameters (missing)")
	}

	if r.N == 0 || r.Modulus == nil || NthRoot < 1 {
		panic("error : invalid r parameters (missing)")
	}

	// Checks if each qi is prime and equal to 1 mod NthRoot
	for i, qi := range r.Modulus {
		if !IsPrime(qi) {
			return fmt.Errorf("invalid modulus (Modulus[%d] is not prime)", i)
		}

		if qi&(NthRoot-1) != 1 {
			r.AllowsNTT = false
			return fmt.Errorf("invalid modulus (Modulus[%d] != 1 mod NthRoot)", i)
		}
	}

	r.NthRoot = NthRoot

	r.RescaleParams = make([][]uint64, len(r.Modulus)-1)

	for j := len(r.Modulus) - 1; j > 0; j-- {

		r.RescaleParams[j-1] = make([]uint64, j)

		for i := 0; i < j; i++ {

			r.RescaleParams[j-1][i] = MForm(r.Modulus[i]-ModExp(r.Modulus[j], r.Modulus[i]-2, r.Modulus[i]), r.Modulus[i], r.BredParams[i])
		}
	}

	r.PsiMont = make([]uint64, len(r.Modulus))
	r.PsiInvMont = make([]uint64, len(r.Modulus))
	r.NttPsi = make([][]uint64, len(r.Modulus))
	r.NttPsiInv = make([][]uint64, len(r.Modulus))
	r.NttNInv = make([]uint64, len(r.Modulus))

	logNthRoot := uint64(bits.Len64(NthRoot>>1) - 1)

	for i, qi := range r.Modulus {

		// 1.1 Computes N^(-1) mod Q in Montgomery form
		r.NttNInv[i] = MForm(ModExp(NthRoot>>1, qi-2, qi), qi, r.BredParams[i])

		// 1.2 Computes Psi and PsiInv in Montgomery form
		r.NttPsi[i] = make([]uint64, NthRoot>>1)
		r.NttPsiInv[i] = make([]uint64, NthRoot>>1)

		// Finds a 2N-th primitive Root
		g := primitiveRoot(qi)

		power := (qi - 1) / NthRoot
		powerInv := (qi - 1) - power

		// Computes Psi and PsiInv in Montgomery form
		PsiMont := MForm(ModExp(g, power, qi), qi, r.BredParams[i])
		PsiInvMont := MForm(ModExp(g, powerInv, qi), qi, r.BredParams[i])

		r.PsiMont[i] = PsiMont
		r.PsiInvMont[i] = PsiInvMont

		r.NttPsi[i][0] = MForm(1, qi, r.BredParams[i])
		r.NttPsiInv[i][0] = MForm(1, qi, r.BredParams[i])

		// Computes nttPsi[j] = nttPsi[j-1]*Psi and nttPsiInv[j] = nttPsiInv[j-1]*PsiInv
		for j := uint64(1); j < NthRoot>>1; j++ {

			indexReversePrev := utils.BitReverse64(uint64(j-1), logNthRoot)
			indexReverseNext := utils.BitReverse64(uint64(j), logNthRoot)

			r.NttPsi[i][indexReverseNext] = MRed(r.NttPsi[i][indexReversePrev], PsiMont, qi, r.MredParams[i])
			r.NttPsiInv[i][indexReverseNext] = MRed(r.NttPsiInv[i][indexReversePrev], PsiInvMont, qi, r.MredParams[i])
		}
	}

	r.AllowsNTT = true

	return nil
}

// Minimal required information to recover the full ring. Used to import and export the ring.
type ringParams struct {
	N       int
	NthRoot uint64
	Modulus []uint64
}

// MarshalBinary encodes the target ring on a slice of bytes.
func (r *Ring) MarshalBinary() ([]byte, error) {

	parameters := ringParams{r.N, r.NthRoot, r.Modulus}

	var buf bytes.Buffer
	enc := gob.NewEncoder(&buf)
	if err := enc.Encode(parameters); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// UnmarshalBinary decodes a slice of bytes on the target Ring.
func (r *Ring) UnmarshalBinary(data []byte) error {

	parameters := ringParams{}

	reader := bytes.NewReader(data)
	dec := gob.NewDecoder(reader)
	if err := dec.Decode(&parameters); err != nil {
		return err
	}

	if err := r.setParameters(parameters.N, parameters.Modulus); err != nil {
		return err
	}
	if err := r.genNTTParams(parameters.NthRoot); err != nil {
		return err
	}

	return nil
}

// NewPoly creates a new polynomial with all coefficients set to 0.
func (r *Ring) NewPoly() *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, len(r.Modulus))
	for i := 0; i < len(r.Modulus); i++ {
		p.Coeffs[i] = make([]uint64, r.N)
	}

	return p
}

// NewPolyLvl creates a new polynomial with all coefficients set to 0.
func (r *Ring) NewPolyLvl(level int) *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, level+1)
	for i := 0; i < level+1; i++ {
		p.Coeffs[i] = make([]uint64, r.N)
	}

	return p
}

// SetCoefficientsInt64 sets the coefficients of p1 from an int64 array.
func (r *Ring) SetCoefficientsInt64(coeffs []int64, p1 *Poly) {
	for i, coeff := range coeffs {
		for j, Qi := range r.Modulus {
			p1.Coeffs[j][i] = CRed(uint64((coeff%int64(Qi) + int64(Qi))), Qi)
		}
	}
}

// SetCoefficientsUint64 sets the coefficients of p1 from an uint64 array.
func (r *Ring) SetCoefficientsUint64(coeffs []uint64, p1 *Poly) {
	for i, coeff := range coeffs {
		for j, Qi := range r.Modulus {
			p1.Coeffs[j][i] = coeff % Qi
		}
	}
}

// SetCoefficientsString parses an array of string as Int variables, and sets the
// coefficients of p1 with these Int variables.
func (r *Ring) SetCoefficientsString(coeffs []string, p1 *Poly) {
	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, Qi := range r.Modulus {
		QiBigint.SetUint64(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(NewIntFromString(coeff), QiBigint).Uint64()
		}
	}
}

// SetCoefficientsBigint sets the coefficients of p1 from an array of Int variables.
func (r *Ring) SetCoefficientsBigint(coeffs []*big.Int, p1 *Poly) {
	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, Qi := range r.Modulus {
		QiBigint.SetUint64(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}
}

// SetCoefficientsBigintLvl sets the coefficients of p1 from an array of Int variables.
func (r *Ring) SetCoefficientsBigintLvl(level int, coeffs []*big.Int, p1 *Poly) {

	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i := 0; i < level+1; i++ {
		QiBigint.SetUint64(r.Modulus[i])
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}
}

// PolyToString reconstructs p1 and returns the result in an array of string.
func (r *Ring) PolyToString(p1 *Poly) []string {

	coeffsBigint := make([]*big.Int, r.N)
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
	r.PolyToBigintLvl(p1.Level(), p1, gap, coeffsBigint)
}

// PolyToBigintLvl reconstructs p1 and returns the result in an array of Int.
// gap defines coefficients X^{i*gap} that will be reconstructed.
// For example, if gap = 1, then all coefficients are reconstructed, while
// if gap = 2 then only coefficients X^{2*i} are reconstructed.
func (r *Ring) PolyToBigintLvl(level int, p1 *Poly, gap int, coeffsBigint []*big.Int) {

	crtReconstruction := make([]*big.Int, level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := r.ModulusAtLevel[level]

	for i := 0; i < level+1; i++ {
		QiB.SetUint64(r.Modulus[i])
		crtReconstruction[i] = new(big.Int).Quo(modulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	for i, j := 0, 0; j < r.N; i, j = i+1, j+gap {

		tmp.SetUint64(0)
		coeffsBigint[i] = new(big.Int)

		for k := 0; k < level+1; k++ {
			coeffsBigint[i].Add(coeffsBigint[i], tmp.Mul(NewUint(p1.Coeffs[k][j]), crtReconstruction[k]))
		}

		coeffsBigint[i].Mod(coeffsBigint[i], modulusBigint)
	}
}

// PolyToBigintCenteredLvl reconstructs p1 and returns the result in an array of Int.
// Coefficients are centered around Q/2
// gap defines coefficients X^{i*gap} that will be reconstructed.
// For example, if gap = 1, then all coefficients are reconstructed, while
// if gap = 2 then only coefficients X^{2*i} are reconstructed.
func (r *Ring) PolyToBigintCenteredLvl(level int, p1 *Poly, gap int, coeffsBigint []*big.Int) {

	crtReconstruction := make([]*big.Int, level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := r.ModulusAtLevel[level]

	for i := 0; i < level+1; i++ {
		QiB.SetUint64(r.Modulus[i])
		crtReconstruction[i] = new(big.Int).Quo(modulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	modulusBigintHalf := new(big.Int)
	modulusBigintHalf.Rsh(modulusBigint, 1)

	var sign int
	for i, j := 0, 0; j < r.N; i, j = i+1, j+gap {

		tmp.SetUint64(0)
		coeffsBigint[i].SetUint64(0)

		for k := 0; k < level+1; k++ {
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

	for i := 0; i < len(r.Modulus); i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	r.Reduce(p1, p1)
	r.Reduce(p2, p2)

	for i := 0; i < len(r.Modulus); i++ {
		for j := 0; j < r.N; j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}

// EqualLvl checks if p1 = p2 in the given Ring, up to a given level.
func (r *Ring) EqualLvl(level int, p1, p2 *Poly) bool {

	for i := 0; i < level+1; i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	r.ReduceLvl(level, p1, p1)
	r.ReduceLvl(level, p2, p2)

	for i := 0; i < level+1; i++ {
		for j := 0; j < r.N; j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}
