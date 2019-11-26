// Package ring implelents a RNS-accelerated modular arithmetic operations for polynomials, including: RNS basis extension; RNS rescaling;  number theoretic transform (NTT); uniform, Gaussian and ternary sampling.
package ring

import (
	"bytes"
	"encoding/gob"
	"errors"
	"github.com/ldsec/lattigo/utils"
	"math/big"
	"math/bits"
)

//==============================
//===== POLYNOMIAL CONTEXT =====
//==============================

// Context is a structure keeping all the variables required to operate on a polynomial represented in this context.
type Context struct {

	// Polynomial nb.Coefficients
	N uint64

	// Moduli
	Modulus []uint64

	// 2^bit_length(Qi) - 1
	mask []uint64

	// Determines if NTT can be used with the current context.
	allowsNTT bool

	// Product of the Moduli
	ModulusBigint *big.Int

	// Fast reduction parameters
	bredParams [][]uint64
	mredParams []uint64

	rescaleParams [][]uint64

	matrixTernary           [][]uint64
	matrixTernaryMontgomery [][]uint64

	//NTT Parameters
	psiMont    []uint64 //2nth primitive root in Montgomery form
	psiInvMont []uint64 //2nth inverse primitive root in Montgomery form

	nttPsi    [][]uint64 //powers of the inverse of the 2nth primitive root in Montgomery form (in bitreversed order)
	nttPsiInv [][]uint64 //powers of the inverse of the 2nth primitive root in Montgomery form (in bitreversed order)
	nttNInv   []uint64   //[N^-1] mod Qi in Montgomery form
}

// NewContext generates a new empty context.
func NewContext() *Context {
	return new(Context)
}

// SetParameters initializes the parameters of an empty context with N and the provided moduli.
// Only checks that N is a power of 2 and computes all the variables that aren't used for the NTT.
func (context *Context) SetParameters(N uint64, Modulus []uint64) {

	// Checks if N is a power of 2
	if (N&(N-1)) != 0 && N != 0 {
		panic("invalid ring degree (must be a power of 2)")
	}

	context.allowsNTT = false

	context.N = N

	context.Modulus = make([]uint64, len(Modulus))
	context.mask = make([]uint64, len(Modulus))

	for i, qi := range Modulus {
		context.Modulus[i] = qi
		context.mask[i] = (1 << uint64(bits.Len64(qi))) - 1
	}

	//Computes the bigQ
	context.ModulusBigint = NewInt(1)
	for _, qi := range context.Modulus {
		context.ModulusBigint.Mul(context.ModulusBigint, NewUint(qi))
	}

	// Computes the fast reduction parameters
	context.bredParams = make([][]uint64, len(context.Modulus))
	context.mredParams = make([]uint64, len(context.Modulus))

	for i, qi := range context.Modulus {

		//Computes the fast modular reduction parameters for the Context
		context.bredParams[i] = BRedParams(qi)

		// If qi is not a power of 2, we can compute the MRedParams (else it should not
		// because it will return an error and there is no valid Montgomery form mod a power of 2)
		if (qi&(qi-1)) != 0 && qi != 0 {
			context.mredParams[i] = MRedParams(qi)
		}
	}

	context.matrixTernary = make([][]uint64, len(context.Modulus))
	context.matrixTernaryMontgomery = make([][]uint64, len(context.Modulus))

	for i, Qi := range context.Modulus {

		context.matrixTernary[i] = make([]uint64, 3)
		context.matrixTernary[i][0] = 0
		context.matrixTernary[i][1] = 1
		context.matrixTernary[i][2] = Qi - 1

		context.matrixTernaryMontgomery[i] = make([]uint64, 3)
		context.matrixTernaryMontgomery[i][0] = 0
		context.matrixTernaryMontgomery[i][1] = MForm(1, Qi, context.bredParams[i])
		context.matrixTernaryMontgomery[i][2] = MForm(Qi-1, Qi, context.bredParams[i])
	}
}

// GenNTTParams checks that N has been correctly initialized, and checks that each moduli is a prime congruent to 1 mod 2N (i.e. allowing NTT).
// Then it computes the variables required for the NTT. ValidateParameters purpose is to validate that the moduli allow the NTT and compute the
// NTT parameters.
func (context *Context) GenNTTParams() error {

	if context.allowsNTT {
		return nil
	}

	if context.N == 0 || context.Modulus == nil {
		panic("error : invalid context parameters (missing)")
	}

	// CHECKS IF VALIDE NTT
	// Checks if each qi is Prime and if qi = 1 mod 2n
	for _, qi := range context.Modulus {
		if IsPrime(qi) == false || qi&((context.N<<1)-1) != 1 {
			context.allowsNTT = false
			return errors.New("warning : provided modulus does not allow NTT")
		}
	}

	context.rescaleParams = make([][]uint64, len(context.Modulus)-1)

	for j := len(context.Modulus) - 1; j > 0; j-- {

		context.rescaleParams[j-1] = make([]uint64, j)

		for i := 0; i < j; i++ {

			context.rescaleParams[j-1][i] = MForm(ModExp(context.Modulus[j], context.Modulus[i]-2, context.Modulus[i]), context.Modulus[i], context.bredParams[i])
		}
	}

	context.psiMont = make([]uint64, len(context.Modulus))
	context.psiInvMont = make([]uint64, len(context.Modulus))
	context.nttPsi = make([][]uint64, len(context.Modulus))
	context.nttPsiInv = make([][]uint64, len(context.Modulus))
	context.nttNInv = make([]uint64, len(context.Modulus))

	bitLenofN := uint64(bits.Len64(context.N) - 1)

	for i, qi := range context.Modulus {

		//2.1 Computes N^(-1) mod Q in Montgomery form
		context.nttNInv[i] = MForm(ModExp(context.N, qi-2, qi), qi, context.bredParams[i])

		//2.2 Computes Psi and PsiInv in Montgomery form
		context.nttPsi[i] = make([]uint64, context.N)
		context.nttPsiInv[i] = make([]uint64, context.N)

		//Finds a 2nth primitive Root
		g := primitiveRoot(qi)

		_2n := uint64(context.N << 1)

		power := (qi - 1) / _2n
		powerInv := (qi - 1) - power

		//Computes Psi and PsiInv in Montgomery Form
		PsiMont := MForm(ModExp(g, power, qi), qi, context.bredParams[i])
		PsiInvMont := MForm(ModExp(g, powerInv, qi), qi, context.bredParams[i])

		context.psiMont[i] = PsiMont
		context.psiInvMont[i] = PsiInvMont

		context.nttPsi[i][0] = MForm(1, qi, context.bredParams[i])
		context.nttPsiInv[i][0] = MForm(1, qi, context.bredParams[i])

		// Computes nttPsi[j] = nttPsi[j-1]*Psi and nttPsiInv[j] = nttPsiInv[j-1]*PsiInv
		for j := uint64(1); j < context.N; j++ {

			indexReversePrev := utils.BitReverse64(j-1, bitLenofN)
			indexReverseNext := utils.BitReverse64(j, bitLenofN)

			context.nttPsi[i][indexReverseNext] = MRed(context.nttPsi[i][indexReversePrev], PsiMont, qi, context.mredParams[i])
			context.nttPsiInv[i][indexReverseNext] = MRed(context.nttPsiInv[i][indexReversePrev], PsiInvMont, qi, context.mredParams[i])
		}
	}

	context.allowsNTT = true

	return nil
}

// Used to export the context. Minimal information to recover the full context.
type smallContext struct {
	N       uint64
	Modulus []uint64
}

func (context *Context) MarshalBinary() ([]byte, error) {

	parameters := smallContext{context.N, context.Modulus}

	var buf bytes.Buffer
	enc := gob.NewEncoder(&buf)
	if err := enc.Encode(parameters); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

func (context *Context) UnMarshalBinary(data []byte) error {

	parameters := smallContext{}

	reader := bytes.NewReader(data)
	dec := gob.NewDecoder(reader)
	if err := dec.Decode(&parameters); err != nil {
		return err
	}

	context.SetParameters(parameters.N, parameters.Modulus)
	context.GenNTTParams()

	return nil
}

// AllowsNTT returns true if the context allows NTT, else false.
func (context *Context) AllowsNTT() bool {
	return context.allowsNTT
}

// GetBRedParams returns the Barret reduction parameters of the context.
func (context *Context) GetBredParams() [][]uint64 {
	return context.bredParams
}

// GetMredParams returns the Montgomery reduction parameters of the context.
func (context *Context) GetMredParams() []uint64 {
	return context.mredParams
}

// GetPsi returns the primitive root used to compute the NTT parameters of the context.
func (context *Context) GetPsi() []uint64 {
	return context.psiMont
}

// GetPsi returns the primitive root used to compute the InvNTT parameters of the context.
func (context *Context) GetPsiInv() []uint64 {
	return context.psiInvMont
}

// GetNttPsi returns the NTT parameters of the context.
func (context *Context) GetNttPsi() [][]uint64 {
	return context.nttPsi
}

//GetNttPsiInv returns the InvNTT parameters of the context.
func (context *Context) GetNttPsiInv() [][]uint64 {
	return context.nttPsiInv
}

// GetNttNInv returns 1/N mod each moduli.
func (context *Context) GetNttNInv() []uint64 {
	return context.nttNInv
}

// NewPoly create a new polynomial with all coefficients set to 0.
func (context *Context) NewPoly() *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, len(context.Modulus))
	for i := 0; i < len(context.Modulus); i++ {
		p.Coeffs[i] = make([]uint64, context.N)
	}

	return p
}

func (context *Context) NewPolyLvl(level uint64) *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, level+1)
	for i := uint64(0); i < level+1; i++ {
		p.Coeffs[i] = make([]uint64, context.N)
	}

	return p
}

// SetCoefficientsInt64 sets the coefficients of p1 from an int64 array.
func (context *Context) SetCoefficientsInt64(coeffs []int64, p1 *Poly) {
	for i, coeff := range coeffs {
		for j, Qi := range context.Modulus {
			p1.Coeffs[j][i] = CRed(uint64((coeff%int64(Qi) + int64(Qi))), Qi)
		}
	}
}

// SetCoefficientsUint64 sets the coefficients of p1 from an uint64 array.
func (context *Context) SetCoefficientsUint64(coeffs []uint64, p1 *Poly) {
	for i, coeff := range coeffs {
		for j, Qi := range context.Modulus {
			p1.Coeffs[j][i] = coeff % Qi
		}
	}
}

// SetCoefficientsString parses an array of string as Int variables, and sets the
// coefficients of p1 with this Int variables.
func (context *Context) SetCoefficientsString(coeffs []string, p1 *Poly) {
	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, Qi := range context.Modulus {
		QiBigint.SetUint64(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(NewIntFromString(coeff), QiBigint).Uint64()
		}
	}
}

// SetCoefficientsBigint sets the coefficients of p1 from an array of Int variables.
func (context *Context) SetCoefficientsBigint(coeffs []*big.Int, p1 *Poly) {
	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, Qi := range context.Modulus {
		QiBigint.SetUint64(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}
}

func (context *Context) SetCoefficientsBigintLvl(level uint64, coeffs []*big.Int, p1 *Poly) {

	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i := uint64(0); i < level+1; i++ {
		QiBigint.SetUint64(context.Modulus[i])
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}
}

//PolyToString reconstructs p1 and returns the result in an array of string.
func (context *Context) PolyToString(p1 *Poly) []string {

	coeffsBigint := make([]*big.Int, context.N)
	context.PolyToBigint(p1, coeffsBigint)
	coeffsString := make([]string, len(coeffsBigint))

	for i := range coeffsBigint {
		coeffsString[i] = coeffsBigint[i].String()
	}

	return coeffsString
}

//PolyToBigint reconstructs p1 and returns the result in an array of Int.
func (context *Context) PolyToBigint(p1 *Poly, coeffsBigint []*big.Int) {

	var qi, level uint64

	level = uint64(len(p1.Coeffs) - 1)

	crtReconstruction := make([]*big.Int, level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := NewUint(1)

	for i := uint64(0); i < level+1; i++ {

		qi = context.Modulus[i]
		QiB.SetUint64(qi)

		modulusBigint.Mul(modulusBigint, QiB)

		crtReconstruction[i] = new(big.Int)
		crtReconstruction[i].Quo(context.ModulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	for x := uint64(0); x < context.N; x++ {

		tmp.SetUint64(0)
		coeffsBigint[x] = new(big.Int)

		for i := uint64(0); i < level+1; i++ {
			coeffsBigint[x].Add(coeffsBigint[x], tmp.Mul(NewUint(p1.Coeffs[i][x]), crtReconstruction[i]))
		}

		coeffsBigint[x].Mod(coeffsBigint[x], modulusBigint)
	}
}

// Equal checks if p1 = p2 in the given context.
func (context *Context) Equal(p1, p2 *Poly) bool {

	for i := 0; i < len(context.Modulus); i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	context.Reduce(p1, p1)
	context.Reduce(p2, p2)

	for i := 0; i < len(context.Modulus); i++ {
		for j := uint64(0); j < context.N; j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}

func (context *Context) EqualLvl(level uint64, p1, p2 *Poly) bool {

	for i := uint64(0); i < level+1; i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	context.ReduceLvl(level, p1, p1)
	context.ReduceLvl(level, p2, p2)

	for i := uint64(0); i < level+1; i++ {
		for j := uint64(0); j < context.N; j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}
