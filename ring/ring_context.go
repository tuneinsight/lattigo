// Package ring implelents a RNS-accelerated modular arithmetic operations for polynomials, including: RNS basis extension; RNS rescaling;  number theoretic transform (NTT); uniform, Gaussian and ternary sampling.
package ring

import (
	"bytes"
	"crypto/rand"
	"encoding/binary"
	"encoding/gob"
	"errors"
	"math/bits"
)

//==============================
//===== POLYNOMIAL CONTEXT =====
//==============================

// Context is a structure keeping all the variable required to operate on a polynomial represented in this context.
// This include its moduli, crt reconstruction, modular reduction and ntt transformation.
type Context struct {

	// Polynomial nb.Coefficients
	N uint64

	// Modulies
	Modulus []uint64

	// 2^bit_length(Qi) - 1
	mask []uint64

	// Determines if NTT can be used with the current context.
	allowsNTT bool

	// Product of the modulies
	ModulusBigint *Int

	// Parameters for the CRT reconstruction
	CrtReconstruction []*Int

	// Fast reduction parameters
	bredParams [][]uint64
	mredParams []uint64

	//NTT Parameters
	psiMont    []uint64 //2nth primitive root in montgomery form
	psiInvMont []uint64 //2nth inverse primitive root in montgomery form

	nttPsi    [][]uint64 //powers of the inverse of the 2nth primitive root in montgomery form (in bitreversed order)
	nttPsiInv [][]uint64 //powers of the inverse of the 2nth primitive root in montgomery form (in bitreversed order)
	nttNInv   []uint64   //[N^-1] mod Qi in montgomery form
}

// NewContext generates a new empty context.
func NewContext() *Context {
	return new(Context)
}

// TODO : comment
func (context *Context) CopyNew() (contextCopy *Context) {
	contextCopy = new(Context)
	contextCopy.N = context.N
	contextCopy.Modulus = make([]uint64, len(context.Modulus))
	contextCopy.mask = make([]uint64, len(context.Modulus))
	contextCopy.allowsNTT = context.allowsNTT
	contextCopy.ModulusBigint = Copy(context.ModulusBigint)
	contextCopy.CrtReconstruction = make([]*Int, len(context.Modulus))
	contextCopy.bredParams = make([][]uint64, len(context.Modulus))
	contextCopy.mredParams = make([]uint64, len(context.Modulus))
	contextCopy.psiMont = make([]uint64, len(context.Modulus))
	contextCopy.psiInvMont = make([]uint64, len(context.Modulus))
	contextCopy.nttPsi = make([][]uint64, len(context.Modulus))
	contextCopy.nttPsiInv = make([][]uint64, len(context.Modulus))
	contextCopy.nttNInv = make([]uint64, len(context.Modulus))

	for i := range context.Modulus {
		contextCopy.Modulus[i] = context.Modulus[i]
		contextCopy.mask[i] = context.mask[i]
		contextCopy.CrtReconstruction[i] = Copy(context.CrtReconstruction[i])
		contextCopy.bredParams[i] = make([]uint64, 2)
		contextCopy.bredParams[i][0] = context.bredParams[i][0]
		contextCopy.bredParams[i][1] = context.bredParams[i][1]
		contextCopy.mredParams[i] = context.mredParams[i]
		contextCopy.psiMont[i] = context.psiMont[i]
		contextCopy.psiInvMont[i] = context.psiInvMont[i]
		contextCopy.nttPsi[i] = make([]uint64, context.N)
		contextCopy.nttPsiInv[i] = make([]uint64, context.N)
		for j := uint64(0); j < context.N; j++ {
			contextCopy.nttPsi[i][j] = context.nttPsi[i][j]
			contextCopy.nttPsiInv[i][j] = context.nttPsiInv[i][j]
		}
		contextCopy.nttNInv[i] = context.nttNInv[i]
	}

	return
}

// SetParameters initialize the parameters of an empty context with N and the provided moduli.
// Only checks that N is a power of 2 and computes all the variable that aren't used for the NTT.
func (context *Context) SetParameters(N uint64, Modulus []uint64) error {

	// Checks if N is a power of 2
	if (N&(N-1)) != 0 && N != 0 {
		return errors.New("invalid ring degree (must be a power of 2)")
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
		// because it will return an error and there is no valid montgomery form mod a power of 2)
		if (qi&(qi-1)) != 0 && qi != 0 {
			context.mredParams[i] = MRedParams(qi)
		}
	}

	return nil
}

// GenNTTParams checks that N has beed correctly initialized, and checks that each moduli is a prime congruent to 1 mod 2N (i.e. allowing NTT).
// Then it computes the variables required for the NTT. ValidateParameters purpose is to validate that the moduli allow the NTT and compute the
// NTT parameters.
func (context *Context) GenNTTParams() error {

	if context.allowsNTT {
		return nil
	}

	if context.N == 0 || context.Modulus == nil {
		return errors.New("error : invalid context parameters (missing)")
	}

	// CHECKS IF VALIDE NTT
	// Checks if each qi is Prime and if qi = 1 mod 2n
	for _, qi := range context.Modulus {
		if IsPrime(qi) == false || qi&((context.N<<1)-1) != 1 {
			context.allowsNTT = false
			return errors.New("warning : provided modulus does not allow NTT")
		}
	}

	context.CrtReconstruction = make([]*Int, len(context.Modulus))

	context.psiMont = make([]uint64, len(context.Modulus))
	context.psiInvMont = make([]uint64, len(context.Modulus))
	context.nttPsi = make([][]uint64, len(context.Modulus))
	context.nttPsiInv = make([][]uint64, len(context.Modulus))
	context.nttNInv = make([]uint64, len(context.Modulus))

	bitLenofN := uint64(bits.Len64(context.N) - 1)

	QiB := new(Int)
	tmp := new(Int)

	for i, qi := range context.Modulus {

		//1.0 CRT reconstruction parameters
		QiB.SetUint(qi)
		context.CrtReconstruction[i] = new(Int)
		context.CrtReconstruction[i].Div(context.ModulusBigint, QiB)
		tmp.Inv(context.CrtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		context.CrtReconstruction[i].Mul(context.CrtReconstruction[i], tmp)

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

			indexReversePrev := bitReverse64(j-1, bitLenofN)
			indexReverseNext := bitReverse64(j, bitLenofN)

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

// Merge merges two context by appending all the element from contextP to the elements of contextQ
// Will return an error if contextQ or contextP do not both agree on the flat allowsNTT. It
// however requires to re-compute the crt reconstruction parameters.
func (context *Context) Merge(contextQ, contextP *Context) error {

	if contextQ.N != contextP.N {
		return errors.New("contexts ring degree to not match")
	}

	context.N = contextQ.N

	context.Modulus = append(contextQ.Modulus, contextP.Modulus...)
	context.mask = append(contextQ.mask, contextP.mask...)

	if context != contextQ && context != contextP {
		context.ModulusBigint = NewUint(0)
	}

	context.ModulusBigint.Mul(contextQ.ModulusBigint, contextP.ModulusBigint)

	// For this part we need to recompute, since each element is a function of all the other modulus
	context.CrtReconstruction = append(contextQ.CrtReconstruction, contextP.CrtReconstruction...)
	QiB := new(Int)
	tmp := new(Int)
	for i, qi := range context.Modulus {
		QiB.SetUint(qi)
		context.CrtReconstruction[i] = new(Int)
		context.CrtReconstruction[i].Div(context.ModulusBigint, QiB)
		tmp.Inv(context.CrtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		context.CrtReconstruction[i].Mul(context.CrtReconstruction[i], tmp)
	}

	context.bredParams = append(contextQ.bredParams, contextP.bredParams...)
	context.mredParams = append(contextQ.mredParams, contextP.mredParams...)

	context.psiMont = append(contextQ.psiMont, contextP.psiMont...)
	context.psiInvMont = append(contextQ.psiInvMont, contextP.psiInvMont...)

	if contextQ.allowsNTT == false && contextP.allowsNTT == false {

		context.allowsNTT = false

	} else if contextQ.allowsNTT && contextP.allowsNTT {

		context.nttPsi = append(contextQ.nttPsi, contextP.nttPsi...)
		context.nttPsiInv = append(contextQ.nttPsiInv, contextP.nttPsiInv...)
		context.nttNInv = append(contextQ.nttNInv, contextP.nttNInv...)
		context.allowsNTT = true

	} else {

		return errors.New("context need both to be allowsNTT or not allowsNTT")
	}

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

// NewUniformPoly generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1]
func (context *Context) UniformPoly(Pol *Poly) {

	var randomBytes []byte
	var randomUint, mask, qi uint64

	n := context.N
	if n < 8 {
		n = 8
	}

	randomBytes = make([]byte, n)
	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	level := uint64(len(Pol.Coeffs) - 1)

	for j := uint64(0); j < level+1; j++ {

		qi = context.Modulus[j]

		// Starts by computing the mask
		mask = (1 << uint64(bits.Len64(qi))) - 1

		// Iterates for each modulus over each coefficient
		for i := uint64(0); i < context.N; i++ {

			// Samples an integer between [0, qi-1]
			for {

				// Replenishes the pool if it runs empty
				if len(randomBytes) < 8 {
					randomBytes = make([]byte, n)
					if _, err := rand.Read(randomBytes); err != nil {
						panic("crypto rand error")
					}
				}

				// Reads bytes from the pool
				randomUint = binary.BigEndian.Uint64(randomBytes[:8]) & mask
				randomBytes = randomBytes[8:] // Discard the used bytes

				// If the integer is between [0, qi-1], breaks the loop
				if randomUint < qi {
					break
				}
			}

			Pol.Coeffs[j][i] = randomUint
		}
	}

	return
}

// NewUniformPoly generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1]
func (context *Context) NewUniformPoly() (Pol *Poly) {

	Pol = context.NewPoly()

	context.UniformPoly(Pol)

	return
}

func (context *Context) NewUniformPolyLvl(level uint64) (Pol *Poly) {

	Pol = context.NewPolyLvl(level)

	context.UniformPoly(Pol)

	return
}

// SetCoefficientsInt64 sets the coefficients of p1 from an int64 array.
func (context *Context) SetCoefficientsInt64(coeffs []int64, p1 *Poly) error {
	if len(coeffs) != int(context.N) {
		return errors.New("error : invalid ring degree (does not match context)")
	}
	for i, coeff := range coeffs {
		for j, Qi := range context.Modulus {
			p1.Coeffs[j][i] = CRed(uint64((coeff%int64(Qi) + int64(Qi))), Qi)
		}
	}
	return nil
}

// SetCoefficientsUint64 sets the coefficient of Pol from an uint64 array.
func (context *Context) SetCoefficientsUint64(coeffs []uint64, p1 *Poly) error {
	if len(coeffs) != int(context.N) {
		return errors.New("error : invalid ring degree (does not match context)")
	}
	for i, coeff := range coeffs {
		for j, Qi := range context.Modulus {
			p1.Coeffs[j][i] = coeff % Qi
		}
	}
	return nil
}

// SetCoefficientsString parses an array of string as Int variables, and sets the
// coefficients of p1 with this Int variables.
func (context *Context) SetCoefficientsString(coeffs []string, p1 *Poly) error {

	if len(coeffs) != int(context.N) {
		return errors.New("error : invalid ring degree (does not match context)")
	}

	QiBigint := new(Int)
	coeffTmp := new(Int)
	for i, Qi := range context.Modulus {
		QiBigint.SetUint(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(NewIntFromString(coeff), QiBigint).Uint64()
		}
	}
	return nil
}

// SetCoefficientsBigint sets the coefficients of p1 from an array of Int variables.
func (context *Context) SetCoefficientsBigint(coeffs []*Int, p1 *Poly) error {

	if len(coeffs) != int(context.N) {
		return errors.New("error : invalid ring degree (does not match context)")
	}

	QiBigint := new(Int)
	coeffTmp := new(Int)
	for i, Qi := range context.Modulus {
		QiBigint.SetUint(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}

	return nil
}

func (context *Context) SetCoefficientsBigintLvl(level uint64, coeffs []*Int, p1 *Poly) error {

	if len(coeffs) != int(context.N) {
		return errors.New("error : invalid ring degree (does not match context)")
	}

	QiBigint := new(Int)
	coeffTmp := new(Int)
	for i := uint64(0); i < level+1; i++ {
		QiBigint.SetUint(context.Modulus[i])
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}

	return nil
}

//PolyToString reconstructs p1 and returns the result in an array of string.
func (context *Context) PolyToString(p1 *Poly) []string {

	coeffsBigint := make([]*Int, context.N)
	context.PolyToBigint(p1, coeffsBigint)
	coeffsString := make([]string, len(coeffsBigint))

	for i := range coeffsBigint {
		coeffsString[i] = coeffsBigint[i].String()
	}

	return coeffsString
}

//PolyToBigint reconstructs p1 and returns the result in an array of Int.
func (context *Context) PolyToBigint(p1 *Poly, coeffsBigint []*Int) {

	var qi, level uint64

	level = uint64(len(p1.Coeffs) - 1)

	crtReconstruction := make([]*Int, level+1)

	QiB := new(Int)
	tmp := new(Int)
	modulusBigint := NewUint(1)

	for i := uint64(0); i < level+1; i++ {

		qi = context.Modulus[i]
		QiB.SetUint(qi)

		modulusBigint.Mul(modulusBigint, QiB)

		crtReconstruction[i] = new(Int)
		crtReconstruction[i].Div(context.ModulusBigint, QiB)
		tmp.Inv(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	for x := uint64(0); x < context.N; x++ {

		tmp.SetUint(0)
		coeffsBigint[x] = new(Int)

		for i := uint64(0); i < level+1; i++ {
			coeffsBigint[x].Add(coeffsBigint[x], tmp.Mul(NewUint(p1.Coeffs[i][x]), crtReconstruction[i]))
		}

		coeffsBigint[x].Mod(coeffsBigint[x], modulusBigint)
	}
}

// GetCenteredCoefficients returns an array containing the coefficients of p1 centered arount each (-Qi/2, Qi/2].
func (context *Context) GetCenteredCoefficients(p1 *Poly) [][]int64 {

	coeffs := make([][]int64, len(context.Modulus))
	var qiHalf int64

	for i, qi := range context.Modulus {
		qiHalf = int64(qi >> 1)
		coeffs[i] = make([]int64, context.N)

		for j := uint64(0); j < context.N; j++ {
			coeffs[i][j] = int64(p1.Coeffs[i][j])

			if coeffs[i][j] > qiHalf {
				coeffs[i][j] -= int64(qi)
			}
		}
	}

	return coeffs
}

// Equal checks if p1 = p2 in the given context.
func (context *Context) Equal(p1, p2 *Poly) bool {

	//if len(p1.Coeffs) != len(context.Modulus) || len(p2.Coeffs) != len(context.Modulus){
	//	return false
	//}

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
