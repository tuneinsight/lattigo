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

type Context struct {

	// Polynomial nb.Coefficients
	N uint64

	// Modulies
	Modulus []uint64

	// 2^bit_length(Qi) - 1
	mask []uint64

	// Determines if NTT can be used with the current context.
	validated bool

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

func NewContext() *Context {
	return new(Context)
}

func (context *Context) SetParameters(N uint64, Modulus []uint64) error {

	// Checks if N is a power of 2
	if (N&(N-1)) != 0 && N != 0 {
		return errors.New("invalid ring degree (must be a power of 2)")
	}

	context.validated = false

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

func (context *Context) ValidateParameters() error {

	if context.validated {
		return nil
	}

	if context.N == 0 || context.Modulus == nil {
		return errors.New("error : invalid context parameters (missing)")
	}

	// CHECKS IF VALIDE NTT
	// Checks if each qi is Prime and if qi = 1 mod 2n
	for _, qi := range context.Modulus {
		if IsPrime(qi) == false || qi&((context.N<<1)-1) != 1 {
			context.validated = false
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
		context.nttNInv[i] = MForm(modexp(context.N, qi-2, qi), qi, context.bredParams[i])

		//2.2 Computes Psi and PsiInv in Montgomery form
		context.nttPsi[i] = make([]uint64, context.N)
		context.nttPsiInv[i] = make([]uint64, context.N)

		//Finds a 2nth primitive Root
		g := primitiveRoot(qi)

		_2n := uint64(context.N << 1)

		power := (qi - 1) / _2n
		powerInv := (qi - 1) - power

		//Computes Psi and PsiInv in Montgomery Form
		PsiMont := MForm(modexp(g, power, qi), qi, context.bredParams[i])
		PsiInvMont := MForm(modexp(g, powerInv, qi), qi, context.bredParams[i])

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

	context.validated = true

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
	context.ValidateParameters()

	return nil
}

// Merge merges two context by appending all the element from contextP to the elements of contextQ
// Will return an error if contextQ or contextP are not both not validated or both validated
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

	if contextQ.validated == false && contextP.validated == false {

		context.validated = false

	} else if contextQ.validated && contextP.validated {

		context.nttPsi = append(contextQ.nttPsi, contextP.nttPsi...)
		context.nttPsiInv = append(contextQ.nttPsiInv, contextP.nttPsiInv...)
		context.nttNInv = append(contextQ.nttNInv, contextP.nttNInv...)
		context.validated = true

	} else {

		return errors.New("context need both to be validated or not validated")
	}

	return nil
}

func (context *Context) IsValidated() bool {
	return context.validated
}

func (context *Context) GetBredParams() [][]uint64 {
	return context.bredParams
}

func (context *Context) GetMredParams() []uint64 {
	return context.mredParams
}

func (context *Context) GetPsi() []uint64 {
	return context.psiMont
}

func (context *Context) GetPsiInv() []uint64 {
	return context.psiInvMont
}

func (context *Context) GetNttPsi() [][]uint64 {
	return context.nttPsi
}

func (context *Context) GetNttPsiInv() [][]uint64 {
	return context.nttPsiInv
}

func (context *Context) GetNttNInv() []uint64 {
	return context.nttNInv
}

// Create a new ring with all coefficients set to 0 from the given context
func (context *Context) NewPoly() *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, len(context.Modulus))
	for i := 0; i < len(context.Modulus); i++ {
		p.Coeffs[i] = make([]uint64, context.N)
	}

	return p
}

// Generates a ring with coefficients following a
// discrete gaussian distribution with derivation sigma
func (context *Context) NewGaussPoly(sigma float64) *Poly {

	Pol := context.NewPoly()

	var coeff int64

	boundMultiplier := float64(6)                   // this parameter is the same as the SEAL library
	positiveBound := int64(boundMultiplier * sigma) // the suggested sigma from SEAL is 3.19
	negativeBound := -int64(boundMultiplier * sigma)

	for i := uint64(0); i < context.N; i++ {

		for {

			coeff = GaussSampling(sigma)

			// Samples values until one falls in the desired range
			if (coeff > positiveBound) || (coeff < negativeBound) {
				continue
			}

			// Once the value falls in the desired range, assign it to the coeff, breaks and jumps to next coeff
			for j, qi := range context.Modulus {
				Pol.Coeffs[j][i] = CRed(uint64(int64(qi)+coeff), qi)
			}

			break
		}
	}

	return Pol
}

// Generates a ring with coefficients following a
// uniform distribution over [0, Qi-1]
func (context *Context) NewUniformPoly() (Pol *Poly) {

	var randomBytes []byte
	var randomUint, mask uint64

	Pol = context.NewPoly()

	n := context.N
	if n < 8 {
		n = 8
	}

	randomBytes = make([]byte, n)
	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for j, qi := range context.Modulus {

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

// Generates a ring with coefficients following a
// uniform distribution over [-1,0,1] mod each Qi
func (context *Context) NewTernaryPoly() *Poly {

	Pol := context.NewPoly()

	var coeff uint64
	var sign uint64

	for i := uint64(0); i < context.N; i++ {

		coeff, sign = randInt3()

		for j, qi := range context.Modulus {
			Pol.Coeffs[j][i] = (coeff & (sign * 0xFFFFFFFFFFFFFFFF)) | (((qi * coeff) - coeff) & ((sign ^ 1) * 0xFFFFFFFFFFFFFFFF))
		}
	}

	return Pol
}

// Generates a ring with coefficients following a
// uniform distribution over [-1,0,1] mod each Qi and applies
// the NTT.
func (context *Context) NewTernaryPolyNTT() *Poly {

	Pol := context.NewPoly()

	var coeff uint64
	var sign uint64

	for i := uint64(0); i < context.N; i++ {

		coeff, sign = randInt3()

		for j, qi := range context.Modulus {

			Pol.Coeffs[j][i] = (coeff & (sign * 0xFFFFFFFFFFFFFFFF)) | ((qi * coeff) & ((sign ^ 1) * 0xFFFFFFFFFFFFFFFF))
		}
	}

	context.NTT(Pol, Pol)

	return Pol
}

// Generates a ring with coefficients following a
// uniform distribution over [-1,0,1] mod each Qi and applies
// the NTT and Montgomery form.
func (context *Context) NewTernaryPolyMontgomeryNTT() *Poly {

	Pol := context.NewPoly()

	var coeff uint64
	var sign uint64

	for i := uint64(0); i < context.N; i++ {

		coeff, sign = randInt3()

		for j, qi := range context.Modulus {

			Pol.Coeffs[j][i] = (coeff & (sign * 0xFFFFFFFFFFFFFFFF)) | ((qi - coeff) & ((sign ^ 1) * 0xFFFFFFFFFFFFFFFF))
		}
	}

	context.MForm(Pol, Pol)
	context.NTT(Pol, Pol)

	return Pol
}

// Sets the coefficients of Pol from an int64 array
// If a coefficient is negative it will instead assay its
// positive inverse mod each Qi
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

// Sets the coefficient of Pol from an uint64 array
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

// Sets the coefficients of Pol from a array of strings
// It will import them as bigint variables and reduce them mod each Qi.
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

// Sets the coefficients of Pol from a array of strings
// It will import them as bigint variables and
// reduce them mod each Qi
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

//Returns an string array containing the reconstructed coefficients
func (context *Context) PolyToString(p1 *Poly) []string {

	coeffsBigint := make([]*Int, context.N)
	context.PolyToBigint(p1, coeffsBigint)
	coeffsString := make([]string, len(coeffsBigint))

	for i := range coeffsBigint {
		coeffsString[i] = coeffsBigint[i].String()
	}

	return coeffsString
}

//Returns an string array containing the reconstructed coefficients
func (context *Context) PolyToBigint(p1 *Poly, coeffsBigint []*Int) {

	tmp := NewInt(0)

	for x := uint64(0); x < context.N; x++ {

		tmp.SetUint(0)
		coeffsBigint[x] = NewUint(0)

		for i := 0; i < len(context.Modulus); i++ {
			coeffsBigint[x].Add(coeffsBigint[x], tmp.Mul(NewUint(p1.Coeffs[i][x]), context.CrtReconstruction[i]))
		}

		coeffsBigint[x].Mod(coeffsBigint[x], context.ModulusBigint)
	}
}

// Returns an array containing the coefficients of Pol
// centered arount each [-Qi/2, Qi/2]
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

// Checks if two polynomials are equal in the given context
// Also applies an exact modular reduction on the two polynomials.
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
