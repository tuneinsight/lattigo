//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.package ckks
package ckks

import (
	"fmt"
	"github.com/hhcho/mpc-core"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"math/big"
)

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen int = 5

var pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"

// Encoder is an interface implenting the encoding algorithms.
type Encoder interface {
	Encode(plaintext *Plaintext, values []complex128, logSlots int)
	EncodeNew(values []complex128, logSlots int) (plaintext *Plaintext)
	EncodeAtLvlNew(level int, values []complex128, logSlots int) (plaintext *Plaintext)

	EncodeNTT(plaintext *Plaintext, values []complex128, logSlots int)
	EncodeNTTNew(values []complex128, logSlots int) (plaintext *Plaintext)
	EncodeNTTAtLvlNew(level int, values []complex128, logSlots int) (plaintext *Plaintext)

	EncodeDiagMatrixAtLvl(level int, vector map[int][]complex128, scale, maxM1N2Ratio float64, logSlots int) (matrix *PtDiagMatrix)

	Decode(plaintext *Plaintext, logSlots int) (res []complex128)
	DecodePublic(plaintext *Plaintext, logSlots int, sigma float64) []complex128

	Embed(values []complex128, logSlots int)
	ScaleUp(pol *ring.Poly, scale float64, moduli []uint64)

	WipeInternalMemory()

	EncodeCoeffs(values []float64, plaintext *Plaintext)
	DecodeCoeffs(plaintext *Plaintext) (res []float64)
	DecodeCoeffsPublic(plaintext *Plaintext, bound float64) (res []float64)

	GetErrSTDTimeDom(valuesWant, valuesHave []complex128, scale float64) (std float64)
	GetErrSTDFreqDom(valuesWant, valuesHave []complex128, scale float64) (std float64)

	EncodeRVecNew(values mpc_core.RVec, slots uint64, fracBits int) (plaintext *Plaintext)
	DecodeRVec(plaintext *Plaintext, slots uint64, fracBits int) (res mpc_core.RVec)
}

// EncoderBigComplex is an interface implenting the encoding algorithms with arbitrary precision.
type EncoderBigComplex interface {
	Encode(plaintext *Plaintext, values []*ring.Complex, logSlots int)
	EncodeNew(values []*ring.Complex, logSlots int) (plaintext *Plaintext)
	EncodeAtLvlNew(level int, values []*ring.Complex, logSlots int) (plaintext *Plaintext)
	EncodeNTT(plaintext *Plaintext, values []*ring.Complex, logSlots int)
	EncodeNTTAtLvlNew(level int, values []*ring.Complex, logSlots int) (plaintext *Plaintext)
	Decode(plaintext *Plaintext, logSlots int) (res []*ring.Complex)
	FFT(values []*ring.Complex, N int)
	InvFFT(values []*ring.Complex, N int)

	//EncodeCoeffs(values []*big.Float, plaintext *Plaintext)
	//DecodeCoeffs(plaintext *Plaintext) (res []*big.Float)
}

// encoder is a struct storing the necessary parameters to encode a slice of complex number on a Plaintext.
type encoder struct {
	params       *Parameters
	ringQ        *ring.Ring
	ringP        *ring.Ring
	bigintChain  []*big.Int
	bigintCoeffs []*big.Int
	qHalf        *big.Int
	polypool     *ring.Poly
	m            int
	rotGroup     []int

	gaussianSampler *ring.GaussianSampler
}

type encoderComplex128 struct {
	encoder
	values      []complex128
	valuesfloat []float64
	roots       []complex128
	rootsBig 	[]bigComplex
}

func newEncoder(params *Parameters) encoder {

	m := 2 * params.N()

	var q *ring.Ring
	var err error
	if q, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	var p *ring.Ring
	if params.PiCount() != 0 {
		if p, err = ring.NewRing(params.N(), params.pi); err != nil {
			panic(err)
		}
	}

	rotGroup := make([]int, m>>1)
	fivePows := 1
	for i := 0; i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= GaloisGen
		fivePows &= (m - 1)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	gaussianSampler := ring.NewGaussianSampler(prng)

	return encoder{
		params:          params.Copy(),
		ringQ:           q,
		ringP:           p,
		bigintChain:     genBigIntChain(params.qi),
		bigintCoeffs:    make([]*big.Int, m>>1),
		qHalf:           ring.NewUint(0),
		polypool:        q.NewPoly(),
		m:               m,
		rotGroup:        rotGroup,
		gaussianSampler: gaussianSampler,
	}
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a Plaintext.
func NewEncoder(params *Parameters) Encoder {

	encoder := newEncoder(params)

	var angle float64
	roots := make([]complex128, encoder.m+1)
	rootsBig := make([]bigComplex, encoder.m+1)
	anglePart := new(big.Float).Mul(piBig(128), big.NewFloat(2))
	anglePart.Quo(anglePart, big.NewFloat(float64(encoder.m)))

	for i := 0; i < encoder.m; i++ {
		angle = 2 * 3.141592653589793 * float64(i) / float64(encoder.m)

		roots[i] = complex(math.Cos(angle), math.Sin(angle))

		angleBig := new(big.Float)
		angleBig.Mul(anglePart, big.NewFloat(float64(i)))
		rootsBig[i] = bigComplex{cosBig(angleBig), sinBig(angleBig)}
	}
	roots[encoder.m] = roots[0]
	rootsBig[encoder.m] = rootsBig[0]

	return &encoderComplex128{
		encoder:     encoder,
		roots:       roots,
		rootsBig: 	 rootsBig,
		values:      make([]complex128, encoder.m>>2),
		valuesfloat: make([]float64, encoder.m>>1),
	}
}

// EncodeNew encodes a slice of complex128 of length slots = 2^{logSlots} on new plaintext at the maximum level.
func (encoder *encoderComplex128) EncodeNew(values []complex128, logSlots int) (plaintext *Plaintext) {
	return encoder.EncodeAtLvlNew(encoder.params.MaxLevel(), values, logSlots)
}

// EncodeAtLvlNew encodes a slice of complex128 of length slots = 2^{logSlots} on new plaintext at the desired level.
func (encoder *encoderComplex128) EncodeAtLvlNew(level int, values []complex128, logSlots int) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, level, encoder.params.scale)
	encoder.Encode(plaintext, values, logSlots)
	return
}

// Encode encodes a slice of complex128 of length slots = 2^{logSlots} on the input plaintext.
func (encoder *encoderComplex128) Encode(plaintext *Plaintext, values []complex128, logSlots int) {
	encoder.Embed(values, logSlots)
	encoder.ScaleUp(plaintext.value, plaintext.scale, encoder.ringQ.Modulus[:plaintext.Level()+1])
	encoder.WipeInternalMemory()
	plaintext.isNTT = false
}

// EncodeNTTNew encodes a slice of complex128 of length slots = 2^{logSlots} on new plaintext at the maximum level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderComplex128) EncodeNTTNew(values []complex128, logSlots int) (plaintext *Plaintext) {
	return encoder.EncodeNTTAtLvlNew(encoder.params.MaxLevel(), values, logSlots)
}

// EncodeNTTAtLvlNew encodes a slice of complex128 of length slots = 2^{logSlots} on new plaintext at the desired level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderComplex128) EncodeNTTAtLvlNew(level int, values []complex128, logSlots int) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.scale)
	encoder.EncodeNTT(plaintext, values, logSlots)
	return
}

// EncodeNTT encodes a slice of complex128 of length slots = 2^{logSlots} on the input plaintext.
// Returns a plaintext in the NTT domain.
func (encoder *encoderComplex128) EncodeNTT(plaintext *Plaintext, values []complex128, logSlots int) {
	encoder.Encode(plaintext, values, logSlots)
	encoder.ringQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// Embed encodes a vector and stores internally the encoded values.
// To be used in conjunction with ScaleUp.
func (encoder *encoderComplex128) Embed(values []complex128, logSlots int) {

	slots := 1 << logSlots

	if len(values) > encoder.params.N()/2 || len(values) > slots || logSlots > encoder.params.LogN()-1 {
		panic("cannot Encode: too many values/slots for the given ring degree")
	}

	for i := range values {
		encoder.values[i] = values[i]
	}

	invfft(encoder.values, slots, encoder.m, encoder.rotGroup, encoder.roots)

	gap := (encoder.ringQ.N >> 1) / slots

	for i, jdx, idx := 0, encoder.ringQ.N>>1, 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx] = real(encoder.values[i])
		encoder.valuesfloat[jdx] = imag(encoder.values[i])
	}
}

// GetErrSTDFreqDom returns the scaled standard deviation of the difference between two complex vectors in the slot domains
func (encoder *encoderComplex128) GetErrSTDFreqDom(valuesWant, valuesHave []complex128, scale float64) (std float64) {

	var err complex128
	for i := range valuesWant {
		err = valuesWant[i] - valuesHave[i]
		encoder.valuesfloat[2*i] = real(err)
		encoder.valuesfloat[2*i+1] = imag(err)
	}

	return StandardDeviation(encoder.valuesfloat[:len(valuesWant)*2], scale)
}

// GetErrSTDTimeDom returns the scaled standard deviation of the coefficient domain of the difference between two complex vectors in the slot domains
func (encoder *encoderComplex128) GetErrSTDTimeDom(valuesWant, valuesHave []complex128, scale float64) (std float64) {

	for i := range valuesHave {
		encoder.values[i] = (valuesWant[i] - valuesHave[i])
	}

	invfft(encoder.values, len(valuesWant), encoder.m, encoder.rotGroup, encoder.roots)

	for i := range valuesWant {
		encoder.valuesfloat[2*i] = real(encoder.values[i])
		encoder.valuesfloat[2*i+1] = imag(encoder.values[i])
	}

	return StandardDeviation(encoder.valuesfloat[:len(valuesWant)*2], scale)

}

// ScaleUp writes the internaly stored encoded values on a polynomial.
func (encoder *encoderComplex128) ScaleUp(pol *ring.Poly, scale float64, moduli []uint64) {
	scaleUpVecExact(encoder.valuesfloat, scale, moduli, pol.Coeffs)
}

// WipeInternalMemory sets the internally stored encoded values of the encoder to zero.
func (encoder *encoderComplex128) WipeInternalMemory() {
	for i := range encoder.values {
		encoder.values[i] = 0
	}

	for i := range encoder.valuesfloat {
		encoder.valuesfloat[i] = 0
	}
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
// Rounds the decimal part of the output (the bits under the scale) to "logPrecision" bits of precision.
func (encoder *encoderComplex128) DecodePublic(plaintext *Plaintext, logSlots int, bound float64) (res []complex128) {
	return encoder.decodePublic(plaintext, logSlots, bound)
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoderComplex128) Decode(plaintext *Plaintext, logSlots int) (res []complex128) {
	return encoder.decodePublic(plaintext, logSlots, 0)
}

func polyToComplexNoCRT(coeffs []uint64, values []complex128, scale float64, logSlots int, Q uint64) {

	slots := 1 << logSlots
	maxSlots := len(coeffs) >> 1
	gap := maxSlots / slots

	var real, imag float64
	for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {

		if coeffs[idx] >= Q>>1 {
			real = -float64(Q - coeffs[idx])
		} else {
			real = float64(coeffs[idx])
		}

		if coeffs[idx+maxSlots] >= Q>>1 {
			imag = -float64(Q - coeffs[idx+maxSlots])
		} else {
			imag = float64(coeffs[idx+maxSlots])
		}

		values[i] = complex(real, imag) / complex(scale, 0)
	}
}

func polyToComplexCRT(poly *ring.Poly, bigintCoeffs []*big.Int, values []complex128, scale float64, logSlots int, ringQ *ring.Ring, Q *big.Int) {

	ringQ.PolyToBigint(poly, bigintCoeffs)

	maxSlots := ringQ.N >> 1
	slots := 1 << logSlots
	gap := maxSlots / slots

	qHalf := new(big.Int)
	qHalf.Set(Q)
	qHalf.Rsh(qHalf, 1)

	var sign int

	for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {

		// Centers the value around the current modulus
		bigintCoeffs[idx].Mod(bigintCoeffs[idx], Q)
		sign = bigintCoeffs[idx].Cmp(qHalf)
		if sign == 1 || sign == 0 {
			bigintCoeffs[idx].Sub(bigintCoeffs[idx], Q)
		}

		// Centers the value around the current modulus
		bigintCoeffs[idx+maxSlots].Mod(bigintCoeffs[idx+maxSlots], Q)
		sign = bigintCoeffs[idx+maxSlots].Cmp(qHalf)
		if sign == 1 || sign == 0 {
			bigintCoeffs[idx+maxSlots].Sub(bigintCoeffs[idx+maxSlots], Q)
		}

		values[i] = complex(scaleDown(bigintCoeffs[idx], scale), scaleDown(bigintCoeffs[idx+maxSlots], scale))
	}
}

func (encoder *encoderComplex128) plaintextToComplex(level int, scale float64, logSlots int, p *ring.Poly, values []complex128) {
	if level == 0 {
		polyToComplexNoCRT(p.Coeffs[0], encoder.values, scale, logSlots, encoder.ringQ.Modulus[0])
	} else {
		polyToComplexCRT(p, encoder.bigintCoeffs, values, scale, logSlots, encoder.ringQ, encoder.bigintChain[level])
	}
}

func roundComplexVector(values []complex128, bound float64) {
	for i := range values {
		a := math.Round(real(values[i])*bound) / bound
		b := math.Round(imag(values[i])*bound) / bound
		values[i] = complex(a, b)
	}
}

func polyToFloatNoCRT(coeffs []uint64, values []float64, scale float64, Q uint64) {

	for i, c := range coeffs {

		if c >= Q>>1 {
			values[i] = -float64(Q-c) / scale
		} else {
			values[i] = float64(c) / scale
		}
	}
}

// PtDiagMatrix is a struct storing a plaintext diagonalized matrix
// ready to be evaluated on a ciphertext using evaluator.MultiplyByDiagMatrice.
type PtDiagMatrix struct {
	LogSlots   int                   // Log of the number of slots of the plaintext (needed to compute the appropriate rotation keys)
	N1         int                   // N1 is the number of inner loops of the baby-step giant-step algo used in the evaluation.
	Level      int                   // Level is the level at which the matrix is encoded (can be circuit dependant)
	Scale      float64               // Scale is the scale at which the matrix is encoded (can be circuit dependant)
	Vec        map[int][2]*ring.Poly // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non zero diagonal.
	naive      bool
	rotOnly    bool // Only one diagonal of the form [1, ..., 1]
	isGaussian bool // Each diagonal of the matrix is of the form [k, ..., k] for k a gaussian integer
}

// EncodeDiagMatrixAtLvl encodes a diagonalized plaintext matrix into PtDiagMatrix struct.
// It can then be evaluated on a ciphertext using evaluator.MultiplyByDiagMatrice.
// maxM1N2Ratio is the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.MultiplyByDiagMatrice.
// Optimal maxM1N2Ratio value is between 4 and 16 depending on the sparsity of the matrix.
func (encoder *encoderComplex128) EncodeDiagMatrixAtLvl(level int, vector map[int][]complex128, scale, maxM1N2Ratio float64, logSlots int) (matrix *PtDiagMatrix) {

	matrix = new(PtDiagMatrix)
	matrix.LogSlots = logSlots
	matrix.N1 = findbestbabygiantstepsplit(vector, 1<<logSlots, maxM1N2Ratio)

	var N, N1 int

	ringQ := encoder.ringQ

	// N1*N2 = N
	N = ringQ.N
	N1 = matrix.N1

	if len(vector) > 2 {

		index := make(map[int][]int)
		for key := range vector {
			idx1 := key / N1
			idx2 := key & (N1 - 1)
			if index[idx1] == nil {
				index[idx1] = []int{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
		}

		matrix.Vec = make(map[int][2]*ring.Poly)

		matrix.Level = level
		matrix.Scale = scale

		for j := range index {

			for _, i := range index[j] {
				matrix.Vec[N1*j+i] = encoder.encodeDiagonal(logSlots, level, scale, vector[N1*j+i], (N>>1)-(N1*j))
			}
		}
	} else {

		matrix.Vec = make(map[int][2]*ring.Poly)

		matrix.Level = level
		matrix.Scale = scale

		for i := range vector {
			matrix.Vec[i] = encoder.encodeDiagonal(logSlots, level, scale, vector[i], i)
		}
	}

	return
}

func (encoder *encoderComplex128) encodeDiagonal(logSlots, level int, scale float64, m []complex128, k int) [2]*ring.Poly {

	ringQ := encoder.ringQ
	ringP := encoder.ringP

	encoder.Embed(rotate(m, k), logSlots)

	mQ := ringQ.NewPolyLvl(level + 1)
	encoder.ScaleUp(mQ, scale, ringQ.Modulus[:level+1])
	ringQ.NTTLvl(level, mQ, mQ)
	ringQ.MFormLvl(level, mQ, mQ)

	mP := ringP.NewPoly()
	encoder.ScaleUp(mP, scale, ringP.Modulus)
	ringP.NTT(mP, mP)
	ringP.MForm(mP, mP)

	encoder.WipeInternalMemory()

	return [2]*ring.Poly{mQ, mP}
}

// Finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func findbestbabygiantstepsplit(vector map[int][]complex128, maxN int, maxRatio float64) (minN int) {

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		index := make(map[int][]int)

		for key := range vector {

			idx1 := key / N1
			idx2 := key & (N1 - 1)

			if index[idx1] == nil {
				index[idx1] = []int{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
		}

		if len(index[0]) > 0 {

			hoisted := len(index[0]) - 1
			normal := len(index) - 1

			// The matrice is very sparse already
			if normal == 0 {
				return N1 / 2
			}

			if hoisted > normal {
				// Finds the next split that has a ratio hoisted/normal greater or equal to maxRatio
				for float64(hoisted)/float64(normal) < maxRatio {

					if normal/2 == 0 {
						break
					}
					N1 *= 2
					hoisted = hoisted*2 + 1
					normal = normal / 2
				}
				return N1
			}
		}
	}

	return 1
}

func (encoder *encoderComplex128) decodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []complex128) {

	if logSlots > encoder.params.LogN()-1 {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	slots := 1 << logSlots

	if plaintext.isNTT {
		encoder.ringQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	} else {
		encoder.ringQ.CopyLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	}

	// B = floor(sigma * sqrt(2*pi))
	encoder.gaussianSampler.ReadAndAddLvl(plaintext.Level(), encoder.polypool, encoder.ringQ, sigma, int(2.5066282746310002*sigma))

	encoder.plaintextToComplex(plaintext.Level(), plaintext.Scale(), logSlots, encoder.polypool, encoder.values)

	fft(encoder.values, slots, encoder.m, encoder.rotGroup, encoder.roots)

	res = make([]complex128, slots)

	for i := range res {
		res[i] = encoder.values[i]
	}

	for i := range encoder.values {
		encoder.values[i] = 0
	}

	return
}

func invfft(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	var lenh, lenq, gap, idx int
	var u, v complex128

	for len := N; len >= 1; len >>= 1 {
		for i := 0; i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = M / lenq
			for j := 0; j < lenh; j++ {
				idx = (lenq - (rotGroup[j] % lenq)) * gap
				u = values[i+j] + values[i+j+lenh]
				v = values[i+j] - values[i+j+lenh]
				v *= roots[idx]
				values[i+j] = u
				values[i+j+lenh] = v

			}
		}
	}

	for i := 0; i < N; i++ {
		values[i] /= complex(float64(N), 0)
	}

	sliceBitReverseInPlaceComplex128(values, N)
}

func fft(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	var lenh, lenq, gap, idx int
	var u, v complex128

	sliceBitReverseInPlaceComplex128(values, N)

	for len := 2; len <= N; len <<= 1 {
		for i := 0; i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = M / lenq
			for j := 0; j < lenh; j++ {
				idx = (rotGroup[j] % lenq) * gap
				u = values[i+j]
				v = values[i+j+lenh]
				v *= roots[idx]
				values[i+j] = u + v
				values[i+j+lenh] = u - v
			}
		}
	}
}

// EncodeCoeffs takes as input a polynomial a0 + a1x + a2x^2 + ... + an-1x^n-1 with float coefficient
// and returns a scaled integer plaintext polynomial. Encodes at the input plaintext level.
func (encoder *encoderComplex128) EncodeCoeffs(values []float64, plaintext *Plaintext) {

	if len(values) > encoder.params.N() {
		panic("cannot EncodeCoeffs : too many values (maximum is N)")
	}

	scaleUpVecExact(values, plaintext.scale, encoder.ringQ.Modulus[:plaintext.Level()+1], plaintext.value.Coeffs)

	plaintext.isNTT = false
}

// EncodeCoeffsNTT takes as input a polynomial a0 + a1x + a2x^2 + ... + an-1x^n-1 with float coefficient
// and returns a scaled integer plaintext polynomial in NTT. Encodes at the input plaintext level.
func (encoder *encoderComplex128) EncodeCoeffsNTT(values []float64, plaintext *Plaintext) {
	encoder.EncodeCoeffs(values, plaintext)
	encoder.ringQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// DecodeCoeffsPublic takes as input a plaintext and returns the scaled down coefficient of the plaintext in float64.
// Rounds the decimal part of the output (the bits under the scale) to "logPrecision" bits of precision.
func (encoder *encoderComplex128) DecodeCoeffsPublic(plaintext *Plaintext, sigma float64) (res []float64) {
	return encoder.decodeCoeffsPublic(plaintext, sigma)
}

func (encoder *encoderComplex128) DecodeCoeffs(plaintext *Plaintext) (res []float64) {
	return encoder.decodeCoeffsPublic(plaintext, 0)
}

// DecodeCoeffs takes as input a plaintext and returns the scaled down coefficient of the plaintext in float64.
func (encoder *encoderComplex128) decodeCoeffsPublic(plaintext *Plaintext, sigma float64) (res []float64) {

	if plaintext.isNTT {
		encoder.ringQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	} else {
		encoder.ringQ.CopyLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	}

	if sigma != 0 {
		// B = floor(sigma * sqrt(2*pi))
		encoder.gaussianSampler.ReadAndAddLvl(plaintext.Level(), encoder.polypool, encoder.ringQ, sigma, int(2.5066282746310002*sigma))
	}

	res = make([]float64, encoder.params.N())

	// We have more than one moduli and need the CRT reconstruction
	if plaintext.Level() > 0 {

		encoder.ringQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

		Q := encoder.bigintChain[plaintext.Level()]

		encoder.qHalf.Set(Q)
		encoder.qHalf.Rsh(encoder.qHalf, 1)

		var sign int

		for i := range res {

			// Centers the value around the current modulus
			encoder.bigintCoeffs[i].Mod(encoder.bigintCoeffs[i], Q)

			sign = encoder.bigintCoeffs[i].Cmp(encoder.qHalf)
			if sign == 1 || sign == 0 {
				encoder.bigintCoeffs[i].Sub(encoder.bigintCoeffs[i], Q)
			}

			res[i] = scaleDown(encoder.bigintCoeffs[i], plaintext.scale)
		}
		// We can directly get the coefficients
	} else {

		Q := encoder.ringQ.Modulus[0]
		coeffs := encoder.polypool.Coeffs[0]

		for i := range res {

			if coeffs[i] >= Q>>1 {
				res[i] = -float64(Q - coeffs[i])
			} else {
				res[i] = float64(coeffs[i])
			}

			res[i] /= plaintext.scale
		}
	}

	return
}

func (encoder *encoderComplex128) EncodeRVecNew(values mpc_core.RVec, slots uint64, fracBits int) (plaintext *Plaintext) {
	if len(values) > encoder.params.Slots() || uint64(len(values)) > slots {
		panic("cannot EncodeRVecNew: too many values for the given number of slots")
	}

	if slots == 0 {
		return NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.Scale())
	}

	rtype := mpc_core.LElem128Zero

	if slots&(slots-1) != 0 { // slots not a power of two
		closestPow := uint64(math.Pow(2, math.Ceil(math.Log2(float64(slots)))))
		newValues := mpc_core.InitRVec(rtype.Zero(), int(closestPow))
		for i := range values {
			newValues[i] = values[i]
		}
		values = newValues
		slots = closestPow
	}

	if uint64(len(values)) != slots {
		panic("cannot EncodeRVecNew: number of values must be equal to slots")
	}

	if values.Type().TypeID() != rtype.TypeID() {
		panic("cannot EncodeRVecNew: only LElem128 supported")
	}

	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.Scale())

	slice := make([]bigComplex, len(encoder.values))
	zeroFloat := big.NewFloat(0)
	for i := range slice {
		if uint64(i) < slots {
			slice[i] = bigComplex{values[i].(mpc_core.LElem128).ToSignedBigFloat(fracBits), big.NewFloat(0)}
		} else {
			slice[i] = bigComplex{zeroFloat, zeroFloat}
		}
	}

	encoder.invfftBig(slice, slots)

	gap := uint64(encoder.params.Slots()) / slots

	floatSlice := make([]*big.Float, len(encoder.valuesfloat))
	for i := range floatSlice {
		floatSlice[i] = zeroFloat
	}

	for i, jdx, idx := uint64(0), uint64(encoder.params.Slots()), uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		floatSlice[idx] = slice[i].real
		floatSlice[jdx] = slice[i].imag
	}

	moduli := encoder.ringQ.Modulus[:plaintext.Level()+1]
	xInt := new(big.Int)
	xFlo := new(big.Float)
	scaleBig := big.NewFloat(plaintext.scale)
	tmp := new(big.Int)
	for i := range floatSlice {
		xFlo.Mul(floatSlice[i], scaleBig)
		xInt = arithRound(xFlo)

		for j := range moduli {
			tmp.Mod(xInt, ring.NewUint(moduli[j]))
			plaintext.value.Coeffs[j][i] = tmp.Uint64()
		}
	}

	encoder.ringQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	return
}


func (encoder *encoderComplex128) DecodeRVec(plaintext *Plaintext, slots uint64, fracBits int) (res mpc_core.RVec) {
	rtype := mpc_core.LElem128Zero

	encoder.ringQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	encoder.ringQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

	Q := encoder.bigintChain[plaintext.Level()]

	maxSlots := uint64(encoder.params.Slots())

	encoder.qHalf.Set(Q)
	encoder.qHalf.Rsh(encoder.qHalf, 1)

	gap := uint64(encoder.params.Slots()) / slots

	values := make([]bigComplex, len(encoder.values))

	var sign int

	for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx].Mod(encoder.bigintCoeffs[idx], Q)
		sign = encoder.bigintCoeffs[idx].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx].Sub(encoder.bigintCoeffs[idx], Q)
		}

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx+maxSlots].Mod(encoder.bigintCoeffs[idx+maxSlots], Q)
		sign = encoder.bigintCoeffs[idx+maxSlots].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx+maxSlots].Sub(encoder.bigintCoeffs[idx+maxSlots], Q)
		}

		values[i] = bigComplex{scaleDownBig(encoder.bigintCoeffs[idx], plaintext.scale), scaleDownBig(encoder.bigintCoeffs[idx+maxSlots], plaintext.scale)}
	}

	encoder.fftBig(values, slots)

	res = mpc_core.InitRVec(rtype.Zero(), int(slots))

	tmpFloat := new(big.Float)
	scale := big.NewFloat(0)
	scale.SetInt(new(big.Int).Lsh(big.NewInt(1), uint(fracBits)))
	for i := uint64(0); i < slots; i++ {
		v := values[i].real
		if v.Sign() < 0 {
			tmpFloat.Mul(tmpFloat.Neg(v), scale)
			res[i] = res[i].Sub(rtype.FromBigInt(arithRound(tmpFloat)))
		} else {
			tmpFloat.Mul(v, scale)
			res[i] = rtype.FromBigInt(arithRound(tmpFloat))
		}
	}

	return res
}

func (encoder *encoderComplex128) invfftlazyBig(values []bigComplex, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v bigComplex

	for len := N; len >= 1; len >>= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = uint64(encoder.m) / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (lenq - (uint64(encoder.rotGroup[j]) % lenq)) * gap
				//u = values[i+j] + values[i+j+lenh]
				u = values[i+j].copy()
				u.add(values[i+j+lenh])

				//v = values[i+j] - values[i+j+lenh]
				v = values[i+j].copy()
				v.sub(values[i+j+lenh])

				//v *= encoder.roots[idx]
				v.mul(encoder.rootsBig[idx])

				values[i+j] = u
				values[i+j+lenh] = v

			}
		}
	}

	sliceBitReverseInPlaceBigComplex(values, N)
}

func (encoder *encoderComplex128) invfftBig(values []bigComplex, N uint64) {
	encoder.invfftlazyBig(values, N)

	for i := uint64(0); i < N; i++ {
		//values[i] /= complex(float64(N), 0)
		values[i].real.Quo(values[i].real, big.NewFloat(float64(N)))
		values[i].imag.Quo(values[i].imag, big.NewFloat(float64(N)))
	}
}

func (encoder *encoderComplex128) fftBig(values []bigComplex, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v bigComplex

	sliceBitReverseInPlaceBigComplex(values, N)

	for len := uint64(2); len <= N; len <<= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = uint64(encoder.m) / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (uint64(encoder.rotGroup[j]) % lenq) * gap
				u = values[i+j]
				v = values[i+j+lenh].copy()

				//v *= encoder.roots[idx]
				v.mul(encoder.rootsBig[idx])

				//values[i+j] = u + v
				values[i+j] = u.copy()
				values[i+j].add(v)

				//values[i+j+lenh] = u - v
				values[i+j+lenh] = u.copy()
				values[i+j+lenh].sub(v)
			}
		}
	}
}

type encoderBigComplex struct {
	encoder
	zero            *big.Float
	cMul            *ring.ComplexMultiplier
	logPrecision    int
	values          []*ring.Complex
	valuesfloat     []*big.Float
	roots           []*ring.Complex
	gaussianSampler *ring.GaussianSampler
}

// NewEncoderBigComplex creates a new encoder using arbitrary precision complex arithmetic.
func NewEncoderBigComplex(params *Parameters, logPrecision int) EncoderBigComplex {
	encoder := newEncoder(params)

	var PI = new(big.Float)
	PI.SetPrec(uint(logPrecision))
	PI.SetString(pi)

	var PIHalf = new(big.Float)
	PIHalf.SetPrec(uint(logPrecision))
	PIHalf.SetString(pi)
	PIHalf.Quo(PIHalf, ring.NewFloat(2, logPrecision))

	var angle *big.Float
	roots := make([]*ring.Complex, encoder.m+1)
	for i := 0; i < encoder.m; i++ {
		angle = ring.NewFloat(2, logPrecision)
		angle.Mul(angle, PI)
		angle.Mul(angle, ring.NewFloat(float64(i), logPrecision))
		angle.Quo(angle, ring.NewFloat(float64(encoder.m), logPrecision))

		real := ring.Cos(angle)
		angle.Sub(PIHalf, angle)
		imag := ring.Cos(angle)

		roots[i] = ring.NewComplex(real, imag)
	}

	roots[encoder.m] = roots[0].Copy()

	values := make([]*ring.Complex, encoder.m>>2)
	valuesfloat := make([]*big.Float, encoder.m>>1)

	for i := 0; i < encoder.m>>2; i++ {

		values[i] = ring.NewComplex(ring.NewFloat(0, logPrecision), ring.NewFloat(0, logPrecision))
		valuesfloat[i*2] = ring.NewFloat(0, logPrecision)
		valuesfloat[(i*2)+1] = ring.NewFloat(0, logPrecision)
	}

	return &encoderBigComplex{
		encoder:      encoder,
		zero:         ring.NewFloat(0, logPrecision),
		cMul:         ring.NewComplexMultiplier(),
		logPrecision: logPrecision,
		roots:        roots,
		values:       values,
		valuesfloat:  valuesfloat,
	}
}

// EncodeNew encodes a slice of ring.Complex of length slots = 2^{logSlots} on a new plaintext at the maximum level.
func (encoder *encoderBigComplex) EncodeNew(values []*ring.Complex, logSlots int) (plaintext *Plaintext) {
	return encoder.EncodeAtLvlNew(encoder.params.MaxLevel(), values, logSlots)
}

// EncodeAtLvlNew encodes a slice of ring.Complex of length slots = 2^{logSlots} on a new plaintext at the desired level.
func (encoder *encoderBigComplex) EncodeAtLvlNew(level int, values []*ring.Complex, logSlots int) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, level, encoder.params.scale)
	encoder.Encode(plaintext, values, logSlots)
	return
}

// EncodeNTTNew encodes a slice of ring.Complex of length slots = 2^{logSlots} on a plaintext at the maximum level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderBigComplex) EncodeNTTNew(values []*ring.Complex, logSlots int) (plaintext *Plaintext) {
	return encoder.EncodeNTTAtLvlNew(encoder.params.MaxLevel(), values, logSlots)
}

// EncodeNTTAtLvlNew encodes a slice of ring.Complex of length slots = 2^{logSlots} on a plaintext at the desired level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderBigComplex) EncodeNTTAtLvlNew(level int, values []*ring.Complex, logSlots int) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.scale)
	encoder.EncodeNTT(plaintext, values, logSlots)
	return
}

// Encode encodes a slice of ring.Complex of length slots = 2^{logSlots} on a plaintext at the input plaintext level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderBigComplex) EncodeNTT(plaintext *Plaintext, values []*ring.Complex, logSlots int) {
	encoder.Encode(plaintext, values, logSlots)
	encoder.ringQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// Encode encodes a slice of ring.Complex of length slots = 2^{logSlots} on a plaintext at the input plaintext level.
func (encoder *encoderBigComplex) Encode(plaintext *Plaintext, values []*ring.Complex, logSlots int) {

	slots := 1 << logSlots

	if len(values) > encoder.params.N()/2 || len(values) > slots || logSlots > encoder.params.LogN()-1 {
		panic("cannot Encode: too many values/slots for the given ring degree")
	}

	if len(values) != slots {
		panic("cannot Encode: number of values must be equal to slots")
	}

	for i := 0; i < slots; i++ {
		encoder.values[i].Set(values[i])
	}

	encoder.InvFFT(encoder.values, slots)

	gap := (encoder.ringQ.N >> 1) / slots

	for i, jdx, idx := 0, (encoder.ringQ.N >> 1), 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx].Set(encoder.values[i].Real())
		encoder.valuesfloat[jdx].Set(encoder.values[i].Imag())
	}

	scaleUpVecExactBigFloat(encoder.valuesfloat, plaintext.scale, encoder.ringQ.Modulus[:plaintext.Level()+1], plaintext.value.Coeffs)

	coeffsBigInt := make([]*big.Int, encoder.params.N())

	encoder.ringQ.PolyToBigint(plaintext.value, coeffsBigInt)

	for i := 0; i < (encoder.ringQ.N >> 1); i++ {
		encoder.values[i].Real().Set(encoder.zero)
		encoder.values[i].Imag().Set(encoder.zero)
	}

	for i := 0; i < encoder.ringQ.N; i++ {
		encoder.valuesfloat[i].Set(encoder.zero)
	}
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
// Rounds the decimal part of the output (the bits under the scale) to "logPrecision" bits of precision.
func (encoder *encoderBigComplex) DecodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []*ring.Complex) {
	return encoder.decodePublic(plaintext, logSlots, sigma)
}

func (encoder *encoderBigComplex) Decode(plaintext *Plaintext, logSlots int) (res []*ring.Complex) {
	return encoder.decodePublic(plaintext, logSlots, 0)
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoderBigComplex) decodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []*ring.Complex) {

	slots := 1 << logSlots

	if logSlots > encoder.params.LogN()-1 {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	encoder.ringQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)

	if sigma != 0 {
		// B = floor(sigma * sqrt(2*pi))
		encoder.gaussianSampler.ReadAndAddLvl(plaintext.Level(), encoder.polypool, encoder.ringQ, sigma, int(2.5066282746310002*sigma+0.5))
	}

	encoder.ringQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

	Q := encoder.bigintChain[plaintext.Level()]

	maxSlots := encoder.ringQ.N >> 1

	scaleFlo := ring.NewFloat(plaintext.Scale(), encoder.logPrecision)

	encoder.qHalf.Set(Q)
	encoder.qHalf.Rsh(encoder.qHalf, 1)

	gap := maxSlots / slots

	var sign int

	for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx].Mod(encoder.bigintCoeffs[idx], Q)
		sign = encoder.bigintCoeffs[idx].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx].Sub(encoder.bigintCoeffs[idx], Q)
		}

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx+maxSlots].Mod(encoder.bigintCoeffs[idx+maxSlots], Q)
		sign = encoder.bigintCoeffs[idx+maxSlots].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx+maxSlots].Sub(encoder.bigintCoeffs[idx+maxSlots], Q)
		}

		encoder.values[i].Real().SetInt(encoder.bigintCoeffs[idx])
		encoder.values[i].Real().Quo(encoder.values[i].Real(), scaleFlo)

		encoder.values[i].Imag().SetInt(encoder.bigintCoeffs[idx+maxSlots])
		encoder.values[i].Imag().Quo(encoder.values[i].Imag(), scaleFlo)
	}

	encoder.FFT(encoder.values, slots)

	res = make([]*ring.Complex, slots)

	for i := range res {
		res[i] = encoder.values[i].Copy()
	}

	for i := 0; i < maxSlots; i++ {
		encoder.values[i].Real().Set(encoder.zero)
		encoder.values[i].Imag().Set(encoder.zero)
	}

	return
}

// InvFFT evaluates the encoding matrix on a slice fo ring.Complex values.
func (encoder *encoderBigComplex) InvFFT(values []*ring.Complex, N int) {

	var lenh, lenq, gap, idx int
	u := ring.NewComplex(nil, nil)
	v := ring.NewComplex(nil, nil)

	for len := N; len >= 1; len >>= 1 {
		for i := 0; i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := 0; j < lenh; j++ {
				idx = (lenq - (encoder.rotGroup[j] % lenq)) * gap
				u.Add(values[i+j], values[i+j+lenh])
				v.Sub(values[i+j], values[i+j+lenh])
				encoder.cMul.Mul(v, encoder.roots[idx], v)
				values[i+j].Set(u)
				values[i+j+lenh].Set(v)
			}
		}
	}

	NBig := ring.NewFloat(float64(N), encoder.logPrecision)
	for i := range values {
		values[i][0].Quo(values[i][0], NBig)
		values[i][1].Quo(values[i][1], NBig)
	}

	sliceBitReverseInPlaceRingComplex(values, N)
}

// FFT evaluates the decoding matrix on a slice fo ring.Complex values.
func (encoder *encoderBigComplex) FFT(values []*ring.Complex, N int) {

	var lenh, lenq, gap, idx int

	u := ring.NewComplex(nil, nil)
	v := ring.NewComplex(nil, nil)

	sliceBitReverseInPlaceRingComplex(values, N)

	for len := 2; len <= N; len <<= 1 {
		for i := 0; i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := 0; j < lenh; j++ {
				idx = (encoder.rotGroup[j] % lenq) * gap
				u.Set(values[i+j])
				v.Set(values[i+j+lenh])
				encoder.cMul.Mul(v, encoder.roots[idx], v)
				values[i+j].Add(u, v)
				values[i+j+lenh].Sub(u, v)
			}
		}
	}
}

type bigComplex struct {
	real *big.Float
	imag *big.Float
}

var piBase, _ = new(big.Float).SetPrec(200).SetString("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194912")
var piTable = make(map[uint]*big.Float)

func piBig(prec uint) *big.Float {
	if prec > piBase.Prec() {
		//fmt.Println("Warning: requested precision for Pi ", prec, " exceeds constant ", piBase.Prec())
		prec = piBase.Prec()
	}

	v, flag := piTable[prec]
	if flag {
		return v
	}

	fmt.Println("Pi", prec)

	v = new(big.Float)
	v.Copy(piBase)
	v.SetPrec(prec)
	piTable[prec] = v

	return v
}
func sinBig(x *big.Float) *big.Float {
	tmp := new(big.Float)
	tmp.Mul(piBig(x.Prec()), big.NewFloat(0.5))
	tmp.Sub(x, tmp)
	return cosBig(tmp) // six(x) = cos(x - pi/2)
}
func cosBig(x *big.Float) *big.Float {
	x = new(big.Float).Abs(x)
	tmp := new(big.Float)
	pi := piBig(x.Prec())
	flag2pi := x.Cmp(tmp.Mul(pi, big.NewFloat(2)))
	flagpi := x.Cmp(pi)
	flagpi2 := x.Cmp(tmp.Mul(pi, big.NewFloat(0.5)))
	if flag2pi > 0 {
		panic("cosBig input outside of [-2*pi, 2*pi]")
	} else if flag2pi == 0 {
		return big.NewFloat(1)
	} else if flagpi >= 0 {
		return tmp.Neg(cosBig(tmp.Sub(x, pi)))
	} else if flagpi2 == 0 {
		return big.NewFloat(0)
	} else if flagpi2 > 0 {
		return tmp.Neg(cosBig(tmp.Sub(x, pi)))
	}

	// x is now in [0, pi/2]

	// Choose the number of steps k of the Cordic algorithm
	// k=250 gives about 150 decimals, k=1000 about 600
	steps := 250 // TODO: choose based on the precision of x

	// xˆ2/2ˆ(2k)
	t := new(big.Float)
	t.SetMantExp(big.NewFloat(1), -2*steps)
	s := new(big.Float)
	tmp.Mul(x, x)
	s.Mul(tmp, t)

	four := big.NewFloat(4)
	for i := 1; i <= steps; i++ {
		tmp.Sub(four, s)
		s.Mul(s, tmp)
	}

	out := new(big.Float)
	tmp.Quo(s, big.NewFloat(2))
	out.Sub(big.NewFloat(1), tmp)
	return out
}

/*func ToSignedBigInt(a mpc_core.LElem128) *big.Int {
	tmp := a.ToBigInt()
	half := new(big.Int).Rsh(mpc_core.LElem128ModBig, 1)
	if tmp.Cmp(half) >= 0 {
		tmp.Sub(tmp, mpc_core.LElem128ModBig)
	}
	return tmp
}*/
/*func ToSignedBigFloat(a mpc_core.LElem128, fracBits int) *big.Float {
	out := new(big.Float).SetInt(ToSignedBigInt(a))
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}*/

func (bc bigComplex) copy() (out bigComplex) {
	out.imag = new(big.Float)
	out.real = new(big.Float)
	out.imag.Copy(bc.imag)
	out.real.Copy(bc.real)
	return out
}

func (bc *bigComplex) fromComplex128(v complex128) {
	bc.imag = big.NewFloat(imag(v))
	bc.real = big.NewFloat(real(v))
}

func (bc *bigComplex) add(v bigComplex) {
	bc.real.Add(bc.real, v.real)
	bc.imag.Add(bc.imag, v.imag)
}

func (bc *bigComplex) sub(v bigComplex) {
	bc.real.Sub(bc.real, v.real)
	bc.imag.Sub(bc.imag, v.imag)
}

func (bc *bigComplex) mul(v bigComplex) {
	t1 := new(big.Float)
	t2 := new(big.Float)
	outr := new(big.Float)
	outi := new(big.Float)

	t1.Mul(bc.real, v.real)
	t2.Mul(bc.imag, v.imag)
	outr.Sub(t1, t2)

	t1.Mul(bc.real, v.imag)
	t2.Mul(bc.imag, v.real)
	outi.Add(t1, t2)

	bc.real = outr
	bc.imag = outi
}

func sliceBitReverseInPlaceBigComplex(slice []bigComplex, N uint64) {

	var bit, j uint64

	for i := uint64(1); i < N; i++ {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}
	}
}

func scaleDownBig(coeff *big.Int, n float64) (x *big.Float) {

	x = new(big.Float).SetInt(coeff)
	x.Quo(x, big.NewFloat(n))

	return
}

var halfFloat = big.NewFloat(0.5)
func arithRound(a *big.Float) *big.Int {
	var i *big.Int
	if a.Signbit() {
		i, _ = new(big.Float).Sub(a, halfFloat).Int(nil)
	} else {
		i, _ = new(big.Float).Add(a, halfFloat).Int(nil)
	}
	return i
}
