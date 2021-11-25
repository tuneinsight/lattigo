//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.package ckks
package ckks

import (
	"math"
	"math/big"
	"math/bits"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

var pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"

// Encoder is an interface implementing the encoding and decoding operations. It provides method to embed []complex128 and []float64 types
// into plaintexts and the inverse operations.
// Two different encodings are provided:
//
//     - Slots: The coefficients are first subjected to a special Fourier transform before being embedded on the plaintext. This encoding can embed
//              []complex128 and []float64 slices of size at most N/2 (N being the ring degree) and leverages the convolution property of the DFT
//              to preserve point-wise complex multiplication in the plaintext domain, i.e. a ciphertext multiplication will result in an element-
//              wise multiplication in the plaintext domain. It also enables plaintext slots cyclic rotations. Other operations, like addition or
//              constant multiplication behave as usual.
//
//     - Coeffs: The coefficients are directly embedded on the plaintext. This encoding only allows to encode []float64 slices, but of size up to N
//               (N being the ring degree) and does not preserve the point-wise multiplication. A ciphertext multiplication will result in a nega-
//               cyclic polynomial convolution in the plaintext domain. This encoding does not provide native slot cyclic rotation.
//               Other operations, like addition or constant multiplication behave as usual.
//
// The figure bellow illustrates the relationship between those two encoding:
//
//                                              Real^{N}          Z_Q[X]/(X^N+1)
// EncodeCoeffs: ----------------------------->[]float64 ---------> Plaintext
//                                                 |
//                    Complex^{N/2}                |
// EncodeSlots:  []complex128/[]float64 -> iDFT ---┘
type Encoder interface {

	// Slots Encoding

	// Encode encodes a set of values on the target plaintext.
	// This method is identical to "EncodeSlots".
	// Encoding is done at the level and scale of the plaintext.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
	// values.(type) can be either []complex128 of []float64.
	// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
	// Returned plaintext is always in the NTT domain.
	Encode(plaintext *Plaintext, values interface{}, logSlots int)

	// EncodeSlots encodes a set of values on a new plaintext.
	// This method is identical to "EncodeSlotsNew".
	// Encoding is done at the provided level and with the provided scale.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
	// values.(type) can be either []complex128 of []float64.
	// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
	// Returned plaintext is always in the NTT domain.
	EncodeNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext)

	// Encode encodes a set of values on the target plaintext.
	// Encoding is done at the level and scale of the plaintext.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
	// values.(type) can be either []complex128 of []float64.
	// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
	// Returned plaintext is always in the NTT domain.
	EncodeSlots(plaintext *Plaintext, values interface{}, logSlots int)

	// EncodeSlotsNew encodes a set of values on a new plaintext.
	// Encoding is done at the provided level and with the provided scale.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
	// values.(type) can be either []complex128 of []float64.
	// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
	// Returned plaintext is always in the NTT domain.
	EncodeSlotsNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext)

	// EncodeDiagMatrixBSGS encodes a diagonalized plaintext matrix into PtDiagMatrix struct.
	// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
	// Evaluation will use the optimized approach (double hoisting and baby-step giant-step).
	// Faster if there is more than a few non-zero diagonals.
	// maxM1N2Ratio is the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.LinearTransform.
	// Optimal maxM1N2Ratio value is between 4 and 16 depending on the sparsity of the matrix.
	EncodeDiagMatrixBSGS(level int, diagMatrix map[int][]complex128, scale, maxM1N2Ratio float64, logSlots int) (matrix PtDiagMatrix)

	// EncodeDiagMatrix encodes a diagonalized plaintext matrix into PtDiagMatrix struct.
	// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
	// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
	// Faster if there is only a few non-zero diagonals but uses more keys.
	EncodeDiagMatrix(level int, vector map[int][]complex128, scale float64, logSlots int) (matrix PtDiagMatrix)

	// Decode decodes the input plaintext on a new slice of complex128.
	Decode(plaintext *Plaintext, logSlots int) (res []complex128)

	// DecodePublic decodes the input plaintext on a new slice of complex128.
	// Adds, before the decoding step, an error with standard deviation sigma.
	// If the underlying ringType is ConjugateInvariant, the imaginary part (and
	// its related error) are zero.
	DecodePublic(plaintext *Plaintext, logSlots int, sigma float64) []complex128

	// Coeffs Encoding

	// EncodeCoeffs encodes the values on the coefficient of the plaintext.
	// Encoding is done at the level and scale of the plaintext.
	// User must ensure that 1<= len(values) <= 2^LogN
	EncodeCoeffs(values []float64, plaintext *Plaintext)

	// EncodeCoeffsNew encodes the values on the coefficient of a new plaintext.
	// Encoding is done at the provided level and with the provided scale.
	// User must ensure that 1<= len(values) <= 2^LogN
	EncodeCoeffsNew(values []float64, level int, scale float64) (plaintext *Plaintext)

	// DecodeCoeffs reconstructs the RNS coefficients of the plaintext on a slice of float64.
	DecodeCoeffs(plaintext *Plaintext) (res []float64)

	// DecodeCoeffsPublic reconstructs the RNS coefficients of the plaintext on a slice of float64.
	// Adds an error with standard deviation sigma.
	DecodeCoeffsPublic(plaintext *Plaintext, bound float64) (res []float64)

	// GetErrSTDCoeffDomain returns StandardDeviation(Encode(valuesWant-valuesHave))*scale
	// Which is the scaled standard deviation in the coefficient domain of the difference
	// of two complex vector in the slot domain.
	GetErrSTDCoeffDomain(valuesWant, valuesHave []complex128, scale float64) (std float64)

	// GetErrSTDSlotDomain returns StandardDeviation(valuesWant-valuesHave)*scale
	// Which is the scaled standard deviation of two complex vectors.
	GetErrSTDSlotDomain(valuesWant, valuesHave []complex128, scale float64) (std float64)
}

// EncoderBigComplex is an interface implementing the encoding algorithms with arbitrary precision.
type EncoderBigComplex interface {

	// Encode encodes a set of values on the target plaintext.
	// Encoding is done at the level and scale of the plaintext.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^LogN.
	Encode(plaintext *Plaintext, values []*ring.Complex, logSlots int)

	// EncodeNew encodes a set of values on a new plaintext.
	// Encoding is done at the provided level and with the provided scale.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^LogN.
	EncodeNew(values []*ring.Complex, level int, scale float64, logSlots int) (plaintext *Plaintext)

	// Decode decodes the input plaintext on a new slice of ring.Complex.
	Decode(plaintext *Plaintext, logSlots int) (res []*ring.Complex)

	// FFT evaluates the decoding matrix on a slice fo ring.Complex values.
	FFT(values []*ring.Complex, N int)

	// InvFFT evaluates the encoding matrix on a slice fo ring.Complex values.
	InvFFT(values []*ring.Complex, N int)
}

// encoder is a struct storing the necessary parameters to encode a slice of complex number on a Plaintext.
type encoder struct {
	params       Parameters
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
	valuesFloat []float64
	roots       []complex128
}

func newEncoder(params Parameters) encoder {

	m := int(params.RingQ().NthRoot)

	rotGroup := make([]int, m>>1)
	fivePows := 1
	for i := 0; i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= int(GaloisGen)
		fivePows &= (m - 1)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	gaussianSampler := ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))

	return encoder{
		params:          params,
		bigintChain:     genBigIntChain(params.Q()),
		bigintCoeffs:    make([]*big.Int, m>>1),
		qHalf:           ring.NewUint(0),
		polypool:        params.RingQ().NewPoly(),
		m:               m,
		rotGroup:        rotGroup,
		gaussianSampler: gaussianSampler,
	}
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a Plaintext.
func NewEncoder(params Parameters) Encoder {

	encoder := newEncoder(params)

	var angle float64
	roots := make([]complex128, encoder.m+1)
	for i := 0; i < encoder.m; i++ {
		angle = 2 * 3.141592653589793 * float64(i) / float64(encoder.m)

		roots[i] = complex(math.Cos(angle), math.Sin(angle))
	}
	roots[encoder.m] = roots[0]

	return &encoderComplex128{
		encoder:     encoder,
		roots:       roots,
		values:      make([]complex128, encoder.m>>2),
		valuesFloat: make([]float64, encoder.m>>1),
	}
}

func (encoder *encoderComplex128) EncodeNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, level, scale)
	encoder.Encode(plaintext, values, logSlots)
	return
}

func (encoder *encoderComplex128) Encode(plaintext *Plaintext, values interface{}, logSlots int) {
	encoder.embed(values, logSlots)
	scaleUpVecExact(encoder.valuesFloat[:encoder.params.N()], plaintext.Scale, encoder.params.RingQ().Modulus[:plaintext.Level()+1], plaintext.Value.Coeffs)
	encoder.params.RingQ().NTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
	plaintext.Value.IsNTT = true
}

func (encoder *encoderComplex128) EncodeSlotsNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext) {
	return encoder.EncodeSlotsNew(values, level, scale, logSlots)
}

func (encoder *encoderComplex128) EncodeSlots(plaintext *Plaintext, values interface{}, logSlots int) {
	encoder.Encode(plaintext, values, logSlots)
}

func (encoder *encoderComplex128) embed(values interface{}, logSlots int) {

	slots := 1 << logSlots

	N := encoder.params.RingQ().N

	var gap int
	switch encoder.params.RingType() {
	case ring.Standard:
		gap = (N >> 1) / slots
	case ring.ConjugateInvariant:
		gap = N / slots
	default:
		panic("invalid ring type")
	}

	// First checks the type of input values
	switch values := values.(type) {

	// If complex
	case []complex128:

		// Checks that the number of values is with the possible range
		if len(values) > int(encoder.params.RingQ().NthRoot>>1) || len(values) > slots || slots > int(encoder.params.RingQ().NthRoot>>2) {
			panic("cannot Encode: too many values/slots for the given ring degree")
		}

		switch encoder.params.RingType() {

		case ring.Standard:

			copy(encoder.values[:len(values)], values)

			for i := len(values); i < slots; i++ {
				encoder.values[i] = 0
			}

			invfft(encoder.values, slots, encoder.m, encoder.rotGroup, encoder.roots)

			for i, idx, jdx := 0, 0, N>>1; i < slots; i, idx, jdx = i+1, idx+gap, jdx+gap {
				encoder.valuesFloat[idx] = real(encoder.values[i])
				encoder.valuesFloat[jdx] = imag(encoder.values[i])
			}

		case ring.ConjugateInvariant:

			// Discards the imaginary part
			for i := range values {
				encoder.values[i] = complex(real(values[i]), 0)
			}

			for i := len(values); i < slots; i++ {
				encoder.values[i] = 0
			}

			invfft(encoder.values, slots, encoder.m, encoder.rotGroup, encoder.roots)

			for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
				encoder.valuesFloat[idx] = real(encoder.values[i])
			}

		// Else panics
		default:
			panic("unsuported ringType")
		}

	// If floats only
	case []float64:

		if len(values) > int(encoder.params.RingQ().NthRoot>>1) || len(values) > slots || slots > int(encoder.params.RingQ().NthRoot>>2) {
			panic("cannot Encode: too many values/slots for the given ring degree")
		}

		for i := range values {
			encoder.values[i] = complex(values[i], 0)
		}

		for i := len(values); i < slots; i++ {
			encoder.values[i] = 0
		}

		invfft(encoder.values, slots, encoder.m, encoder.rotGroup, encoder.roots)

		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			encoder.valuesFloat[idx] = real(encoder.values[i])
		}

		if encoder.params.RingType() == ring.Standard {
			for i, jdx := 0, N>>1; i < slots; i, jdx = i+1, jdx+gap {
				encoder.valuesFloat[jdx] = imag(encoder.values[i])
			}
		}

	default:
		panic("values must be []complex128 or []float64")
	}
}

func (encoder *encoderComplex128) GetErrSTDSlotDomain(valuesWant, valuesHave []complex128, scale float64) (std float64) {
	var err complex128
	for i := range valuesWant {
		err = valuesWant[i] - valuesHave[i]
		encoder.valuesFloat[2*i] = real(err)
		encoder.valuesFloat[2*i+1] = imag(err)
	}
	return StandardDeviation(encoder.valuesFloat[:len(valuesWant)*2], scale)
}

func (encoder *encoderComplex128) GetErrSTDCoeffDomain(valuesWant, valuesHave []complex128, scale float64) (std float64) {

	for i := range valuesHave {
		encoder.values[i] = (valuesWant[i] - valuesHave[i])
	}

	for i := len(valuesHave); i < len(encoder.values); i++ {
		encoder.values[i] = complex(0, 0)
	}

	// Runs FFT^-1 with the smallest power of two length that is greater than the input size
	invfft(encoder.values, 1<<bits.Len64(uint64(len(valuesHave)-1)), encoder.m, encoder.rotGroup, encoder.roots)

	for i := range valuesWant {
		encoder.valuesFloat[2*i] = real(encoder.values[i])
		encoder.valuesFloat[2*i+1] = imag(encoder.values[i])
	}

	return StandardDeviation(encoder.valuesFloat[:len(valuesWant)*2], scale)

}

// DecodePublic decodes the Plaintext values to a slice of complex128 values of size at most N/2.
// Adds a Gaussian error to the plaintext of variance sigma and bound floor(sqrt(2*pi)*sigma) before decoding
func (encoder *encoderComplex128) DecodePublic(plaintext *Plaintext, logSlots int, bound float64) (res []complex128) {
	return encoder.decodePublic(plaintext, logSlots, bound)
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoderComplex128) Decode(plaintext *Plaintext, logSlots int) (res []complex128) {
	return encoder.decodePublic(plaintext, logSlots, 0)
}

func polyToComplexNoCRT(coeffs []uint64, values []complex128, scale float64, logSlots int, isreal bool, ringQ *ring.Ring) {

	slots := 1 << logSlots
	maxSlots := int(ringQ.NthRoot >> 2)
	gap := maxSlots / slots
	Q := ringQ.Modulus[0]

	var real float64
	for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {

		if coeffs[idx] >= Q>>1 {
			real = -float64(Q - coeffs[idx])
		} else {
			real = float64(coeffs[idx])
		}

		values[i] = complex(real, 0)
	}

	if !isreal {
		var imag float64
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {

			if coeffs[idx+maxSlots] >= Q>>1 {
				imag = -float64(Q - coeffs[idx+maxSlots])
			} else {
				imag = float64(coeffs[idx+maxSlots])
			}

			values[i] += complex(0, imag)
		}
	}

	for i := 0; i < slots; i++ {
		values[i] /= complex(scale, 0)
	}
}

func polyToComplexCRT(poly *ring.Poly, bigintCoeffs []*big.Int, values []complex128, scale float64, logSlots int, isreal bool, ringQ *ring.Ring, Q *big.Int) {

	ringQ.PolyToBigint(poly, bigintCoeffs)

	maxSlots := int(ringQ.NthRoot >> 2)
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

		values[i] = complex(scaleDown(bigintCoeffs[idx], scale), 0)
	}

	if !isreal {
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			// Centers the value around the current modulus
			bigintCoeffs[idx+maxSlots].Mod(bigintCoeffs[idx+maxSlots], Q)
			sign = bigintCoeffs[idx+maxSlots].Cmp(qHalf)
			if sign == 1 || sign == 0 {
				bigintCoeffs[idx+maxSlots].Sub(bigintCoeffs[idx+maxSlots], Q)
			}
			values[i] += complex(0, scaleDown(bigintCoeffs[idx+maxSlots], scale))
		}
	}
}

func (encoder *encoderComplex128) plaintextToComplex(level int, scale float64, logSlots int, p *ring.Poly, values []complex128) {

	isreal := encoder.params.RingType() == ring.ConjugateInvariant
	if level == 0 {
		polyToComplexNoCRT(p.Coeffs[0], values, scale, logSlots, isreal, encoder.params.RingQ())
	} else {
		polyToComplexCRT(p, encoder.bigintCoeffs, values, scale, logSlots, isreal, encoder.params.RingQ(), encoder.bigintChain[level])
	}

	if isreal { // [X]/(X^N+1) to [X+X^-1]/(X^N+1)
		tmp := encoder.values
		slots := 1 << logSlots
		for i := 1; i < slots; i++ {
			tmp[i] -= complex(0, real(tmp[slots-i]))
		}
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
	LogSlots   int                 // Log of the number of slots of the plaintext (needed to compute the appropriate rotation keys)
	N1         int                 // N1 is the number of inner loops of the baby-step giant-step algo used in the evaluation.
	Level      int                 // Level is the level at which the matrix is encoded (can be circuit dependent)
	Scale      float64             // Scale is the scale at which the matrix is encoded (can be circuit dependent)
	Vec        map[int]rlwe.PolyQP // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non zero diagonal.
	Naive      bool
	isGaussian bool // Each diagonal of the matrix is of the form [k, ..., k] for k a Gaussian integer
}

// BsgsIndex returns the index map and needed rotation for the BSGS matrix-vector multiplication algorithm.
func BsgsIndex(el interface{}, slots, N1 int) (index map[int][]int, rotations []int) {
	index = make(map[int][]int)
	rotations = []int{}
	switch element := el.(type) {
	case map[int][]complex128:
		for key := range element {
			key &= (slots - 1)
			idx1 := key / N1
			idx2 := key & (N1 - 1)
			if index[idx1] == nil {
				index[idx1] = []int{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
			if !utils.IsInSliceInt(idx2, rotations) {
				rotations = append(rotations, idx2)
			}
		}
	case map[int]bool:
		for key := range element {
			key &= (slots - 1)
			idx1 := key / N1
			idx2 := key & (N1 - 1)
			if index[idx1] == nil {
				index[idx1] = []int{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
			if !utils.IsInSliceInt(idx2, rotations) {
				rotations = append(rotations, idx2)
			}
		}
	case map[int]rlwe.PolyQP:
		for key := range element {
			key &= (slots - 1)
			idx1 := key / N1
			idx2 := key & (N1 - 1)
			if index[idx1] == nil {
				index[idx1] = []int{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
			if !utils.IsInSliceInt(idx2, rotations) {
				rotations = append(rotations, idx2)
			}
		}
	case []int:
		for key := range element {
			key &= (slots - 1)
			idx1 := key / N1
			idx2 := key & (N1 - 1)
			if index[idx1] == nil {
				index[idx1] = []int{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
			if !utils.IsInSliceInt(idx2, rotations) {
				rotations = append(rotations, idx2)
			}
		}
	}

	return
}

func (encoder *encoderComplex128) EncodeDiagMatrixBSGS(level int, diagMatrix map[int][]complex128, scale, maxM1N2Ratio float64, logSlots int) (matrix PtDiagMatrix) {

	slots := 1 << logSlots

	// N1*N2 = N
	n1 := FindBestBSGSSplit(diagMatrix, slots, maxM1N2Ratio)

	index, _ := BsgsIndex(diagMatrix, slots, n1)

	vec := make(map[int]rlwe.PolyQP)

	for j := range index {

		for _, i := range index[j] {

			// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
			v := diagMatrix[n1*j+i]
			if len(v) == 0 {
				v = diagMatrix[(n1*j+i)-slots]
			}

			if len(v) != slots {
				panic("diagMatrix []complex slices mut have len '1<<logSlots'")
			}

			vec[n1*j+i] = encoder.encodeDiagonalSingle(logSlots, level, scale, utils.RotateComplex128Slice(v, -n1*j))
		}
	}

	return PtDiagMatrix{LogSlots: logSlots, N1: n1, Vec: vec, Level: level, Scale: scale}
}

func (encoder *encoderComplex128) EncodeDiagMatrix(level int, diagMatrix map[int][]complex128, scale float64, logSlots int) (matrix PtDiagMatrix) {

	vec := make(map[int]rlwe.PolyQP)
	slots := 1 << logSlots
	for i := range diagMatrix {

		idx := i
		if idx < 0 {
			idx += slots
		}
		vec[idx] = encoder.encodeDiagonalSingle(logSlots, level, scale, diagMatrix[i])
	}

	return PtDiagMatrix{LogSlots: logSlots, N1: 0, Vec: vec, Level: level, Scale: scale, Naive: true}
}

func (encoder *encoderComplex128) encodeDiagonalSingle(logSlots, level int, scale float64, m []complex128) (vecQP rlwe.PolyQP) {

	levelQ := level
	levelP := encoder.params.PCount() - 1
	ringQP := encoder.params.RingQP()

	encoder.embed(m, logSlots)

	vecQP = ringQP.NewPolyLvl(levelQ, levelP)
	scaleUpVecExact(encoder.valuesFloat[:encoder.params.N()], scale, encoder.params.RingQ().Modulus[:level+1], vecQP.Q.Coeffs)
	scaleUpVecExact(encoder.valuesFloat[:encoder.params.N()], scale, encoder.params.RingP().Modulus, vecQP.P.Coeffs)
	ringQP.NTTLvl(levelQ, levelP, vecQP, vecQP)
	ringQP.MFormLvl(levelQ, levelP, vecQP, vecQP)

	return
}

// FindBestBSGSSplit finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func FindBestBSGSSplit(diagMatrix interface{}, maxN int, maxRatio float64) (minN int) {

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		index, _ := BsgsIndex(diagMatrix, maxN, N1)

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

	slots := 1 << logSlots

	if slots > int(encoder.params.RingQ().NthRoot>>2) {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	if plaintext.Value.IsNTT {
		encoder.params.RingQ().InvNTTLvl(plaintext.Level(), plaintext.Value, encoder.polypool)
	} else {
		ring.CopyValuesLvl(plaintext.Level(), plaintext.Value, encoder.polypool)
	}

	// B = floor(sigma * sqrt(2*pi))
	if sigma != 0 {
		encoder.gaussianSampler.ReadAndAddFromDistLvl(plaintext.Level(), encoder.polypool, encoder.params.RingQ(), sigma, int(2.5066282746310002*sigma))
	}

	encoder.plaintextToComplex(plaintext.Level(), plaintext.Scale, logSlots, encoder.polypool, encoder.values)

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

	SliceBitReverseInPlaceComplex128(values, N)
}

func fft(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	var lenh, lenq, gap, idx int
	var u, v complex128

	SliceBitReverseInPlaceComplex128(values, N)

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
	scaleUpVecExact(values, plaintext.Scale, encoder.params.RingQ().Modulus[:plaintext.Level()+1], plaintext.Value.Coeffs)
	encoder.params.RingQ().NTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
	plaintext.Value.IsNTT = true
}

func (encoder *encoderComplex128) EncodeCoeffsNew(values []float64, level int, scale float64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, level, scale)
	encoder.EncodeCoeffs(values, plaintext)
	return
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

	if plaintext.Value.IsNTT {
		encoder.params.RingQ().InvNTTLvl(plaintext.Level(), plaintext.Value, encoder.polypool)
	} else {
		ring.CopyValuesLvl(plaintext.Level(), plaintext.Value, encoder.polypool)
	}

	if sigma != 0 {
		// B = floor(sigma * sqrt(2*pi))
		encoder.gaussianSampler.ReadAndAddFromDistLvl(plaintext.Level(), encoder.polypool, encoder.params.RingQ(), sigma, int(2.5066282746310002*sigma))
	}

	res = make([]float64, encoder.params.N())

	// We have more than one moduli and need the CRT reconstruction
	if plaintext.Level() > 0 {

		encoder.params.RingQ().PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

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

			res[i] = scaleDown(encoder.bigintCoeffs[i], plaintext.Scale)
		}
		// We can directly get the coefficients
	} else {

		Q := encoder.params.RingQ().Modulus[0]
		coeffs := encoder.polypool.Coeffs[0]

		for i := range res {

			if coeffs[i] >= Q>>1 {
				res[i] = -float64(Q - coeffs[i])
			} else {
				res[i] = float64(coeffs[i])
			}

			res[i] /= plaintext.Scale
		}
	}

	return
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
func NewEncoderBigComplex(params Parameters, logPrecision int) EncoderBigComplex {
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

func (encoder *encoderBigComplex) EncodeNew(values []*ring.Complex, level int, scale float64, logSlots int) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, level, scale)
	encoder.Encode(plaintext, values, logSlots)
	return
}

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

	gap := (encoder.params.RingQ().N >> 1) / slots

	for i, jdx, idx := 0, (encoder.params.RingQ().N >> 1), 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx].Set(encoder.values[i].Real())
		encoder.valuesfloat[jdx].Set(encoder.values[i].Imag())
	}

	scaleUpVecExactBigFloat(encoder.valuesfloat, plaintext.Scale, encoder.params.RingQ().Modulus[:plaintext.Level()+1], plaintext.Value.Coeffs)

	coeffsBigInt := make([]*big.Int, encoder.params.N())

	encoder.params.RingQ().PolyToBigint(plaintext.Value, coeffsBigInt)

	for i := 0; i < (encoder.params.RingQ().N >> 1); i++ {
		encoder.values[i].Real().Set(encoder.zero)
		encoder.values[i].Imag().Set(encoder.zero)
	}

	for i := 0; i < encoder.params.RingQ().N; i++ {
		encoder.valuesfloat[i].Set(encoder.zero)
	}

	encoder.params.RingQ().NTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
	plaintext.Value.IsNTT = true
}

func (encoder *encoderBigComplex) DecodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []*ring.Complex) {
	return encoder.decodePublic(plaintext, logSlots, sigma)
}

func (encoder *encoderBigComplex) Decode(plaintext *Plaintext, logSlots int) (res []*ring.Complex) {
	return encoder.decodePublic(plaintext, logSlots, 0)
}

func (encoder *encoderBigComplex) decodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []*ring.Complex) {

	slots := 1 << logSlots

	if logSlots > encoder.params.LogN()-1 {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	encoder.params.RingQ().InvNTTLvl(plaintext.Level(), plaintext.Value, encoder.polypool)

	if sigma != 0 {
		// B = floor(sigma * sqrt(2*pi))
		encoder.gaussianSampler.ReadAndAddFromDistLvl(plaintext.Level(), encoder.polypool, encoder.params.RingQ(), sigma, int(2.5066282746310002*sigma+0.5))
	}

	encoder.params.RingQ().PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

	Q := encoder.bigintChain[plaintext.Level()]

	maxSlots := encoder.params.RingQ().N >> 1

	scaleFlo := ring.NewFloat(plaintext.Scale, encoder.logPrecision)

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

	SliceBitReverseInPlaceRingComplex(values, N)
}

func (encoder *encoderBigComplex) FFT(values []*ring.Complex, N int) {

	var lenh, lenq, gap, idx int

	u := ring.NewComplex(nil, nil)
	v := ring.NewComplex(nil, nil)

	SliceBitReverseInPlaceRingComplex(values, N)

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
