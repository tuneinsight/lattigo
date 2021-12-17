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

// Encoder is an interface implementing the encoding and decoding operations. It provides methods to encode/decode []complex128 and []float64 types
// into/from Plaintext types.
// Two different encodings are provided:
//
//     - Coeffs: The coefficients are directly embedded on the plaintext. This encoding only allows to encode []float64 slices, but of size up to N
//               (N being the ring degree) and does not preserve the point-wise multiplication. A ciphertext multiplication will result in a nega-
//               cyclic polynomial convolution in the plaintext domain. This encoding does not provide native slot cyclic rotation.
//               Other operations, like addition or constant multiplication behave as usual.
//
//     - Slots: The coefficients are first subjected to a special Fourier transform before being embedded in the plaintext by using Coeffs encoding.
//              This encoding can embed []complex128 and []float64 slices of size at most N/2 (N being the ring degree) and leverages the convolution
//              property of the DFT to preserve point-wise complex multiplication in the plaintext domain, i.e. a ciphertext multiplication will result
//              in an element-wise multiplication in the plaintext domain. It also enables plaintext slots cyclic rotations. Other operations, like
//              constant multiplication behave as usual. It is condidered the default encoding method for CKKS.
//
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
	Encode(values interface{}, plaintext *Plaintext, logSlots int)

	// EncodeNew encodes a set of values on a new plaintext.
	// This method is identical to "EncodeSlotsNew".
	// Encoding is done at the provided level and with the provided scale.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
	// values.(type) can be either []complex128 of []float64.
	// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
	// Returned plaintext is always in the NTT domain.
	EncodeNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext)

	// EncodeSlots encodes a set of values on the target plaintext.
	// Encoding is done at the level and scale of the plaintext.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
	// values.(type) can be either []complex128 of []float64.
	// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
	// Returned plaintext is always in the NTT domain.
	EncodeSlots(values interface{}, plaintext *Plaintext, logSlots int)

	// EncodeSlotsNew encodes a set of values on a new plaintext.
	// Encoding is done at the provided level and with the provided scale.
	// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
	// values.(type) can be either []complex128 of []float64.
	// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
	// Returned plaintext is always in the NTT domain.
	EncodeSlotsNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext)

	// Decode decodes the input plaintext on a new slice of complex128.
	Decode(plaintext *Plaintext, logSlots int) (res []complex128)

	// DecodeSlots decodes the input plaintext on a new slice of complex128.
	DecodeSlots(plaintext *Plaintext, logSlots int) (res []complex128)

	// DecodePublic decodes the input plaintext on a new slice of complex128.
	// Adds, before the decoding step, an error with standard deviation sigma.
	// If the underlying ringType is ConjugateInvariant, the imaginary part (and
	// its related error) are zero.
	DecodePublic(plaintext *Plaintext, logSlots int, sigma float64) []complex128

	// DecodeSlotsPublic decodes the input plaintext on a new slice of complex128.
	// Adds, before the decoding step, an error with standard deviation sigma.
	// If the underlying ringType is ConjugateInvariant, the imaginary part (and
	// its related error) are zero.
	DecodeSlotsPublic(plaintext *Plaintext, logSlots int, sigma float64) []complex128

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
	Encode(values []*ring.Complex, plaintext *Plaintext, logSlots int)

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
	encoder.Encode(values, plaintext, logSlots)
	return
}

func (encoder *encoderComplex128) Encode(values interface{}, plaintext *Plaintext, logSlots int) {

	ringQ := encoder.params.RingQ()

	if logSlots == encoder.params.MaxLogSlots() {
		encoder.embed(values, logSlots, plaintext.Scale, plaintext.Value)
		ringQ.NTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
	} else {

		encoder.embed(values, logSlots, plaintext.Scale, encoder.polypool)

		var n int
		var NTT func(coeffsIn, coeffsOut []uint64, N int, nttPsi []uint64, Q, QInv uint64, bredParams []uint64)
		switch encoder.params.RingType() {
		case ring.Standard:
			n = 2 << logSlots
			NTT = ring.NTT
		case ring.ConjugateInvariant:
			n = 1 << logSlots
			NTT = ring.NTTConjugateInvariant
		}

		N := encoder.params.N()
		gap := N / n
		for i := 0; i < plaintext.Level()+1; i++ {
			// NTT in dimension n
			NTT(encoder.polypool.Coeffs[i][:n], encoder.polypool.Coeffs[i][:n], n, ringQ.NttPsi[i], ringQ.Modulus[i], ringQ.MredParams[i], ringQ.BredParams[i])
			// Maps NTT in dimension n to NTT in dimension N
			tmp0 := encoder.polypool.Coeffs[i]
			tmp1 := plaintext.Value.Coeffs[i]
			for j := 0; j < n; j++ {
				coeff := tmp0[j]
				for w := 0; w < gap; w++ {
					tmp1[j*gap+w] = coeff
				}
			}
		}
	}

	plaintext.Value.IsNTT = true
}

func (encoder *encoderComplex128) EncodeSlotsNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext) {
	return encoder.EncodeNew(values, level, scale, logSlots)
}

func (encoder *encoderComplex128) EncodeSlots(values interface{}, plaintext *Plaintext, logSlots int) {
	encoder.Encode(values, plaintext, logSlots)
}

func (encoder *encoderComplex128) embed(values interface{}, logSlots int, scale float64, polyOut interface{}) {

	slots := 1 << logSlots

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

			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				encoder.valuesFloat[i] = real(encoder.values[i])
				encoder.valuesFloat[j] = imag(encoder.values[i])
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

			for i := 0; i < slots; i++ {
				encoder.valuesFloat[i] = real(encoder.values[i])
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

		for i := 0; i < slots; i++ {
			encoder.valuesFloat[i] = real(encoder.values[i])
		}

		if encoder.params.RingType() == ring.Standard {
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				encoder.valuesFloat[j] = imag(encoder.values[i])
			}
		}

	default:
		panic("values must be []complex128 or []float64")
	}

	switch encoder.params.RingType() {
	case ring.Standard:
		encoder.scaleUp(encoder.valuesFloat[:2*slots], scale, polyOut)
	case ring.ConjugateInvariant:
		encoder.scaleUp(encoder.valuesFloat[:slots], scale, polyOut)
	default:
		panic("invalid ring type")
	}

}

func (encoder *encoderComplex128) scaleUp(values []float64, scale float64, polyOut interface{}) {
	switch p := polyOut.(type) {
	case rlwe.PolyQP:
		levelQ := p.Q.Level()
		levelP := p.P.Level()
		ringQP := encoder.params.RingQP()
		scaleUpVecExact(values, scale, ringQP.RingQ.Modulus[:levelQ+1], p.Q.Coeffs)
		scaleUpVecExact(values, scale, ringQP.RingP.Modulus[:levelP+1], p.P.Coeffs)
	case *ring.Poly:
		scaleUpVecExact(values, scale, encoder.params.RingQ().Modulus[:p.Level()+1], p.Coeffs)
	default:
		panic("invalid polyOut type")
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
	return encoder.DecodeSlotsPublic(plaintext, logSlots, bound)
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoderComplex128) Decode(plaintext *Plaintext, logSlots int) (res []complex128) {
	return encoder.DecodeSlotsPublic(plaintext, logSlots, 0)
}

// DecodeSlotsPublic decodes the Plaintext values to a slice of complex128 values of size at most N/2.
// Adds a Gaussian error to the plaintext of variance sigma and bound floor(sqrt(2*pi)*sigma) before decoding
func (encoder *encoderComplex128) DecodeSlotsPublic(plaintext *Plaintext, logSlots int, bound float64) (res []complex128) {
	return encoder.decodePublic(plaintext, logSlots, bound)
}

// DecodeSlots decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoderComplex128) DecodeSlots(plaintext *Plaintext, logSlots int) (res []complex128) {
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
	encoder.Encode(values, plaintext, logSlots)
	return
}

func (encoder *encoderBigComplex) Encode(values []*ring.Complex, plaintext *Plaintext, logSlots int) {

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
