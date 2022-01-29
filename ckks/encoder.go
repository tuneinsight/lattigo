//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.package ckks
package ckks

import (
	"math"
	"math/big"
	"math/bits"
	"unsafe"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

var pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"

// Encoder is an interface that implements the encoding and decoding operations. It provides methods to encode/decode []complex128 and []float64 types
// into/from Plaintext types.
// Two different encodings are provided:
//
//     - Coeffs: The coefficients are directly embedded on the plaintext. This encoding only allows to encode []float64 slices, but of size up to N
//               (N being the ring degree) and does not preserve the point-wise multiplication. A ciphertext multiplication will result in a nega-
//               cyclic polynomial convolution in the plaintext domain. This encoding does not provide native slot cyclic rotation.
//               Other operations, like addition or constant multiplication, behave as usual.
//
//     - Slots: The coefficients are first subjected to a special Fourier transform before being embedded in the plaintext by using Coeffs encoding.
//              This encoding can embed []complex128 and []float64 slices of size at most N/2 (N being the ring degree) and leverages the convolution
//              property of the DFT to preserve point-wise complex multiplication in the plaintext domain, i.e. a ciphertext multiplication will result
//              in an element-wise multiplication in the plaintext domain. It also enables cyclic rotations on plaintext slots. Other operations, like
//              constant multiplication, behave as usual. It is considered the default encoding method for CKKS.
//
//
// The figure bellow illustrates the relationship between these two encodings:
//
//                                              Real^{N}          Z_Q[X]/(X^N+1)
// EncodeCoeffs: ----------------------------->[]float64 ---------> Plaintext
//                                                 |
//                    Complex^{N/2}                |
// EncodeSlots:  []complex128/[]float64 -> iDFT ---┘
type Encoder interface {

	// Slots Encoding
	Encode(values interface{}, plaintext *Plaintext, logSlots int)
	EncodeNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext)
	EncodeSlots(values interface{}, plaintext *Plaintext, logSlots int)
	EncodeSlotsNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext)
	Decode(plaintext *Plaintext, logSlots int) (res []complex128)
	DecodeSlots(plaintext *Plaintext, logSlots int) (res []complex128)
	DecodePublic(plaintext *Plaintext, logSlots int, sigma float64) []complex128
	DecodeSlotsPublic(plaintext *Plaintext, logSlots int, sigma float64) []complex128

	// Coeffs Encoding
	EncodeCoeffs(values []float64, plaintext *Plaintext)
	EncodeCoeffsNew(values []float64, level int, scale float64) (plaintext *Plaintext)
	DecodeCoeffs(plaintext *Plaintext) (res []float64)
	DecodeCoeffsPublic(plaintext *Plaintext, bound float64) (res []float64)

	// Utility
	GetErrSTDCoeffDomain(valuesWant, valuesHave []complex128, scale float64) (std float64)
	GetErrSTDSlotDomain(valuesWant, valuesHave []complex128, scale float64) (std float64)
	ShallowCopy() Encoder
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

// ShallowCopy creates a shallow copy of encoder in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// encoder can be used concurrently.
func (ecd *encoder) ShallowCopy() *encoder {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &encoder{
		params:          ecd.params,
		bigintChain:     ecd.bigintChain,
		bigintCoeffs:    make([]*big.Int, ecd.m>>1),
		qHalf:           ring.NewUint(0),
		polypool:        ecd.params.RingQ().NewPoly(),
		m:               ecd.m,
		rotGroup:        ecd.rotGroup,
		gaussianSampler: ring.NewGaussianSampler(prng, ecd.params.RingQ(), ecd.params.Sigma(), int(6*ecd.params.Sigma())),
	}
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

	ecd := newEncoder(params)

	var angle float64
	roots := make([]complex128, ecd.m+1)
	for i := 0; i < ecd.m; i++ {
		angle = 2 * 3.141592653589793 * float64(i) / float64(ecd.m)

		roots[i] = complex(math.Cos(angle), math.Sin(angle))
	}
	roots[ecd.m] = roots[0]

	return &encoderComplex128{
		encoder:     ecd,
		roots:       roots,
		values:      make([]complex128, ecd.m>>2),
		valuesFloat: make([]float64, ecd.m>>1),
	}
}

// Encode encodes a set of values on the target plaintext.
// This method is identical to "EncodeSlots".
// Encoding is done at the level and scale of the plaintext.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
// values.(type) can be either []complex128 of []float64.
// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
// Returned plaintext is always in the NTT domain.
func (ecd *encoderComplex128) Encode(values interface{}, plaintext *Plaintext, logSlots int) {
	ecd.embed(values, logSlots, plaintext.Scale, plaintext.Value)
	ecd.switchToNTTDomain(logSlots, false, plaintext.Value)
}

// EncodeNew encodes a set of values on a new plaintext.
// This method is identical to "EncodeSlotsNew".
// Encoding is done at the provided level and with the provided scale.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
// values.(type) can be either []complex128 of []float64.
// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
// Returned plaintext is always in the NTT domain.
func (ecd *encoderComplex128) EncodeNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext) {
	plaintext = NewPlaintext(ecd.params, level, scale)
	ecd.Encode(values, plaintext, logSlots)
	return
}

// EncodeSlots encodes a set of values on the target plaintext.
// Encoding is done at the level and scale of the plaintext.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
// values.(type) can be either []complex128 of []float64.
// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
// Returned plaintext is always in the NTT domain.
func (ecd *encoderComplex128) EncodeSlots(values interface{}, plaintext *Plaintext, logSlots int) {
	ecd.Encode(values, plaintext, logSlots)
}

// EncodeSlotsNew encodes a set of values on a new plaintext.
// Encoding is done at the provided level and with the provided scale.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN.
// values.(type) can be either []complex128 of []float64.
// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
// Returned plaintext is always in the NTT domain.
func (ecd *encoderComplex128) EncodeSlotsNew(values interface{}, level int, scale float64, logSlots int) (plaintext *Plaintext) {
	return ecd.EncodeNew(values, level, scale, logSlots)
}

// Decode decodes the input plaintext on a new slice of complex128.
// This method is the same as .DecodeSlots(*).
func (ecd *encoderComplex128) Decode(plaintext *Plaintext, logSlots int) (res []complex128) {
	return ecd.DecodeSlotsPublic(plaintext, logSlots, 0)
}

// DecodeSlots decodes the input plaintext on a new slice of complex128.
func (ecd *encoderComplex128) DecodeSlots(plaintext *Plaintext, logSlots int) (res []complex128) {
	return ecd.decodePublic(plaintext, logSlots, 0)
}

// DecodePublic decodes the input plaintext on a new slice of complex128.
// This method is the same as .DecodeSlotsPublic(*).
// Adds, before the decoding step, an error with standard deviation sigma and bound floor(sqrt(2*pi)*sigma).
// If the underlying ringType is ConjugateInvariant, the imaginary part (and
// its related error) are zero.
func (ecd *encoderComplex128) DecodePublic(plaintext *Plaintext, logSlots int, bound float64) (res []complex128) {
	return ecd.DecodeSlotsPublic(plaintext, logSlots, bound)
}

// DecodeSlotsPublic decodes the input plaintext on a new slice of complex128.
// Adds, before the decoding step, an error with standard deviation sigma and bound floor(sqrt(2*pi)*sigma).
// If the underlying ringType is ConjugateInvariant, the imaginary part (and
// its related error) are zero.
func (ecd *encoderComplex128) DecodeSlotsPublic(plaintext *Plaintext, logSlots int, bound float64) (res []complex128) {
	return ecd.decodePublic(plaintext, logSlots, bound)
}

// EncodeCoeffs encodes the values on the coefficient of the plaintext polynomial.
// Encoding is done at the level and scale of the plaintext.
// User must ensure that 1<= len(values) <= 2^LogN
func (ecd *encoderComplex128) EncodeCoeffs(values []float64, plaintext *Plaintext) {

	if len(values) > ecd.params.N() {
		panic("cannot EncodeCoeffs : too many values (maximum is N)")
	}
	floatToFixedPointCRT(plaintext.Level(), values, plaintext.Scale, ecd.params.RingQ(), plaintext.Value.Coeffs)
	ecd.params.RingQ().NTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
	plaintext.Value.IsNTT = true
}

// EncodeCoeffsNew encodes the values on the coefficient of a new plaintext.
// Encoding is done at the provided level and with the provided scale.
// User must ensure that 1<= len(values) <= 2^LogN
func (ecd *encoderComplex128) EncodeCoeffsNew(values []float64, level int, scale float64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(ecd.params, level, scale)
	ecd.EncodeCoeffs(values, plaintext)
	return
}

// DecodeCoeffs reconstructs the RNS coefficients of the plaintext on a slice of float64.
func (ecd *encoderComplex128) DecodeCoeffs(plaintext *Plaintext) (res []float64) {
	return ecd.decodeCoeffsPublic(plaintext, 0)
}

// DecodeCoeffsPublic reconstructs the RNS coefficients of the plaintext on a slice of float64.
// Adds an error with standard deviation sigma and bound floor(sqrt(2*pi)*sigma).
func (ecd *encoderComplex128) DecodeCoeffsPublic(plaintext *Plaintext, sigma float64) (res []float64) {
	return ecd.decodeCoeffsPublic(plaintext, sigma)
}

// GetErrSTDCoeffDomain returns StandardDeviation(Encode(valuesWant-valuesHave))*scale
// which is the scaled standard deviation in the coefficient domain of the difference
// of two complex vector in the slot domain.
func (ecd *encoderComplex128) GetErrSTDCoeffDomain(valuesWant, valuesHave []complex128, scale float64) (std float64) {

	for i := range valuesHave {
		ecd.values[i] = (valuesWant[i] - valuesHave[i])
	}

	for i := len(valuesHave); i < len(ecd.values); i++ {
		ecd.values[i] = complex(0, 0)
	}

	// Runs FFT^-1 with the smallest power of two length that is greater than the input size
	invfft(ecd.values, 1<<bits.Len64(uint64(len(valuesHave)-1)), ecd.m, ecd.rotGroup, ecd.roots)

	for i := range valuesWant {
		ecd.valuesFloat[2*i] = real(ecd.values[i])
		ecd.valuesFloat[2*i+1] = imag(ecd.values[i])
	}

	return StandardDeviation(ecd.valuesFloat[:len(valuesWant)*2], scale)
}

// GetErrSTDSlotDomain returns StandardDeviation(valuesWant-valuesHave)*scale
// which is the scaled standard deviation of two complex vectors.
func (ecd *encoderComplex128) GetErrSTDSlotDomain(valuesWant, valuesHave []complex128, scale float64) (std float64) {
	var err complex128
	for i := range valuesWant {
		err = valuesWant[i] - valuesHave[i]
		ecd.valuesFloat[2*i] = real(err)
		ecd.valuesFloat[2*i+1] = imag(err)
	}
	return StandardDeviation(ecd.valuesFloat[:len(valuesWant)*2], scale)
}

// ShallowCopy creates a shallow copy of this encoderComplex128 in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encoder can be used concurrently.
func (ecd *encoderComplex128) ShallowCopy() Encoder {
	return &encoderComplex128{
		encoder:     *ecd.encoder.ShallowCopy(),
		values:      make([]complex128, len(ecd.values)),
		valuesFloat: make([]float64, len(ecd.valuesFloat)),
		roots:       ecd.roots,
	}
}

func (ecd *encoderComplex128) embed(values interface{}, logSlots int, scale float64, polyOut interface{}) {

	slots := 1 << logSlots
	var lenValues int

	// First checks the type of input values
	switch values := values.(type) {

	// If complex
	case []complex128:
		// Checks that the number of values is with the possible range
		if len(values) > int(ecd.params.RingQ().NthRoot>>1) || len(values) > slots || slots > int(ecd.params.RingQ().NthRoot>>2) {
			panic("cannot Encode: too many values/slots for the given ring degree")
		}

		lenValues = len(values)

		switch ecd.params.RingType() {

		case ring.Standard:
			copy(ecd.values[:len(values)], values)

		case ring.ConjugateInvariant:
			// Discards the imaginary part
			for i := range values {
				ecd.values[i] = complex(real(values[i]), 0)
			}

		// Else panics
		default:
			panic("unsuported ringType")
		}

	// If floats only
	case []float64:
		if len(values) > int(ecd.params.RingQ().NthRoot>>1) || len(values) > slots || slots > int(ecd.params.RingQ().NthRoot>>2) {
			panic("cannot Encode: too many values/slots for the given ring degree")
		}

		lenValues = len(values)

		for i := range values {
			ecd.values[i] = complex(values[i], 0)
		}

	default:
		panic("values must be []complex128 or []float64")
	}

	for i := lenValues; i < slots; i += 8 {
		v := (*[8]complex128)(unsafe.Pointer(&ecd.values[i]))
		v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7] = 0, 0, 0, 0, 0, 0, 0, 0
	}

	invfft(ecd.values, slots, ecd.m, ecd.rotGroup, ecd.roots)

	isRingStandard := ecd.params.RingType() == ring.Standard

	switch p := polyOut.(type) {
	case rlwe.PolyQP:
		levelQ := p.Q.Level()
		levelP := p.P.Level()
		complexToFixedPointCRT(levelQ, ecd.values[:slots], scale, ecd.params.RingQ(), p.Q.Coeffs, isRingStandard)
		complexToFixedPointCRT(levelP, ecd.values[:slots], scale, ecd.params.RingP(), p.P.Coeffs, isRingStandard)
	case *ring.Poly:
		complexToFixedPointCRT(p.Level(), ecd.values[:slots], scale, ecd.params.RingQ(), p.Coeffs, isRingStandard)
	default:
		panic("invalid polyOut type")
	}
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

	maxSlots := int(ringQ.NthRoot >> 2)
	slots := 1 << logSlots
	gap := maxSlots / slots

	ringQ.PolyToBigint(poly, gap, bigintCoeffs)

	qHalf := new(big.Int)
	qHalf.Set(Q)
	qHalf.Rsh(qHalf, 1)

	var sign int

	for i := 0; i < slots; i++ {

		// Centers the value around the current modulus
		bigintCoeffs[i].Mod(bigintCoeffs[i], Q)
		sign = bigintCoeffs[i].Cmp(qHalf)
		if sign == 1 || sign == 0 {
			bigintCoeffs[i].Sub(bigintCoeffs[i], Q)
		}

		values[i] = complex(scaleDown(bigintCoeffs[i], scale), 0)
	}

	if !isreal {
		for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
			// Centers the value around the current modulus
			bigintCoeffs[j].Mod(bigintCoeffs[j], Q)
			sign = bigintCoeffs[j].Cmp(qHalf)
			if sign == 1 || sign == 0 {
				bigintCoeffs[j].Sub(bigintCoeffs[j], Q)
			}
			values[i] += complex(0, scaleDown(bigintCoeffs[j], scale))
		}
	}
}

func (ecd *encoderComplex128) plaintextToComplex(level int, scale float64, logSlots int, p *ring.Poly, values []complex128) {

	isreal := ecd.params.RingType() == ring.ConjugateInvariant
	if level == 0 {
		polyToComplexNoCRT(p.Coeffs[0], values, scale, logSlots, isreal, ecd.params.RingQ())
	} else {
		polyToComplexCRT(p, ecd.bigintCoeffs, values, scale, logSlots, isreal, ecd.params.RingQ(), ecd.bigintChain[level])
	}

	if isreal { // [X]/(X^N+1) to [X+X^-1]/(X^N+1)
		tmp := ecd.values
		slots := 1 << logSlots
		for i := 1; i < slots; i++ {
			tmp[i] -= complex(0, real(tmp[slots-i]))
		}
	}
}

func (ecd *encoderComplex128) decodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []complex128) {

	slots := 1 << logSlots

	if slots > int(ecd.params.RingQ().NthRoot>>2) {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	if plaintext.Value.IsNTT {
		ecd.params.RingQ().InvNTTLvl(plaintext.Level(), plaintext.Value, ecd.polypool)
	} else {
		ring.CopyValuesLvl(plaintext.Level(), plaintext.Value, ecd.polypool)
	}

	// B = floor(sigma * sqrt(2*pi))
	if sigma != 0 {
		ecd.gaussianSampler.ReadAndAddFromDistLvl(plaintext.Level(), ecd.polypool, ecd.params.RingQ(), sigma, int(2.5066282746310002*sigma))
	}

	ecd.plaintextToComplex(plaintext.Level(), plaintext.Scale, logSlots, ecd.polypool, ecd.values)

	fft(ecd.values, slots, ecd.m, ecd.rotGroup, ecd.roots)

	res = make([]complex128, slots)

	for i := range res {
		res[i] = ecd.values[i]
	}

	for i := range ecd.values {
		ecd.values[i] = 0
	}

	return
}

func invfft(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	var lenh, lenq, mask, logGap int

	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1

	for loglen := logN; loglen > 0; loglen-- {
		len := 1 << loglen
		lenh = len >> 1
		lenq = len << 2
		logGap = logM - 2 - loglen
		mask = lenq - 1

		if lenh > 8 {
			for i := 0; i < N; i += len {
				for j, k := 0, i; j < lenh; j, k = j+8, k+8 {

					u := (*[8]complex128)(unsafe.Pointer(&values[k]))
					v := (*[8]complex128)(unsafe.Pointer(&values[k+lenh]))
					w := (*[8]int)(unsafe.Pointer(&rotGroup[j]))

					u[0], v[0] = u[0]+v[0], (u[0]-v[0])*roots[(lenq-(w[0]&mask))<<logGap]
					u[1], v[1] = u[1]+v[1], (u[1]-v[1])*roots[(lenq-(w[1]&mask))<<logGap]
					u[2], v[2] = u[2]+v[2], (u[2]-v[2])*roots[(lenq-(w[2]&mask))<<logGap]
					u[3], v[3] = u[3]+v[3], (u[3]-v[3])*roots[(lenq-(w[3]&mask))<<logGap]
					u[4], v[4] = u[4]+v[4], (u[4]-v[4])*roots[(lenq-(w[4]&mask))<<logGap]
					u[5], v[5] = u[5]+v[5], (u[5]-v[5])*roots[(lenq-(w[5]&mask))<<logGap]
					u[6], v[6] = u[6]+v[6], (u[6]-v[6])*roots[(lenq-(w[6]&mask))<<logGap]
					u[7], v[7] = u[7]+v[7], (u[7]-v[7])*roots[(lenq-(w[7]&mask))<<logGap]
				}
			}
		} else if lenh == 8 {

			psi0 := roots[(lenq-(rotGroup[0]&mask))<<logGap]
			psi1 := roots[(lenq-(rotGroup[1]&mask))<<logGap]
			psi2 := roots[(lenq-(rotGroup[2]&mask))<<logGap]
			psi3 := roots[(lenq-(rotGroup[3]&mask))<<logGap]
			psi4 := roots[(lenq-(rotGroup[4]&mask))<<logGap]
			psi5 := roots[(lenq-(rotGroup[5]&mask))<<logGap]
			psi6 := roots[(lenq-(rotGroup[6]&mask))<<logGap]
			psi7 := roots[(lenq-(rotGroup[7]&mask))<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[0], u[8] = u[0]+u[8], (u[0]-u[8])*psi0
				u[1], u[9] = u[1]+u[9], (u[1]-u[9])*psi1
				u[2], u[10] = u[2]+u[10], (u[2]-u[10])*psi2
				u[3], u[11] = u[3]+u[11], (u[3]-u[11])*psi3
				u[4], u[12] = u[4]+u[12], (u[4]-u[12])*psi4
				u[5], u[13] = u[5]+u[13], (u[5]-u[13])*psi5
				u[6], u[14] = u[6]+u[14], (u[6]-u[14])*psi6
				u[7], u[15] = u[7]+u[15], (u[7]-u[15])*psi7
			}

		} else if lenh == 4 {

			psi0 := roots[(lenq-(rotGroup[0]&mask))<<logGap]
			psi1 := roots[(lenq-(rotGroup[1]&mask))<<logGap]
			psi2 := roots[(lenq-(rotGroup[2]&mask))<<logGap]
			psi3 := roots[(lenq-(rotGroup[3]&mask))<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[0], u[4] = u[0]+u[4], (u[0]-u[4])*psi0
				u[1], u[5] = u[1]+u[5], (u[1]-u[5])*psi1
				u[2], u[6] = u[2]+u[6], (u[2]-u[6])*psi2
				u[3], u[7] = u[3]+u[7], (u[3]-u[7])*psi3
				u[8], u[12] = u[8]+u[12], (u[8]-u[12])*psi0
				u[9], u[13] = u[9]+u[13], (u[9]-u[13])*psi1
				u[10], u[14] = u[10]+u[14], (u[10]-u[14])*psi2
				u[11], u[15] = u[11]+u[15], (u[11]-u[15])*psi3
			}
		} else if lenh == 2 {

			psi0 := roots[(lenq-(rotGroup[0]&mask))<<logGap]
			psi1 := roots[(lenq-(rotGroup[1]&mask))<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[0], u[2] = u[0]+u[2], (u[0]-u[2])*psi0
				u[1], u[3] = u[1]+u[3], (u[1]-u[3])*psi1
				u[4], u[6] = u[4]+u[6], (u[4]-u[6])*psi0
				u[5], u[7] = u[5]+u[7], (u[5]-u[7])*psi1
				u[8], u[10] = u[8]+u[10], (u[8]-u[10])*psi0
				u[9], u[11] = u[9]+u[11], (u[9]-u[11])*psi1
				u[12], u[14] = u[12]+u[14], (u[12]-u[14])*psi0
				u[13], u[15] = u[13]+u[15], (u[13]-u[15])*psi1
			}
		} else if lenh == 1 {

			psi0 := roots[(lenq-(rotGroup[0]&mask))<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[0], u[1] = u[0]+u[1], (u[0]-u[1])*psi0
				u[2], u[3] = u[2]+u[3], (u[2]-u[3])*psi0
				u[4], u[5] = u[4]+u[5], (u[4]-u[5])*psi0
				u[6], u[7] = u[6]+u[7], (u[6]-u[7])*psi0
				u[8], u[9] = u[8]+u[9], (u[8]-u[9])*psi0
				u[10], u[11] = u[10]+u[11], (u[10]-u[11])*psi0
				u[12], u[13] = u[12]+u[13], (u[12]-u[13])*psi0
				u[14], u[15] = u[14]+u[15], (u[14]-u[15])*psi0
			}
		}
	}

	invN := complex(float64(N), 0)
	for i := 0; i < N; i = i + 8 {

		v := (*[8]complex128)(unsafe.Pointer(&values[i]))

		v[0] /= invN
		v[1] /= invN
		v[2] /= invN
		v[3] /= invN
		v[4] /= invN
		v[5] /= invN
		v[6] /= invN
		v[7] /= invN
	}

	SliceBitReverseInPlaceComplex128(values, N)
}

func fft(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	SliceBitReverseInPlaceComplex128(values, N)

	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1

	for loglen := 1; loglen <= logN; loglen++ {
		len := 1 << loglen
		lenh := len >> 1
		lenq := len << 2
		logGap := logM - 2 - loglen
		mask := lenq - 1

		if lenh > 8 {
			for i := 0; i < N; i += len {

				for j, k := 0, i; j < lenh; j, k = j+8, k+8 {

					u := (*[8]complex128)(unsafe.Pointer(&values[k]))
					v := (*[8]complex128)(unsafe.Pointer(&values[k+lenh]))
					w := (*[8]int)(unsafe.Pointer(&rotGroup[j]))

					v[0] *= roots[(w[0]&mask)<<logGap]
					v[1] *= roots[(w[1]&mask)<<logGap]
					v[2] *= roots[(w[2]&mask)<<logGap]
					v[3] *= roots[(w[3]&mask)<<logGap]
					v[4] *= roots[(w[4]&mask)<<logGap]
					v[5] *= roots[(w[5]&mask)<<logGap]
					v[6] *= roots[(w[6]&mask)<<logGap]
					v[7] *= roots[(w[7]&mask)<<logGap]

					u[0], v[0] = u[0]+v[0], u[0]-v[0]
					u[1], v[1] = u[1]+v[1], u[1]-v[1]
					u[2], v[2] = u[2]+v[2], u[2]-v[2]
					u[3], v[3] = u[3]+v[3], u[3]-v[3]
					u[4], v[4] = u[4]+v[4], u[4]-v[4]
					u[5], v[5] = u[5]+v[5], u[5]-v[5]
					u[6], v[6] = u[6]+v[6], u[6]-v[6]
					u[7], v[7] = u[7]+v[7], u[7]-v[7]
				}
			}
		} else if lenh == 8 {

			psi0 := roots[(rotGroup[0]&mask)<<logGap]
			psi1 := roots[(rotGroup[1]&mask)<<logGap]
			psi2 := roots[(rotGroup[2]&mask)<<logGap]
			psi3 := roots[(rotGroup[3]&mask)<<logGap]
			psi4 := roots[(rotGroup[4]&mask)<<logGap]
			psi5 := roots[(rotGroup[5]&mask)<<logGap]
			psi6 := roots[(rotGroup[6]&mask)<<logGap]
			psi7 := roots[(rotGroup[7]&mask)<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[8] *= psi0
				u[9] *= psi1
				u[10] *= psi2
				u[11] *= psi3
				u[12] *= psi4
				u[13] *= psi5
				u[14] *= psi6
				u[15] *= psi7

				u[0], u[8] = u[0]+u[8], u[0]-u[8]
				u[1], u[9] = u[1]+u[9], u[1]-u[9]
				u[2], u[10] = u[2]+u[10], u[2]-u[10]
				u[3], u[11] = u[3]+u[11], u[3]-u[11]
				u[4], u[12] = u[4]+u[12], u[4]-u[12]
				u[5], u[13] = u[5]+u[13], u[5]-u[13]
				u[6], u[14] = u[6]+u[14], u[6]-u[14]
				u[7], u[15] = u[7]+u[15], u[7]-u[15]

			}
		} else if lenh == 4 {

			psi0 := roots[(rotGroup[0]&mask)<<logGap]
			psi1 := roots[(rotGroup[1]&mask)<<logGap]
			psi2 := roots[(rotGroup[2]&mask)<<logGap]
			psi3 := roots[(rotGroup[3]&mask)<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[4] *= psi0
				u[5] *= psi1
				u[6] *= psi2
				u[7] *= psi3
				u[12] *= psi0
				u[13] *= psi1
				u[14] *= psi2
				u[15] *= psi3

				u[0], u[4] = u[0]+u[4], u[0]-u[4]
				u[1], u[5] = u[1]+u[5], u[1]-u[5]
				u[2], u[6] = u[2]+u[6], u[2]-u[6]
				u[3], u[7] = u[3]+u[7], u[3]-u[7]
				u[8], u[12] = u[8]+u[12], u[8]-u[12]
				u[9], u[13] = u[9]+u[13], u[9]-u[13]
				u[10], u[14] = u[10]+u[14], u[10]-u[14]
				u[11], u[15] = u[11]+u[15], u[11]-u[15]
			}
		} else if lenh == 2 {

			psi0 := roots[(rotGroup[0]&mask)<<logGap]
			psi1 := roots[(rotGroup[1]&mask)<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[2] *= psi0
				u[3] *= psi1
				u[6] *= psi0
				u[7] *= psi1
				u[10] *= psi0
				u[11] *= psi1
				u[14] *= psi0
				u[15] *= psi1

				u[0], u[2] = u[0]+u[2], u[0]-u[2]
				u[1], u[3] = u[1]+u[3], u[1]-u[3]
				u[4], u[6] = u[4]+u[6], u[4]-u[6]
				u[5], u[7] = u[5]+u[7], u[5]-u[7]
				u[8], u[10] = u[8]+u[10], u[8]-u[10]
				u[9], u[11] = u[9]+u[11], u[9]-u[11]
				u[12], u[14] = u[12]+u[14], u[12]-u[14]
				u[13], u[15] = u[13]+u[15], u[13]-u[15]
			}
		} else if lenh == 1 {

			psi0 := roots[(rotGroup[0]&mask)<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[1] *= psi0
				u[3] *= psi0
				u[5] *= psi0
				u[7] *= psi0
				u[9] *= psi0
				u[11] *= psi0
				u[13] *= psi0
				u[15] *= psi0

				u[0], u[1] = u[0]+u[1], u[0]-u[1]
				u[2], u[3] = u[2]+u[3], u[2]-u[3]
				u[4], u[5] = u[4]+u[5], u[4]-u[5]
				u[6], u[7] = u[6]+u[7], u[6]-u[7]
				u[8], u[9] = u[8]+u[9], u[8]-u[9]
				u[10], u[11] = u[10]+u[11], u[10]-u[11]
				u[12], u[13] = u[12]+u[13], u[12]-u[13]
				u[14], u[15] = u[14]+u[15], u[14]-u[15]
			}
		}
	}
}

func (ecd *encoderComplex128) switchToNTTDomain(logSlots int, montgomery bool, polyOut interface{}) {
	switch p := polyOut.(type) {
	case rlwe.PolyQP:
		NttAndMontgomeryLvl(p.Q.Level(), logSlots, ecd.params.RingQ(), montgomery, p.Q)
		NttAndMontgomeryLvl(p.P.Level(), logSlots, ecd.params.RingP(), montgomery, p.P)
	case *ring.Poly:
		NttAndMontgomeryLvl(p.Level(), logSlots, ecd.params.RingQ(), montgomery, p)
	default:
		panic("invalid polyOut type")
	}
}

func (ecd *encoderComplex128) decodeCoeffsPublic(plaintext *Plaintext, sigma float64) (res []float64) {

	if plaintext.Value.IsNTT {
		ecd.params.RingQ().InvNTTLvl(plaintext.Level(), plaintext.Value, ecd.polypool)
	} else {
		ring.CopyValuesLvl(plaintext.Level(), plaintext.Value, ecd.polypool)
	}

	if sigma != 0 {
		// B = floor(sigma * sqrt(2*pi))
		ecd.gaussianSampler.ReadAndAddFromDistLvl(plaintext.Level(), ecd.polypool, ecd.params.RingQ(), sigma, int(2.5066282746310002*sigma))
	}

	res = make([]float64, ecd.params.N())

	// We have more than one moduli and need the CRT reconstruction
	if plaintext.Level() > 0 {

		ecd.params.RingQ().PolyToBigint(ecd.polypool, 1, ecd.bigintCoeffs)

		Q := ecd.bigintChain[plaintext.Level()]

		ecd.qHalf.Set(Q)
		ecd.qHalf.Rsh(ecd.qHalf, 1)

		var sign int

		for i := range res {

			// Centers the value around the current modulus
			ecd.bigintCoeffs[i].Mod(ecd.bigintCoeffs[i], Q)

			sign = ecd.bigintCoeffs[i].Cmp(ecd.qHalf)
			if sign == 1 || sign == 0 {
				ecd.bigintCoeffs[i].Sub(ecd.bigintCoeffs[i], Q)
			}

			res[i] = scaleDown(ecd.bigintCoeffs[i], plaintext.Scale)
		}
		// We can directly get the coefficients
	} else {

		Q := ecd.params.RingQ().Modulus[0]
		coeffs := ecd.polypool.Coeffs[0]

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

// EncoderBigComplex is an interface that implements the encoding algorithms with arbitrary precision.
type EncoderBigComplex interface {
	Encode(values []*ring.Complex, plaintext *Plaintext, logSlots int)
	EncodeNew(values []*ring.Complex, level int, scale float64, logSlots int) (plaintext *Plaintext)
	Decode(plaintext *Plaintext, logSlots int) (res []*ring.Complex)
	DecodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []*ring.Complex)
	FFT(values []*ring.Complex, N int)
	InvFFT(values []*ring.Complex, N int)
	ShallowCopy() EncoderBigComplex
}

type encoderBigComplex struct {
	encoder
	zero         *big.Float
	cMul         *ring.ComplexMultiplier
	logPrecision int
	values       []*ring.Complex
	valuesfloat  []*big.Float
	roots        []*ring.Complex
}

// NewEncoderBigComplex creates a new encoder using arbitrary precision complex arithmetic.
func NewEncoderBigComplex(params Parameters, logPrecision int) EncoderBigComplex {
	ecd := newEncoder(params)

	var PI = new(big.Float)
	PI.SetPrec(uint(logPrecision))
	PI.SetString(pi)

	var PIHalf = new(big.Float)
	PIHalf.SetPrec(uint(logPrecision))
	PIHalf.SetString(pi)
	PIHalf.Quo(PIHalf, ring.NewFloat(2, logPrecision))

	var angle *big.Float
	roots := make([]*ring.Complex, ecd.m+1)
	for i := 0; i < ecd.m; i++ {
		angle = ring.NewFloat(2, logPrecision)
		angle.Mul(angle, PI)
		angle.Mul(angle, ring.NewFloat(float64(i), logPrecision))
		angle.Quo(angle, ring.NewFloat(float64(ecd.m), logPrecision))

		real := ring.Cos(angle)
		angle.Sub(PIHalf, angle)
		imag := ring.Cos(angle)

		roots[i] = ring.NewComplex(real, imag)
	}

	roots[ecd.m] = roots[0].Copy()

	values := make([]*ring.Complex, ecd.m>>2)
	valuesfloat := make([]*big.Float, ecd.m>>1)

	for i := 0; i < ecd.m>>2; i++ {

		values[i] = ring.NewComplex(ring.NewFloat(0, logPrecision), ring.NewFloat(0, logPrecision))
		valuesfloat[i*2] = ring.NewFloat(0, logPrecision)
		valuesfloat[(i*2)+1] = ring.NewFloat(0, logPrecision)
	}

	return &encoderBigComplex{
		encoder:      ecd,
		zero:         ring.NewFloat(0, logPrecision),
		cMul:         ring.NewComplexMultiplier(),
		logPrecision: logPrecision,
		roots:        roots,
		values:       values,
		valuesfloat:  valuesfloat,
	}
}

// Encode encodes a set of values on the target plaintext.
// Encoding is done at the level and scale of the plaintext.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^LogN.
func (ecd *encoderBigComplex) Encode(values []*ring.Complex, plaintext *Plaintext, logSlots int) {

	slots := 1 << logSlots

	if len(values) > ecd.params.N()/2 || len(values) > slots || logSlots > ecd.params.LogN()-1 {
		panic("cannot Encode: too many values/slots for the given ring degree")
	}

	if len(values) != slots {
		panic("cannot Encode: number of values must be equal to slots")
	}

	for i := 0; i < slots; i++ {
		ecd.values[i].Set(values[i])
	}

	ecd.InvFFT(ecd.values, slots)

	gap := (ecd.params.RingQ().N >> 1) / slots

	for i, jdx, idx := 0, (ecd.params.RingQ().N >> 1), 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		ecd.valuesfloat[idx].Set(ecd.values[i].Real())
		ecd.valuesfloat[jdx].Set(ecd.values[i].Imag())
	}

	scaleUpVecExactBigFloat(ecd.valuesfloat, plaintext.Scale, ecd.params.RingQ().Modulus[:plaintext.Level()+1], plaintext.Value.Coeffs)

	for i := 0; i < (ecd.params.RingQ().N >> 1); i++ {
		ecd.values[i].Real().Set(ecd.zero)
		ecd.values[i].Imag().Set(ecd.zero)
	}

	for i := 0; i < ecd.params.RingQ().N; i++ {
		ecd.valuesfloat[i].Set(ecd.zero)
	}

	ecd.params.RingQ().NTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
	plaintext.Value.IsNTT = true
}

// EncodeNew encodes a set of values on a new plaintext.
// Encoding is done at the provided level and with the provided scale.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^LogN.
func (ecd *encoderBigComplex) EncodeNew(values []*ring.Complex, level int, scale float64, logSlots int) (plaintext *Plaintext) {
	plaintext = NewPlaintext(ecd.params, level, scale)
	ecd.Encode(values, plaintext, logSlots)
	return
}

// Decode decodes the input plaintext on a new slice of ring.Complex.
func (ecd *encoderBigComplex) Decode(plaintext *Plaintext, logSlots int) (res []*ring.Complex) {
	return ecd.decodePublic(plaintext, logSlots, 0)
}

func (ecd *encoderBigComplex) DecodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []*ring.Complex) {
	return ecd.decodePublic(plaintext, logSlots, sigma)
}

// FFT evaluates the decoding matrix on a slice of ring.Complex values.
func (ecd *encoderBigComplex) FFT(values []*ring.Complex, N int) {

	var lenh, lenq, gap, idx int

	u := ring.NewComplex(nil, nil)
	v := ring.NewComplex(nil, nil)

	SliceBitReverseInPlaceRingComplex(values, N)

	for len := 2; len <= N; len <<= 1 {
		for i := 0; i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = ecd.m / lenq
			for j := 0; j < lenh; j++ {
				idx = (ecd.rotGroup[j] % lenq) * gap
				u.Set(values[i+j])
				v.Set(values[i+j+lenh])
				ecd.cMul.Mul(v, ecd.roots[idx], v)
				values[i+j].Add(u, v)
				values[i+j+lenh].Sub(u, v)
			}
		}
	}
}

// InvFFT evaluates the encoding matrix on a slice of ring.Complex values.
func (ecd *encoderBigComplex) InvFFT(values []*ring.Complex, N int) {

	var lenh, lenq, gap, idx int
	u := ring.NewComplex(nil, nil)
	v := ring.NewComplex(nil, nil)

	for len := N; len >= 1; len >>= 1 {
		for i := 0; i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = ecd.m / lenq
			for j := 0; j < lenh; j++ {
				idx = (lenq - (ecd.rotGroup[j] % lenq)) * gap
				u.Add(values[i+j], values[i+j+lenh])
				v.Sub(values[i+j], values[i+j+lenh])
				ecd.cMul.Mul(v, ecd.roots[idx], v)
				values[i+j].Set(u)
				values[i+j+lenh].Set(v)
			}
		}
	}

	NBig := ring.NewFloat(float64(N), ecd.logPrecision)
	for i := range values {
		values[i][0].Quo(values[i][0], NBig)
		values[i][1].Quo(values[i][1], NBig)
	}

	SliceBitReverseInPlaceRingComplex(values, N)
}

// ShallowCopy creates a shallow copy of this encoderBigComplex in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// EncoderBigComplex can be used concurrently.
func (ecd *encoderBigComplex) ShallowCopy() EncoderBigComplex {

	values := make([]*ring.Complex, ecd.m>>2)
	valuesfloat := make([]*big.Float, ecd.m>>1)

	for i := 0; i < ecd.m>>2; i++ {

		values[i] = ring.NewComplex(ring.NewFloat(0, ecd.logPrecision), ring.NewFloat(0, ecd.logPrecision))
		valuesfloat[i*2] = ring.NewFloat(0, ecd.logPrecision)
		valuesfloat[(i*2)+1] = ring.NewFloat(0, ecd.logPrecision)
	}

	return &encoderBigComplex{
		encoder:      *ecd.encoder.ShallowCopy(),
		zero:         ring.NewFloat(0, ecd.logPrecision),
		cMul:         ring.NewComplexMultiplier(),
		logPrecision: ecd.logPrecision,
		values:       values,
		valuesfloat:  valuesfloat,
		roots:        ecd.roots,
	}
}

func (ecd *encoderBigComplex) decodePublic(plaintext *Plaintext, logSlots int, sigma float64) (res []*ring.Complex) {

	slots := 1 << logSlots

	if logSlots > ecd.params.LogN()-1 {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	ecd.params.RingQ().InvNTTLvl(plaintext.Level(), plaintext.Value, ecd.polypool)

	if sigma != 0 {
		// B = floor(sigma * sqrt(2*pi))
		ecd.gaussianSampler.ReadAndAddFromDistLvl(plaintext.Level(), ecd.polypool, ecd.params.RingQ(), sigma, int(2.5066282746310002*sigma+0.5))
	}

	Q := ecd.bigintChain[plaintext.Level()]

	maxSlots := ecd.params.RingQ().N >> 1

	scaleFlo := ring.NewFloat(plaintext.Scale, ecd.logPrecision)

	ecd.qHalf.Set(Q)
	ecd.qHalf.Rsh(ecd.qHalf, 1)

	gap := maxSlots / slots

	ecd.params.RingQ().PolyToBigint(ecd.polypool, gap, ecd.bigintCoeffs)

	var sign int

	for i, j := 0, slots; i < slots; i, j = i+1, j+1 {

		// Centers the value around the current modulus
		ecd.bigintCoeffs[i].Mod(ecd.bigintCoeffs[i], Q)
		sign = ecd.bigintCoeffs[i].Cmp(ecd.qHalf)
		if sign == 1 || sign == 0 {
			ecd.bigintCoeffs[i].Sub(ecd.bigintCoeffs[i], Q)
		}

		// Centers the value around the current modulus
		ecd.bigintCoeffs[j].Mod(ecd.bigintCoeffs[j], Q)
		sign = ecd.bigintCoeffs[j].Cmp(ecd.qHalf)
		if sign == 1 || sign == 0 {
			ecd.bigintCoeffs[j].Sub(ecd.bigintCoeffs[j], Q)
		}

		ecd.values[i].Real().SetInt(ecd.bigintCoeffs[i])
		ecd.values[i].Real().Quo(ecd.values[i].Real(), scaleFlo)

		ecd.values[i].Imag().SetInt(ecd.bigintCoeffs[j])
		ecd.values[i].Imag().Quo(ecd.values[i].Imag(), scaleFlo)
	}

	ecd.FFT(ecd.values, slots)

	res = make([]*ring.Complex, slots)

	for i := range res {
		res[i] = ecd.values[i].Copy()
	}

	for i := 0; i < maxSlots; i++ {
		ecd.values[i].Real().Set(ecd.zero)
		ecd.values[i].Imag().Set(ecd.zero)
	}

	return
}
