package ckks

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = ring.GaloisGen

var pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"

// Encoder is an interface that implements the encoding and decoding operations. It provides methods to encode/decode []complex128 and []float64 types
// into/from Plaintext types.
// Two different encodings are provided:
//
//   - Coeffs: The coefficients are directly embedded on the plaintext. This encoding only allows to encode []float64 slices, but of size up to N
//     (N being the ring degree) and does not preserve the point-wise multiplication. A ciphertext multiplication will result in a nega-
//     cyclic polynomial convolution in the plaintext domain. This encoding does not provide native slot cyclic rotation.
//     Other operations, like addition or constant multiplication, behave as usual.
//
//   - Slots: The coefficients are first subjected to a special Fourier transform before being embedded in the plaintext by using Coeffs encoding.
//     This encoding can embed []complex128 and []float64 slices of size at most N/2 (N being the ring degree) and leverages the convolution
//     property of the DFT to preserve point-wise complex multiplication in the plaintext domain, i.e. a ciphertext multiplication will result
//     in an element-wise multiplication in the plaintext domain. It also enables cyclic rotations on plaintext slots. Other operations, like
//     constant multiplication, behave as usual. It is considered the default encoding method for CKKS.
//
// The figure bellow illustrates the relationship between these two encodings:
//
//	                                             Real^{N}          Z_Q[X]/(X^N+1)
//		EncodeCoeffs: ----------------------------->[]float64 ---------> Plaintext
//	                                                |
//	                   Complex^{N/2}                |
//		EncodeSlots:  []complex128/[]float64 -> iDFT ---┘
type Encoder interface {

	// Slots Encoding
	Encode(values interface{}, plaintext *rlwe.Plaintext, logSlots int)
	EncodeNew(values interface{}, level int, scale rlwe.Scale, logSlots int) (plaintext *rlwe.Plaintext)
	EncodeSlots(values interface{}, plaintext *rlwe.Plaintext, logSlots int)
	EncodeSlotsNew(values interface{}, level int, scale rlwe.Scale, logSlots int) (plaintext *rlwe.Plaintext)
	Decode(plaintext *rlwe.Plaintext, logSlots int) (res []complex128)
	DecodeSlots(plaintext *rlwe.Plaintext, logSlots int) (res []complex128)
	DecodePublic(plaintext *rlwe.Plaintext, logSlots int, sigma float64) []complex128
	DecodeSlotsPublic(plaintext *rlwe.Plaintext, logSlots int, sigma float64) []complex128

	FFT(values []complex128, N int)
	IFFT(values []complex128, N int)

	// Coeffs Encoding
	EncodeCoeffs(values []float64, plaintext *rlwe.Plaintext)
	EncodeCoeffsNew(values []float64, level int, scale rlwe.Scale) (plaintext *rlwe.Plaintext)
	DecodeCoeffs(plaintext *rlwe.Plaintext) (res []float64)
	DecodeCoeffsPublic(plaintext *rlwe.Plaintext, bound float64) (res []float64)

	// Utility
	Embed(values interface{}, logSlots int, scale rlwe.Scale, montgomery bool, polyOut interface{})
	GetErrSTDCoeffDomain(valuesWant, valuesHave []complex128, scale rlwe.Scale) (std float64)
	GetErrSTDSlotDomain(valuesWant, valuesHave []complex128, scale rlwe.Scale) (std float64)
	ShallowCopy() Encoder
	Parameters() Parameters
}

// encoder is a struct storing the necessary parameters to encode a slice of complex number on a Plaintext.
type encoder struct {
	params       Parameters
	bigintCoeffs []*big.Int
	qHalf        *big.Int
	buff         *ring.Poly
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

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &encoder{
		params:          ecd.params,
		bigintCoeffs:    make([]*big.Int, ecd.m>>1),
		qHalf:           ring.NewUint(0),
		buff:            ecd.params.RingQ().NewPoly(),
		m:               ecd.m,
		rotGroup:        ecd.rotGroup,
		gaussianSampler: ring.NewGaussianSampler(prng, ecd.params.RingQ(), ecd.params.Sigma(), int(6*ecd.params.Sigma())),
	}
}

// Parameters returns the parameters used by the encoder.
func (ecd *encoder) Parameters() Parameters {
	return ecd.params
}

func newEncoder(params Parameters) encoder {

	m := int(params.RingQ().NthRoot())

	rotGroup := make([]int, m>>2)
	fivePows := 1
	for i := 0; i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= int(GaloisGen)
		fivePows &= (m - 1)
	}

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	gaussianSampler := ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))

	return encoder{
		params:          params,
		bigintCoeffs:    make([]*big.Int, m>>1),
		qHalf:           ring.NewUint(0),
		buff:            params.RingQ().NewPoly(),
		m:               m,
		rotGroup:        rotGroup,
		gaussianSampler: gaussianSampler,
	}
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a Plaintext.
func NewEncoder(params Parameters) Encoder {

	ecd := newEncoder(params)

	return &encoderComplex128{
		encoder:     ecd,
		roots:       GetRootsFloat64(ecd.m),
		values:      make([]complex128, ecd.m>>2),
		valuesFloat: make([]float64, ecd.m>>1),
	}
}

// Encode encodes a set of values on the target plaintext.
// This method is identical to "EncodeSlots".
// Encoding is done at the level and scale of the plaintext.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN and that logSlots >= 3.
// values.(type) can be either []complex128 of []float64.
// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
// Returned plaintext is always in the NTT domain.
func (ecd *encoderComplex128) Encode(values interface{}, plaintext *rlwe.Plaintext, logSlots int) {
	ecd.Embed(values, logSlots, plaintext.Scale, false, plaintext.Value)
}

// EncodeNew encodes a set of values on a new plaintext.
// This method is identical to "EncodeSlotsNew".
// Encoding is done at the provided level and with the provided scale.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN and that logSlots >= 3.
// values.(type) can be either []complex128 of []float64.
// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
// Returned plaintext is always in the NTT domain.
func (ecd *encoderComplex128) EncodeNew(values interface{}, level int, scale rlwe.Scale, logSlots int) (plaintext *rlwe.Plaintext) {
	plaintext = NewPlaintext(ecd.params, level)
	plaintext.Scale = scale
	ecd.Encode(values, plaintext, logSlots)
	return
}

// EncodeSlots encodes a set of values on the target plaintext.
// Encoding is done at the level and scale of the plaintext.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN and that logSlots >= 3.
// values.(type) can be either []complex128 of []float64.
// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
// Returned plaintext is always in the NTT domain.
func (ecd *encoderComplex128) EncodeSlots(values interface{}, plaintext *rlwe.Plaintext, logSlots int) {
	ecd.Encode(values, plaintext, logSlots)
}

// EncodeSlotsNew encodes a set of values on a new plaintext.
// Encoding is done at the provided level and with the provided scale.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^logN and that logSlots >= 3.
// values.(type) can be either []complex128 of []float64.
// The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
// Returned plaintext is always in the NTT domain.
func (ecd *encoderComplex128) EncodeSlotsNew(values interface{}, level int, scale rlwe.Scale, logSlots int) (plaintext *rlwe.Plaintext) {
	return ecd.EncodeNew(values, level, scale, logSlots)
}

// Decode decodes the input plaintext on a new slice of complex128.
// This method is the same as .DecodeSlots(*).
func (ecd *encoderComplex128) Decode(plaintext *rlwe.Plaintext, logSlots int) (res []complex128) {
	return ecd.DecodeSlotsPublic(plaintext, logSlots, 0)
}

// DecodeSlots decodes the input plaintext on a new slice of complex128.
func (ecd *encoderComplex128) DecodeSlots(plaintext *rlwe.Plaintext, logSlots int) (res []complex128) {
	return ecd.decodePublic(plaintext, logSlots, 0)
}

// DecodePublic decodes the input plaintext on a new slice of complex128.
// This method is the same as .DecodeSlotsPublic(*).
// Adds, before the decoding step, an error with standard deviation sigma and bound floor(sqrt(2*pi)*sigma).
// If the underlying ringType is ConjugateInvariant, the imaginary part (and
// its related error) are zero.
func (ecd *encoderComplex128) DecodePublic(plaintext *rlwe.Plaintext, logSlots int, bound float64) (res []complex128) {
	return ecd.DecodeSlotsPublic(plaintext, logSlots, bound)
}

// DecodeSlotsPublic decodes the input plaintext on a new slice of complex128.
// Adds, before the decoding step, an error with standard deviation sigma and bound floor(sqrt(2*pi)*sigma).
// If the underlying ringType is ConjugateInvariant, the imaginary part (and
// its related error) are zero.
func (ecd *encoderComplex128) DecodeSlotsPublic(plaintext *rlwe.Plaintext, logSlots int, bound float64) (res []complex128) {
	return ecd.decodePublic(plaintext, logSlots, bound)
}

// EncodeCoeffs encodes the values on the coefficient of the plaintext polynomial.
// Encoding is done at the level and scale of the plaintext.
// User must ensure that 1<= len(values) <= 2^LogN
func (ecd *encoderComplex128) EncodeCoeffs(values []float64, plaintext *rlwe.Plaintext) {

	if len(values) > ecd.params.N() {
		panic("cannot EncodeCoeffs: too many values (maximum is N)")
	}

	FloatToFixedPointCRT(ecd.params.RingQ().AtLevel(plaintext.Level()), values, plaintext.Scale.Float64(), plaintext.Value.Coeffs)
	ecd.params.RingQ().AtLevel(plaintext.Level()).NTT(plaintext.Value, plaintext.Value)
}

// EncodeCoeffsNew encodes the values on the coefficient of a new plaintext.
// Encoding is done at the provided level and with the provided scale.
// User must ensure that 1<= len(values) <= 2^LogN
func (ecd *encoderComplex128) EncodeCoeffsNew(values []float64, level int, scale rlwe.Scale) (plaintext *rlwe.Plaintext) {
	plaintext = NewPlaintext(ecd.params, level)
	plaintext.Scale = scale
	ecd.EncodeCoeffs(values, plaintext)
	return
}

// DecodeCoeffs reconstructs the RNS coefficients of the plaintext on a slice of float64.
func (ecd *encoderComplex128) DecodeCoeffs(plaintext *rlwe.Plaintext) (res []float64) {
	return ecd.decodeCoeffsPublic(plaintext, 0)
}

// DecodeCoeffsPublic reconstructs the RNS coefficients of the plaintext on a slice of float64.
// Adds an error with standard deviation sigma and bound floor(sqrt(2*pi)*sigma).
func (ecd *encoderComplex128) DecodeCoeffsPublic(plaintext *rlwe.Plaintext, sigma float64) (res []float64) {
	return ecd.decodeCoeffsPublic(plaintext, sigma)
}

// GetErrSTDCoeffDomain returns StandardDeviation(Encode(valuesWant-valuesHave))*scale
// which is the scaled standard deviation in the coefficient domain of the difference
// of two complex vector in the slot domain.
func (ecd *encoderComplex128) GetErrSTDCoeffDomain(valuesWant, valuesHave []complex128, scale rlwe.Scale) (std float64) {

	for i := range valuesHave {
		ecd.values[i] = (valuesWant[i] - valuesHave[i])
	}

	for i := len(valuesHave); i < len(ecd.values); i++ {
		ecd.values[i] = complex(0, 0)
	}

	logSlots := bits.Len64(uint64(len(valuesHave) - 1))

	ecd.IFFT(ecd.values, logSlots)

	for i := range valuesWant {
		ecd.valuesFloat[2*i] = real(ecd.values[i])
		ecd.valuesFloat[2*i+1] = imag(ecd.values[i])
	}

	return StandardDeviation(ecd.valuesFloat[:len(valuesWant)*2], scale.Float64())
}

// GetErrSTDSlotDomain returns StandardDeviation(valuesWant-valuesHave)*scale
// which is the scaled standard deviation of two complex vectors.
func (ecd *encoderComplex128) GetErrSTDSlotDomain(valuesWant, valuesHave []complex128, scale rlwe.Scale) (std float64) {
	var err complex128
	for i := range valuesWant {
		err = valuesWant[i] - valuesHave[i]
		ecd.valuesFloat[2*i] = real(err)
		ecd.valuesFloat[2*i+1] = imag(err)
	}

	return StandardDeviation(ecd.valuesFloat[:len(valuesWant)*2], scale.Float64())
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

// Embed is a generic method to encode a set of values on the target polyOut interface.
// This method it as the core of the slot encoding.
// values: values.(type) can be either []complex128 of []float64.
//
//	The imaginary part of []complex128 will be discarded if ringType == ring.ConjugateInvariant.
//
// logslots: user must ensure that 1 <= len(values) <= 2^logSlots < 2^logN and that logSlots >= 3.
// scale: the scaling factor used do discretize float64 to fixed point integers.
// montgomery: if true then the value written on polyOut are put in the Montgomery domain.
// polyOut: polyOut.(type) can be either ringqp.Poly or *ring.Poly.
//
//	The encoding encoding is done at the level of polyOut.
//
// Values written on  polyOut are always in the NTT domain.
func (ecd *encoderComplex128) Embed(values interface{}, logSlots int, scale rlwe.Scale, montgomery bool, polyOut interface{}) {

	if logSlots < minLogSlots || logSlots > ecd.params.MaxLogSlots() {
		panic(fmt.Sprintf("cannot Embed: logSlots (%d) must be greater or equal to %d and smaller than %d\n", logSlots, minLogSlots, ecd.params.MaxLogSlots()))
	}

	slots := 1 << logSlots
	var lenValues int

	// First checks the type of input values
	switch values := values.(type) {

	// If complex
	case []complex128:
		// Checks that the number of values is with the possible range
		if len(values) > ecd.params.MaxSlots() || len(values) > slots {
			panic(fmt.Sprintf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxSlots (%d)\n", len(values), slots, ecd.params.MaxSlots()))
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
			panic("cannot Embed: ringType must be ring.Standard or ring.ConjugateInvariant")
		}

	// If floats only
	case []float64:
		if len(values) > ecd.params.MaxSlots() || len(values) > slots {
			panic(fmt.Sprintf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxSlots (%d)\n", len(values), slots, ecd.params.MaxSlots()))
		}

		lenValues = len(values)

		for i := range values {
			ecd.values[i] = complex(values[i], 0)
		}

	default:
		panic("cannot Embed: values.(Type) must be []complex128 or []float64")
	}

	for i := lenValues; i < slots; i++ {
		ecd.values[i] = 0
	}

	ecd.IFFT(ecd.values, logSlots)

	switch p := polyOut.(type) {
	case ringqp.Poly:
		ComplexToFixedPointCRT(ecd.params.RingQ().AtLevel(p.Q.Level()), ecd.values[:slots], scale.Float64(), p.Q.Coeffs)
		NttSparseAndMontgomery(ecd.params.RingQ().AtLevel(p.Q.Level()), logSlots, montgomery, p.Q)

		if p.P != nil {
			ComplexToFixedPointCRT(ecd.params.RingP().AtLevel(p.P.Level()), ecd.values[:slots], scale.Float64(), p.P.Coeffs)
			NttSparseAndMontgomery(ecd.params.RingP().AtLevel(p.P.Level()), logSlots, montgomery, p.P)
		}
	case *ring.Poly:
		ComplexToFixedPointCRT(ecd.params.RingQ().AtLevel(p.Level()), ecd.values[:slots], scale.Float64(), p.Coeffs)
		NttSparseAndMontgomery(ecd.params.RingQ().AtLevel(p.Level()), logSlots, montgomery, p)
	default:
		panic("cannot Embed: invalid polyOut.(Type) must be ringqp.Poly or *ring.Poly")
	}
}

func polyToComplexNoCRT(coeffs []uint64, values []complex128, scale rlwe.Scale, logSlots int, isreal bool, ringQ *ring.Ring) {

	slots := 1 << logSlots
	maxSlots := int(ringQ.NthRoot() >> 2)
	gap := maxSlots / slots
	Q := ringQ.SubRings[0].Modulus
	var c uint64
	for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
		c = coeffs[idx]
		if c >= Q>>1 {
			values[i] = complex(-float64(Q-c), 0)
		} else {
			values[i] = complex(float64(c), 0)
		}
	}

	if !isreal {
		for i, idx := 0, maxSlots; i < slots; i, idx = i+1, idx+gap {
			c = coeffs[idx]
			if c >= Q>>1 {
				values[i] += complex(0, -float64(Q-c))
			} else {
				values[i] += complex(0, float64(c))
			}
		}
	}

	divideComplex128SliceVec(values, complex(scale.Float64(), 0))
}

func polyToComplexCRT(poly *ring.Poly, bigintCoeffs []*big.Int, values []complex128, scale rlwe.Scale, logSlots int, isreal bool, ringQ *ring.Ring, Q *big.Int) {

	maxSlots := int(ringQ.NthRoot() >> 2)
	slots := 1 << logSlots
	gap := maxSlots / slots

	ringQ.PolyToBigint(poly, gap, bigintCoeffs)

	qHalf := new(big.Int)
	qHalf.Set(Q)
	qHalf.Rsh(qHalf, 1)

	var sign int

	scalef64 := scale.Float64()

	var c *big.Int
	for i := 0; i < slots; i++ {
		c = bigintCoeffs[i]
		c.Mod(c, Q)
		if sign = c.Cmp(qHalf); sign == 1 || sign == 0 {
			c.Sub(c, Q)
		}
		values[i] = complex(scaleDown(c, scalef64), 0)
	}

	if !isreal {
		for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
			c = bigintCoeffs[j]
			c.Mod(c, Q)
			if sign = c.Cmp(qHalf); sign == 1 || sign == 0 {
				c.Sub(c, Q)
			}
			values[i] += complex(0, scaleDown(c, scalef64))
		}
	}
}

func (ecd *encoderComplex128) plaintextToComplex(level int, scale rlwe.Scale, logSlots int, p *ring.Poly, values []complex128) {

	isreal := ecd.params.RingType() == ring.ConjugateInvariant
	if level == 0 {
		polyToComplexNoCRT(p.Coeffs[0], values, scale, logSlots, isreal, ecd.params.RingQ())
	} else {
		polyToComplexCRT(p, ecd.bigintCoeffs, values, scale, logSlots, isreal, ecd.params.RingQ(), ecd.params.RingQ().ModulusAtLevel[level])
	}

	if isreal { // [X]/(X^N+1) to [X+X^-1]/(X^N+1)
		tmp := ecd.values
		slots := 1 << logSlots
		for i := 1; i < slots; i++ {
			tmp[i] -= complex(0, real(tmp[slots-i]))
		}
	}
}

func (ecd *encoderComplex128) decodePublic(plaintext *rlwe.Plaintext, logSlots int, sigma float64) (res []complex128) {

	if logSlots > ecd.params.MaxLogSlots() || logSlots < minLogSlots {
		panic(fmt.Sprintf("cannot Decode: ensure that %d <= logSlots (%d) <= %d", minLogSlots, logSlots, ecd.params.MaxLogSlots()))
	}

	if plaintext.IsNTT {
		ecd.params.RingQ().AtLevel(plaintext.Level()).INTT(plaintext.Value, ecd.buff)
	} else {
		ring.CopyLvl(plaintext.Level(), plaintext.Value, ecd.buff)
	}

	// B = floor(sigma * sqrt(2*pi))
	if sigma != 0 {
		ecd.gaussianSampler.AtLevel(plaintext.Level()).ReadAndAddFromDist(ecd.buff, ecd.params.RingQ(), sigma, int(2.5066282746310002*sigma))
	}

	ecd.plaintextToComplex(plaintext.Level(), plaintext.Scale, logSlots, ecd.buff, ecd.values)

	ecd.FFT(ecd.values, logSlots)

	res = make([]complex128, 1<<logSlots)
	copy(res, ecd.values)

	return
}

func (ecd *encoderComplex128) decodeCoeffsPublic(plaintext *rlwe.Plaintext, sigma float64) (res []float64) {

	if plaintext.IsNTT {
		ecd.params.RingQ().AtLevel(plaintext.Level()).INTT(plaintext.Value, ecd.buff)
	} else {
		ring.CopyLvl(plaintext.Level(), plaintext.Value, ecd.buff)
	}

	if sigma != 0 {
		// B = floor(sigma * sqrt(2*pi))
		ecd.gaussianSampler.AtLevel(plaintext.Level()).ReadAndAddFromDist(ecd.buff, ecd.params.RingQ(), sigma, int(2.5066282746310002*sigma))
	}

	res = make([]float64, ecd.params.N())

	sf64 := plaintext.Scale.Float64()

	// We have more than one moduli and need the CRT reconstruction
	if plaintext.Level() > 0 {

		ecd.params.RingQ().PolyToBigint(ecd.buff, 1, ecd.bigintCoeffs)

		Q := ecd.params.RingQ().ModulusAtLevel[plaintext.Level()]

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

			res[i] = scaleDown(ecd.bigintCoeffs[i], sf64)
		}
		// We can directly get the coefficients
	} else {

		Q := ecd.params.RingQ().SubRings[0].Modulus
		coeffs := ecd.buff.Coeffs[0]

		for i := range res {

			if coeffs[i] >= Q>>1 {
				res[i] = -float64(Q - coeffs[i])
			} else {
				res[i] = float64(coeffs[i])
			}

			res[i] /= sf64
		}
	}

	return
}

func (ecd *encoderComplex128) IFFT(values []complex128, logN int) {
	if logN < 3 {
		SpecialiFFTVec(values, 1<<logN, ecd.m, ecd.rotGroup, ecd.roots)
	} else {
		SpecialiFFTUL8Vec(values, 1<<logN, ecd.m, ecd.rotGroup, ecd.roots)
	}
}

func (ecd *encoderComplex128) FFT(values []complex128, logN int) {
	if logN < 3 {
		SpecialFFTVec(values, 1<<logN, ecd.m, ecd.rotGroup, ecd.roots)
	} else {
		SpecialFFTUL8Vec(values, 1<<logN, ecd.m, ecd.rotGroup, ecd.roots)
	}
}

// EncoderBigComplex is an interface that implements the encoding algorithms with arbitrary precision.
type EncoderBigComplex interface {
	Encode(values []*ring.Complex, plaintext *rlwe.Plaintext, logSlots int)
	EncodeNew(values []*ring.Complex, level int, scale rlwe.Scale, logSlots int) (plaintext *rlwe.Plaintext)
	Decode(plaintext *rlwe.Plaintext, logSlots int) (res []*ring.Complex)
	DecodePublic(plaintext *rlwe.Plaintext, logSlots int, sigma float64) (res []*ring.Complex)
	FFT(values []*ring.Complex, N int)
	InvFFT(values []*ring.Complex, N int)
	ShallowCopy() EncoderBigComplex
}

type encoderBigComplex struct {
	encoder
	zero        *big.Float
	cMul        *ring.ComplexMultiplier
	prec        uint
	values      []*ring.Complex
	valuesfloat []*big.Float
	roots       []*ring.Complex
}

// NewEncoderBigComplex creates a new encoder using arbitrary precision complex arithmetic.
func NewEncoderBigComplex(params Parameters, prec uint) EncoderBigComplex {

	ecd := newEncoder(params)

	values := make([]*ring.Complex, ecd.m>>2)
	valuesfloat := make([]*big.Float, ecd.m>>1)

	for i := 0; i < ecd.m>>2; i++ {

		values[i] = ring.NewComplex(ring.NewFloat(0, prec), ring.NewFloat(0, prec))
		valuesfloat[i*2] = ring.NewFloat(0, prec)
		valuesfloat[(i*2)+1] = ring.NewFloat(0, prec)
	}

	return &encoderBigComplex{
		encoder:     ecd,
		zero:        ring.NewFloat(0, prec),
		cMul:        ring.NewComplexMultiplier(),
		prec:        prec,
		roots:       GetRootsbigFloat(ecd.m, prec),
		values:      values,
		valuesfloat: valuesfloat,
	}
}

// Encode encodes a set of values on the target plaintext.
// Encoding is done at the level and scale of the plaintext.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^LogN.
func (ecd *encoderBigComplex) Encode(values []*ring.Complex, plaintext *rlwe.Plaintext, logSlots int) {

	slots := 1 << logSlots
	N := ecd.params.N()

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

	gap := (ecd.params.RingQ().N() >> 1) / slots

	for i, jdx, idx := 0, N>>1, 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		ecd.valuesfloat[idx].Set(ecd.values[i].Real())
		ecd.valuesfloat[jdx].Set(ecd.values[i].Imag())
	}

	scaleUpVecExactBigFloat(ecd.valuesfloat, plaintext.Scale.Float64(), ecd.params.RingQ().ModuliChain()[:plaintext.Level()+1], plaintext.Value.Coeffs)

	halfN := N >> 1
	for i := 0; i < halfN; i++ {
		ecd.values[i].Real().Set(ecd.zero)
		ecd.values[i].Imag().Set(ecd.zero)
	}

	for i := 0; i < N; i++ {
		ecd.valuesfloat[i].Set(ecd.zero)
	}

	ecd.params.RingQ().AtLevel(plaintext.Level()).NTT(plaintext.Value, plaintext.Value)
}

// EncodeNew encodes a set of values on a new plaintext.
// Encoding is done at the provided level and with the provided scale.
// User must ensure that 1 <= len(values) <= 2^logSlots < 2^LogN.
func (ecd *encoderBigComplex) EncodeNew(values []*ring.Complex, level int, scale rlwe.Scale, logSlots int) (plaintext *rlwe.Plaintext) {
	plaintext = NewPlaintext(ecd.params, level)
	plaintext.Scale = scale
	ecd.Encode(values, plaintext, logSlots)
	return
}

// Decode decodes the input plaintext on a new slice of ring.Complex.
func (ecd *encoderBigComplex) Decode(plaintext *rlwe.Plaintext, logSlots int) (res []*ring.Complex) {
	return ecd.decodePublic(plaintext, logSlots, 0)
}

func (ecd *encoderBigComplex) DecodePublic(plaintext *rlwe.Plaintext, logSlots int, sigma float64) (res []*ring.Complex) {
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

	NBig := ring.NewFloat(float64(N), ecd.prec)
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

		values[i] = ring.NewComplex(ring.NewFloat(0, ecd.prec), ring.NewFloat(0, ecd.prec))
		valuesfloat[i*2] = ring.NewFloat(0, ecd.prec)
		valuesfloat[(i*2)+1] = ring.NewFloat(0, ecd.prec)
	}

	return &encoderBigComplex{
		encoder:     *ecd.encoder.ShallowCopy(),
		zero:        ring.NewFloat(0, ecd.prec),
		cMul:        ring.NewComplexMultiplier(),
		prec:        ecd.prec,
		values:      values,
		valuesfloat: valuesfloat,
		roots:       ecd.roots,
	}
}

func (ecd *encoderBigComplex) decodePublic(plaintext *rlwe.Plaintext, logSlots int, sigma float64) (res []*ring.Complex) {

	slots := 1 << logSlots

	if logSlots > ecd.params.LogN()-1 {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	ecd.params.RingQ().AtLevel(plaintext.Level()).INTT(plaintext.Value, ecd.buff)

	if sigma != 0 {
		// B = floor(sigma * sqrt(2*pi))
		ecd.gaussianSampler.AtLevel(plaintext.Level()).ReadAndAddFromDist(ecd.buff, ecd.params.RingQ(), sigma, int(2.5066282746310002*sigma+0.5))
	}

	Q := ecd.params.RingQ().ModulusAtLevel[plaintext.Level()]

	maxSlots := ecd.params.N() >> 1

	scaleFlo := plaintext.Scale.Value

	ecd.qHalf.Set(Q)
	ecd.qHalf.Rsh(ecd.qHalf, 1)

	gap := maxSlots / slots

	ecd.params.RingQ().PolyToBigint(ecd.buff, gap, ecd.bigintCoeffs)

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
		ecd.values[i].Real().Quo(ecd.values[i].Real(), &scaleFlo)

		ecd.values[i].Imag().SetInt(ecd.bigintCoeffs[j])
		ecd.values[i].Imag().Quo(ecd.values[i].Imag(), &scaleFlo)
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
