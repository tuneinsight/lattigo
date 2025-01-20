package ckks

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

type Float interface {
	float64 | complex128 | *big.Float | *bignum.Complex
}

// FloatSlice is an empty interface whose goal is to
// indicate that the expected input should be []Float.
// See Float for information on the type constraint.
type FloatSlice interface {
}

// Complex is an empty interface whose goal is to
// indicate that the expected input should be a complex number.
type Complex any

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = ring.GaloisGen

// Encoder is a type that implements the encoding and decoding interface for the CKKS scheme. It provides methods to encode/decode
// []complex128/[]*[bignum.Complex] and []float64/[]*[big.Float] types into/from Plaintext types.
//
// Two different encodings domains are provided:
//
//   - Coefficients: The coefficients are directly embedded on the plaintext. This encoding only allows to encode []float64/[]*[big.Float] slices,
//     but of size up to N (N being the ring degree) and does not preserve the point-wise multiplication. A ciphertext multiplication will result
//     in a negacyclic polynomial convolution in the plaintext domain. This encoding does not provide native slot cyclic rotation.
//     Other operations, like addition or constant multiplication, behave as usual.
//
//   - Slots: The coefficients are first subjected to a special Fourier transform before being embedded in the plaintext by using Coeffs encoding.
//     This encoding can embed []complex128/[]*[bignum.Complex] and []float64/[]*[big.Float] slices of size at most N/2 (N being the ring degree) and
//     leverages the convolution property of the DFT to preserve point-wise complex multiplication in the plaintext domain, i.e. a ciphertext
//     multiplication will result in an element-wise multiplication in the plaintext domain. It also enables cyclic rotations on plaintext slots.
//     Other operations, like constant multiplication, behave as usual. It is considered the default encoding method for CKKS.
//
// The figure bellow illustrates the relationship between these two encodings:
//
//	                                                    Z_Q[X]/(X^N+1)
//		Coefficients: ---------------> Real^{N} ---------> Plaintext
//	                                      |
//	                                      |
//		Slots: Complex^{N/2} -> iDFT -----┘
type Encoder struct {
	parameters Parameters

	prec uint

	// bigintCoeffs []*big.Int
	// qHalf        *big.Int
	// buff     ring.Poly
	m        int
	rotGroup []int

	roots interface{}

	// Pools used to recycle large objects.
	BuffPolyPool    structs.BufferPool[*ring.Poly]
	BuffBigIntPool  structs.BufferPool[*[]*big.Int]
	BuffComplexPool structs.BufferPool[Complex]
}

// NewEncoder creates a new [Encoder] from the target parameters.
// Optional field `precision` can be given. If precision is empty
// or <= 53, then float64 and complex128 types will be used to
// perform the encoding. Else *[big.Float] and *[bignum.Complex] will be used.
func NewEncoder(parameters Parameters, precision ...uint) (ecd *Encoder) {

	m := int(parameters.RingQ().NthRoot())

	rotGroup := make([]int, m>>2)
	fivePows := 1
	for i := 0; i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= int(GaloisGen)
		fivePows &= (m - 1)
	}

	var prec uint
	if len(precision) != 0 && precision[0] != 0 {
		prec = precision[0]
	} else {
		prec = parameters.EncodingPrecision()
	}

	ecd = &Encoder{
		prec:       prec,
		parameters: parameters,
		m:          m,
		rotGroup:   rotGroup,
	}

	ecd.BuffBigIntPool = structs.NewSyncPool(func() *[]*big.Int {
		buff := make([]*big.Int, m>>1)
		return &buff
	})

	ringQ := parameters.RingQ()
	ecd.BuffPolyPool = ringQ.NewBuffFromUintPool()

	if prec <= 53 {

		ecd.roots = GetRootsComplex128(ecd.m)

		ecd.BuffComplexPool = structs.NewSyncPool(func() Complex {
			buff := make([]complex128, ecd.m>>2)
			return &buff
		})

	} else {

		ecd.roots = GetRootsBigComplex(ecd.m, prec)

		ecd.BuffComplexPool = structs.NewSyncPool(func() Complex {
			buff := make([]*bignum.Complex, ecd.m>>2)
			for i := 0; i < ecd.m>>2; i++ {
				buff[i] = &bignum.Complex{bignum.NewFloat(0, prec), bignum.NewFloat(0, prec)}
			}
			return &buff
		})
	}

	return
}

// Prec returns the precision in bits used by the target Encoder.
// A precision <= 53 will use float64, else *[big.Float].
func (ecd Encoder) Prec() uint {
	return ecd.prec
}

func (ecd Encoder) GetParameters() Parameters {
	return ecd.parameters
}

func (ecd Encoder) GetRLWEParameters() rlwe.Parameters {
	return ecd.parameters.Parameters
}

// Encode encodes a [FloatSlice] on the target plaintext.
// Encoding is done at the level and scale of the plaintext.
// Encoding domain is done according to the metadata of the plaintext.
// User must ensure that 1 <= len(values) <= 2^pt.LogMaxDimensions < 2^logN.
// The imaginary part will be discarded if ringType == ring.ConjugateInvariant.
func (ecd Encoder) Encode(values interface{}, pt *rlwe.Plaintext) (err error) {
	if pt.IsBatched {
		return ecd.Embed(values, pt.MetaData, pt.Value)

	} else {

		switch values := values.(type) {
		case []float64:

			if len(values) > ecd.parameters.N() {
				return fmt.Errorf("cannot Encode: maximum number of values is %d but len(values) is %d", ecd.parameters.N(), len(values))
			}

			Float64ToFixedPointCRT(ecd.parameters.RingQ().AtLevel(pt.Level()), values, pt.Scale.Float64(), pt.Value.Coeffs)

		case []*big.Float:

			if len(values) > ecd.parameters.N() {
				return fmt.Errorf("cannot Encode: maximum number of values is %d but len(values) is %d", ecd.parameters.N(), len(values))
			}

			BigFloatToFixedPointCRT(ecd.parameters.RingQ().AtLevel(pt.Level()), values, &pt.Scale.Value, pt.Value.Coeffs)

		default:
			return fmt.Errorf("cannot Encode: supported values.(type) for IsBatched=False is []float64 or []*big.Float, but %T was given", values)
		}

		ecd.parameters.RingQ().AtLevel(pt.Level()).NTT(pt.Value, pt.Value)
	}

	return
}

// Decode decodes the input plaintext on a new FloatSlice.
func (ecd Encoder) Decode(pt *rlwe.Plaintext, values interface{}) (err error) {
	return ecd.DecodePublic(pt, values, 0)
}

// DecodePublic decodes the input plaintext on a [FloatSlice].
// It adds, before the decoding step (i.e. in the Ring) noise that follows the given distribution parameters.
// If the underlying ringType is [ring.ConjugateInvariant], the imaginary part (and its related error) are zero.
func (ecd Encoder) DecodePublic(pt *rlwe.Plaintext, values FloatSlice, logprec float64) (err error) {
	return ecd.decodePublic(pt, values, logprec)
}

// Embed is a generic method to encode a [FloatSlice] on the target polyOut.
// This method it as the core of the slot encoding.
// Values are encoded according to the provided metadata.
// Accepted polyOut.(type) are ringqp.Poly and ring.Poly.
// The imaginary part will be discarded if ringType == ring.ConjugateInvariant.
func (ecd Encoder) Embed(values interface{}, metadata *rlwe.MetaData, polyOut interface{}) (err error) {
	if ecd.prec <= 53 {
		return ecd.embedDouble(values, metadata, polyOut)
	}

	return ecd.embedArbitrary(values, metadata, polyOut)
}

// embedDouble encode a FloatSlice on polyOut using FFT with complex128 arithmetic.
// Values are encoded according to the provided metadata.
// Accepted polyOut.(type) are [ringqp.Poly] and [ring.Poly].
func (ecd Encoder) embedDouble(values FloatSlice, metadata *rlwe.MetaData, polyOut interface{}) (err error) {

	if maxLogCols := ecd.parameters.LogMaxDimensions().Cols; metadata.LogDimensions.Cols < 0 || metadata.LogDimensions.Cols > maxLogCols {
		return fmt.Errorf("cannot Embed: logSlots (%d) must be greater or equal to %d and smaller than %d", metadata.LogDimensions.Cols, 0, maxLogCols)
	}

	slots := 1 << metadata.LogDimensions.Cols
	var lenValues int

	buffRef := ecd.BuffComplexPool.Get().(*[]complex128)
	defer ecd.BuffComplexPool.Put(buffRef)
	buffCmplx := *buffRef

	switch values := values.(type) {

	case []complex128:

		lenValues = len(values)

		if maxCols := ecd.parameters.MaxDimensions().Cols; lenValues > maxCols || lenValues > slots {
			return fmt.Errorf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxCols (%d)", len(values), slots, maxCols)
		}

		if ecd.parameters.RingType() == ring.ConjugateInvariant {
			for i := range values {
				buffCmplx[i] = complex(real(values[i]), 0)
			}
		} else {
			copy(buffCmplx[:len(values)], values)
		}

	case []*bignum.Complex:

		lenValues = len(values)

		if maxCols := ecd.parameters.MaxDimensions().Cols; lenValues > maxCols || lenValues > slots {
			return fmt.Errorf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxCols (%d)", len(values), slots, maxCols)
		}

		if ecd.parameters.RingType() == ring.ConjugateInvariant {
			for i := range values {
				if values[i] != nil {
					f64, _ := values[i][0].Float64()
					buffCmplx[i] = complex(f64, 0)
				} else {
					buffCmplx[i] = 0
				}
			}
		} else {
			for i := range values {
				if values[i] != nil {
					buffCmplx[i] = values[i].Complex128()
				} else {
					buffCmplx[i] = 0
				}
			}
		}

	case []float64:

		lenValues = len(values)

		if maxCols := ecd.parameters.MaxDimensions().Cols; lenValues > maxCols || lenValues > slots {
			return fmt.Errorf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxCols (%d)", len(values), slots, maxCols)
		}

		for i := range values {
			buffCmplx[i] = complex(values[i], 0)
		}

	case []*big.Float:

		lenValues = len(values)

		if maxCols := ecd.parameters.MaxDimensions().Cols; lenValues > maxCols || lenValues > slots {
			return fmt.Errorf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxCols (%d)", len(values), slots, maxCols)
		}

		for i := range values {
			if values[i] != nil {
				f64, _ := values[i].Float64()
				buffCmplx[i] = complex(f64, 0)
			} else {
				buffCmplx[i] = 0
			}
		}
	default:
		return fmt.Errorf("cannot Embed: values.(Type) must be []complex128, []*bignum.Complex, []float64 or []*big.Float, but is %T", values)
	}

	// Zeroes all other values
	for i := lenValues; i < slots; i++ {
		buffCmplx[i] = 0
	}

	// IFFT
	if err = ecd.IFFT(buffCmplx[:slots], metadata.LogDimensions.Cols); err != nil {
		return
	}

	// Maps Y = X^{N/n} -> X and quantizes.
	switch p := polyOut.(type) {
	case ringqp.Poly:
		Complex128ToFixedPointCRT(ecd.parameters.RingQ().AtLevel(p.Q.Level()), buffCmplx[:slots], metadata.Scale.Float64(), p.Q.Coeffs)
		rlwe.NTTSparseAndMontgomery(ecd.parameters.RingQ().AtLevel(p.Q.Level()), metadata, p.Q)

		if p.P.Level() > -1 {
			Complex128ToFixedPointCRT(ecd.parameters.RingP().AtLevel(p.P.Level()), buffCmplx[:slots], metadata.Scale.Float64(), p.P.Coeffs)
			rlwe.NTTSparseAndMontgomery(ecd.parameters.RingP().AtLevel(p.P.Level()), metadata, p.P)
		}
	case ring.Poly:
		Complex128ToFixedPointCRT(ecd.parameters.RingQ().AtLevel(p.Level()), buffCmplx[:slots], metadata.Scale.Float64(), p.Coeffs)
		rlwe.NTTSparseAndMontgomery(ecd.parameters.RingQ().AtLevel(p.Level()), metadata, p)
	default:
		return fmt.Errorf("cannot Embed: invalid polyOut.(Type) must be ringqp.Poly or ring.Poly")
	}

	return
}

// embedArbitrary encode a FloatSlice on polyOut using FFT with *bignum.Complex arithmetic.
// Values are encoded according to the provided metadata.
// Accepted polyOut.(type) are [ringqp.Poly] and [ring.Poly].
func (ecd Encoder) embedArbitrary(values FloatSlice, metadata *rlwe.MetaData, polyOut interface{}) (err error) {

	if maxLogCols := ecd.parameters.LogMaxDimensions().Cols; metadata.LogDimensions.Cols < 0 || metadata.LogDimensions.Cols > maxLogCols {
		return fmt.Errorf("cannot Embed: logSlots (%d) must be greater or equal to %d and smaller than %d", metadata.LogDimensions.Cols, 0, maxLogCols)
	}

	slots := 1 << metadata.LogDimensions.Cols
	var lenValues int

	buffRef := ecd.BuffComplexPool.Get().(*[]*bignum.Complex)
	defer ecd.BuffComplexPool.Put(buffRef)
	buffCmplx := *buffRef

	switch values := values.(type) {

	case []complex128:

		lenValues = len(values)

		if maxCols := ecd.parameters.MaxDimensions().Cols; lenValues > maxCols || lenValues > slots {
			return fmt.Errorf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxCols (%d)", len(values), slots, maxCols)
		}

		if ecd.parameters.RingType() == ring.ConjugateInvariant {
			for i := range values {
				buffCmplx[i][0].SetFloat64(real(values[i]))
				buffCmplx[i][1].SetFloat64(0)
			}
		} else {
			for i := range values {
				buffCmplx[i][0].SetFloat64(real(values[i]))
				buffCmplx[i][1].SetFloat64(imag(values[i]))
			}
		}

	case []*bignum.Complex:

		lenValues = len(values)

		if maxCols := ecd.parameters.MaxDimensions().Cols; lenValues > maxCols || lenValues > slots {
			return fmt.Errorf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxCols (%d)", len(values), slots, maxCols)
		}

		if ecd.parameters.RingType() == ring.ConjugateInvariant {
			for i := range values {
				if values[i] != nil {
					buffCmplx[i][0].Set(values[i][0])
				} else {
					buffCmplx[i][0].SetFloat64(0)
				}

				buffCmplx[i][1].SetFloat64(0)
			}
		} else {
			for i := range values {
				if values[i] != nil {
					buffCmplx[i].Set(values[i])
				} else {
					buffCmplx[i][0].SetFloat64(0)
					buffCmplx[i][1].SetFloat64(0)
				}
			}
		}

	case []float64:

		lenValues = len(values)

		if maxCols := ecd.parameters.MaxDimensions().Cols; lenValues > maxCols || lenValues > slots {
			return fmt.Errorf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxCols (%d)", len(values), slots, maxCols)
		}

		for i := range values {
			buffCmplx[i][0].SetFloat64(values[i])
			buffCmplx[i][1].SetFloat64(0)
		}

	case []*big.Float:

		lenValues = len(values)

		if maxCols := ecd.parameters.MaxDimensions().Cols; lenValues > maxCols || lenValues > slots {
			return fmt.Errorf("cannot Embed: ensure that #values (%d) <= slots (%d) <= maxCols (%d)", len(values), slots, maxCols)
		}

		for i := range values {
			if values[i] != nil {
				buffCmplx[i][0].Set(values[i])
			} else {
				buffCmplx[i][0].SetFloat64(0)
			}

			buffCmplx[i][1].SetFloat64(0)
		}
	default:
		return fmt.Errorf("cannot Embed: values.(Type) must be []complex128, []*bignum.Complex, []float64 or []*big.Float, but is %T", values)
	}

	// Zeroes all other values
	for i := lenValues; i < slots; i++ {
		buffCmplx[i][0].SetFloat64(0)
		buffCmplx[i][1].SetFloat64(0)
	}

	if err = ecd.IFFT(buffCmplx[:slots], metadata.LogDimensions.Cols); err != nil {
		return
	}

	// Maps Y = X^{N/n} -> X and quantizes.
	switch p := polyOut.(type) {

	case ring.Poly:

		ComplexArbitraryToFixedPointCRT(ecd.parameters.RingQ().AtLevel(p.Level()), buffCmplx[:slots], &metadata.Scale.Value, p.Coeffs)
		rlwe.NTTSparseAndMontgomery(ecd.parameters.RingQ().AtLevel(p.Level()), metadata, p)

	case ringqp.Poly:

		ComplexArbitraryToFixedPointCRT(ecd.parameters.RingQ().AtLevel(p.Q.Level()), buffCmplx[:slots], &metadata.Scale.Value, p.Q.Coeffs)
		rlwe.NTTSparseAndMontgomery(ecd.parameters.RingQ().AtLevel(p.Q.Level()), metadata, p.Q)

		if p.P.Level() > -1 {
			ComplexArbitraryToFixedPointCRT(ecd.parameters.RingP().AtLevel(p.P.Level()), buffCmplx[:slots], &metadata.Scale.Value, p.P.Coeffs)
			rlwe.NTTSparseAndMontgomery(ecd.parameters.RingP().AtLevel(p.P.Level()), metadata, p.P)
		}

	default:
		return fmt.Errorf("cannot Embed: invalid polyOut.(Type) must be ringqp.Poly or ring.Poly")
	}

	return
}

// plaintextToComplex maps a CRT polynomial to a complex valued [FloatSlice].
func (ecd Encoder) plaintextToComplex(level int, scale rlwe.Scale, logSlots int, p ring.Poly, values FloatSlice) (err error) {

	isreal := ecd.parameters.RingType() == ring.ConjugateInvariant
	if level == 0 {
		return polyToComplexNoCRT(p.Coeffs[0], values, scale, logSlots, isreal, ecd.parameters.RingQ().AtLevel(level))
	}
	bigintCoeffs := ecd.BuffBigIntPool.Get()
	defer ecd.BuffBigIntPool.Put(bigintCoeffs)
	return polyToComplexCRT(p, *bigintCoeffs, values, scale, logSlots, isreal, ecd.parameters.RingQ().AtLevel(level))
}

// plaintextToFloat maps a CRT polynomial to a real valued [FloatSlice].
func (ecd Encoder) plaintextToFloat(level int, scale rlwe.Scale, logSlots int, p ring.Poly, values FloatSlice) (err error) {
	if level == 0 {
		return ecd.polyToFloatNoCRT(p.Coeffs[0], values, scale, logSlots, ecd.parameters.RingQ().AtLevel(level))
	}
	return ecd.polyToFloatCRT(p, values, scale, logSlots, ecd.parameters.RingQ().AtLevel(level))
}

// decodePublic decode a plaintext to a [FloatSlice].
// The method will add a flooding noise before the decoding process following the defined distribution if it is not nil.
func (ecd Encoder) decodePublic(pt *rlwe.Plaintext, values FloatSlice, logprec float64) (err error) {

	logSlots := pt.LogDimensions.Cols
	slots := 1 << logSlots

	if maxLogCols := ecd.parameters.LogMaxDimensions().Cols; logSlots > maxLogCols || logSlots < 0 {
		return fmt.Errorf("cannot Decode: ensure that %d <= logSlots (%d) <= %d", 0, logSlots, maxLogCols)
	}

	buff := ecd.BuffPolyPool.Get()
	defer ecd.BuffPolyPool.Put(buff)
	if pt.IsNTT {
		ecd.parameters.RingQ().AtLevel(pt.Level()).INTT(pt.Value, *buff)
	} else {
		buff.CopyLvl(pt.Level(), pt.Value)
	}

	switch values.(type) {
	case []complex128, []float64, []*bignum.Complex, []*big.Float:
	default:
		return fmt.Errorf("cannot decode: values.(type) accepted are []complex128, []float64, []*bignum.Complex, []*big.Float but is %T", values)
	}

	if pt.IsBatched {

		if ecd.prec <= 53 {

			buffRef := ecd.BuffComplexPool.Get().(*[]complex128)
			defer ecd.BuffComplexPool.Put(buffRef)
			buffCmplx := *buffRef

			if err = ecd.plaintextToComplex(pt.Level(), pt.Scale, logSlots, *buff, buffCmplx); err != nil {
				return
			}

			if err = ecd.FFT(buffCmplx[:slots], logSlots); err != nil {
				return
			}

			if logprec != 0 {

				scale := math.Exp2(logprec)

				switch values.(type) {
				case []*bignum.Complex, []complex128:
					for i := 0; i < slots; i++ {
						buffCmplx[i] = complex(math.Round(real(buffCmplx[i])*scale)/scale, math.Round(imag(buffCmplx[i])*scale)/scale)
					}
				default:
					for i := 0; i < slots; i++ {
						buffCmplx[i] = complex(math.Round(real(buffCmplx[i])*scale)/scale, 0)
					}
				}
			}

			switch values := values.(type) {
			case []float64:

				slots := utils.Min(len(values), slots)

				for i := 0; i < slots; i++ {
					values[i] = real(buffCmplx[i])
				}

			case []complex128:
				copy(values, buffCmplx)
			case []*big.Float:

				slots := utils.Min(len(values), slots)

				for i := 0; i < slots; i++ {

					if values[i] == nil {
						values[i] = new(big.Float)
					}

					values[i].SetFloat64(real(buffCmplx[i]))
				}

			case []*bignum.Complex:

				slots := utils.Min(len(values), slots)

				for i := 0; i < slots; i++ {

					if values[i] == nil {
						values[i] = &bignum.Complex{
							new(big.Float),
							new(big.Float),
						}
					} else {
						if values[i][0] == nil {
							values[i][0] = new(big.Float)
						}

						if values[i][1] == nil {
							values[i][1] = new(big.Float)
						}
					}

					values[i][0].SetFloat64(real(buffCmplx[i]))
					values[i][1].SetFloat64(imag(buffCmplx[i]))
				}
			}
		} else {

			buffRef := ecd.BuffComplexPool.Get().(*[]*bignum.Complex)
			defer ecd.BuffComplexPool.Put(buffRef)
			buffCmplx := *buffRef

			if err = ecd.plaintextToComplex(pt.Level(), pt.Scale, logSlots, *buff, buffCmplx[:slots]); err != nil {
				return
			}

			if err = ecd.FFT(buffCmplx[:slots], logSlots); err != nil {
				return
			}

			var scale, half, zero *big.Float
			var tmp *big.Int
			if logprec != 0 {

				// 2^logprec
				scale = new(big.Float).SetPrec(ecd.Prec()).SetFloat64(logprec)
				scale.Mul(scale, bignum.Log2(ecd.Prec()))
				scale = bignum.Exp(scale)

				tmp = new(big.Int)
				half = new(big.Float).SetFloat64(0.5)
				zero = new(big.Float)
			}

			switch values := values.(type) {
			case []float64:

				slots := utils.Min(len(values), slots)

				for i := 0; i < slots; i++ {
					values[i], _ = buffCmplx[i][0].Float64()
				}

				if logprec != 0 {

					scaleF64, _ := scale.Float64()

					for i := 0; i < slots; i++ {
						values[i] = math.Round(values[i]*scaleF64) / scaleF64
					}
				}

			case []complex128:

				slots := utils.Min(len(values), slots)

				for i := 0; i < slots; i++ {
					values[i] = buffCmplx[i].Complex128()
				}

				if logprec != 0 {

					scaleF64, _ := scale.Float64()

					for i := 0; i < slots; i++ {
						values[i] = complex(math.Round(real(values[i])*scaleF64)/scaleF64, math.Round(imag(values[i])*scaleF64)/scaleF64)
					}
				}

			case []*big.Float:
				slots := utils.Min(len(values), slots)

				for i := 0; i < slots; i++ {

					if values[i] == nil {
						values[i] = new(big.Float)
					}

					values[i].Set(buffCmplx[i][0])
				}
				if logprec != 0 {
					for i := range values {
						values[i].Mul(values[i], scale)

						// Adds/Subtracts 0.5
						if values[i].Cmp(zero) >= 0 {
							values[i].Add(values[i], half)
						} else {
							values[i].Sub(values[i], half)
						}

						// Round = floor +/- 0.5
						values[i].Int(tmp)

						values[i].SetInt(tmp)

						values[i].Quo(values[i], scale)
					}
				}

			case []*bignum.Complex:

				slots := utils.Min(len(values), slots)

				for i := 0; i < slots; i++ {

					if values[i] == nil {
						values[i] = &bignum.Complex{
							new(big.Float),
							new(big.Float),
						}
					} else {
						if values[i][0] == nil {
							values[i][0] = new(big.Float)
						}

						if values[i][1] == nil {
							values[i][1] = new(big.Float)
						}
					}

					values[i][0].Set(buffCmplx[i][0])
					values[i][1].Set(buffCmplx[i][1])
				}

				if logprec != 0 {
					for i := range values {

						// Real
						values[i][0].Mul(values[i][0], scale)

						// Adds/Subtracts 0.5
						if values[i][0].Cmp(zero) >= 0 {
							values[i][0].Add(values[i][0], half)
						} else {
							values[i][0].Sub(values[i][0], half)
						}

						// Round = floor +/- 0.5
						values[i][0].Int(tmp)
						values[i][0].SetInt(tmp)
						values[i][0].Quo(values[i][0], scale)

						// Imag
						values[i][1].Mul(values[i][1], scale)

						// Adds/Subtracts 0.5
						if values[i][1].Cmp(zero) >= 0 {
							values[i][1].Add(values[i][1], half)
						} else {
							values[i][1].Sub(values[i][1], half)
						}

						// Round = floor +/- 0.5
						values[i][1].Int(tmp)
						values[i][1].SetInt(tmp)
						values[i][1].Quo(values[i][1], scale)
					}
				}
			}
		}

	} else {
		return ecd.plaintextToFloat(pt.Level(), pt.Scale, logSlots, *buff, values)
	}

	return
}

// IFFT evaluates the special 2^{LogN}-th encoding discrete Fourier transform on [FloatSlice].
func (ecd Encoder) IFFT(values FloatSlice, logN int) (err error) {
	switch values := values.(type) {
	case []complex128:
		switch roots := ecd.roots.(type) {
		case []complex128:
			if logN < 4 {
				SpecialIFFTDouble(values, 1<<logN, ecd.m, ecd.rotGroup, ecd.roots.([]complex128))
			} else {
				SpecialiFFTDoubleUnrolled8(values, 1<<logN, ecd.m, ecd.rotGroup, ecd.roots.([]complex128))
			}
		default:
			return fmt.Errorf("cannot IFFT: values.(type)=%T doesn't roots.(type) = %T", values, roots)
		}

	case []*bignum.Complex:
		switch roots := ecd.roots.(type) {
		case []*bignum.Complex:
			SpecialIFFTArbitrary(values, 1<<logN, ecd.m, ecd.rotGroup, ecd.roots.([]*bignum.Complex))
		default:
			return fmt.Errorf("cannot IFFT: values.(type)=%T doesn't roots.(type) = %T", values, roots)
		}
	default:
		return fmt.Errorf("cannot IFFT: invalid values.(type), accepted types are []complex128 and []*bignum.Complex but is %T", values)
	}
	return

}

// FFT evaluates the special 2^{LogN}-th decoding discrete Fourier transform on [FloatSlice].
func (ecd Encoder) FFT(values FloatSlice, logN int) (err error) {
	switch values := values.(type) {
	case []complex128:
		switch roots := ecd.roots.(type) {
		case []complex128:
			if logN < 4 {
				SpecialFFTDouble(values, 1<<logN, ecd.m, ecd.rotGroup, roots)
			} else {
				SpecialFFTDoubleUL8(values, 1<<logN, ecd.m, ecd.rotGroup, roots)
			}
		default:
			return fmt.Errorf("cannot IFFT: values.(type)=%T doesn't roots.(type) = %T", values, roots)
		}

	case []*bignum.Complex:
		switch roots := ecd.roots.(type) {
		case []*bignum.Complex:
			SpecialFFTArbitrary(values, 1<<logN, ecd.m, ecd.rotGroup, roots)
		default:
			return fmt.Errorf("cannot IFFT: values.(type)=%T doesn't roots.(type) = %T", values, roots)
		}

	default:
		return fmt.Errorf("cannot IFFT: invalid values.(type), accepted types are []complex128 and []*bignum.Complex but is %T", values)
	}
	return
}

// polyToComplexNoCRT decodes a single-level CRT poly on a complex valued [FloatSlice].
func polyToComplexNoCRT(coeffs []uint64, values FloatSlice, scale rlwe.Scale, logSlots int, isreal bool, ringQ *ring.Ring) (err error) {

	slots := 1 << logSlots
	maxCols := int(ringQ.NthRoot() >> 2)
	gap := maxCols / slots
	Q := ringQ.SubRings[0].Modulus
	var c uint64

	switch values := values.(type) {
	case []complex128:
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			c = coeffs[idx]
			if c >= Q>>1 {
				values[i] = complex(-float64(Q-c), 0)
			} else {
				values[i] = complex(float64(c), 0)
			}
		}

		if !isreal {
			for i, idx := 0, maxCols; i < slots; i, idx = i+1, idx+gap {
				c = coeffs[idx]
				if c >= Q>>1 {
					values[i] += complex(0, -float64(Q-c))
				} else {
					values[i] += complex(0, float64(c))
				}
			}
		} else {
			// [X]/(X^N+1) to [X+X^-1]/(X^N+1)
			slots := 1 << logSlots
			for i := 1; i < slots; i++ {
				values[i] -= complex(0, real(values[slots-i]))
			}
		}

		divideComplex128SliceUnrolled8(values, complex(scale.Float64(), 0))

	case []*bignum.Complex:

		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {

			if values[i] == nil {
				values[i] = &bignum.Complex{
					new(big.Float),
					nil,
				}
			} else {
				if values[i][0] == nil {
					values[i][0] = new(big.Float)
				}
			}

			if c = coeffs[idx]; c >= Q>>1 {
				values[i][0].SetInt64(-int64(Q - c))
			} else {
				values[i][0].SetInt64(int64(c))
			}
		}

		if !isreal {
			for i, idx := 0, maxCols; i < slots; i, idx = i+1, idx+gap {

				if values[i][1] == nil {
					values[i][1] = new(big.Float)
				}

				if c = coeffs[idx]; c >= Q>>1 {
					values[i][1].SetInt64(-int64(Q - c))
				} else {
					values[i][1].SetInt64(int64(c))
				}
			}
		} else {
			slots := 1 << logSlots

			for i := 1; i < slots; i++ {
				values[i][1].Sub(values[i][1], values[slots-i][0])
			}
		}

		s := &scale.Value

		for i := range values {
			values[i][0].Quo(values[i][0], s)
			values[i][1].Quo(values[i][1], s)
		}

	default:
		return fmt.Errorf("cannot polyToComplexNoCRT: values.(Type) must be []complex128 or []*bignum.Complex but is %T", values)
	}

	return
}

// polyToComplexNoCRT decodes a multiple-level CRT poly on a complex valued [FloatSlice].
func polyToComplexCRT(poly ring.Poly, bigintCoeffs []*big.Int, values FloatSlice, scale rlwe.Scale, logSlots int, isreal bool, ringQ *ring.Ring) (err error) {

	maxCols := int(ringQ.NthRoot() >> 2)
	slots := 1 << logSlots
	gap := maxCols / slots

	ringQ.PolyToBigint(poly, gap, bigintCoeffs)

	Q := ringQ.ModulusAtLevel[ringQ.Level()]

	qHalf := new(big.Int)
	qHalf.Set(Q)
	qHalf.Rsh(qHalf, 1)

	var sign int

	switch values := values.(type) {

	case []complex128:
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
		} else {
			// [X]/(X^N+1) to [X+X^-1]/(X^N+1)
			slots := 1 << logSlots
			for i := 1; i < slots; i++ {
				values[i] -= complex(0, real(values[slots-i]))
			}
		}
	case []*bignum.Complex:

		var c *big.Int
		for i := 0; i < slots; i++ {
			c = bigintCoeffs[i]
			c.Mod(c, Q)
			if sign = c.Cmp(qHalf); sign == 1 || sign == 0 {
				c.Sub(c, Q)
			}

			if values[i] == nil {
				values[i] = &bignum.Complex{
					new(big.Float),
					nil,
				}
			} else {
				if values[i][0] == nil {
					values[i][0] = new(big.Float)
				}
			}

			values[i][0].SetInt(c)
		}

		if !isreal {
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				c = bigintCoeffs[j]
				c.Mod(c, Q)
				if sign = c.Cmp(qHalf); sign == 1 || sign == 0 {
					c.Sub(c, Q)
				}

				if values[i][1] == nil {
					values[i][1] = new(big.Float)
				}

				values[i][1].SetInt(c)
			}
		} else {
			// [X]/(X^N+1) to [X+X^-1]/(X^N+1)
			slots := 1 << logSlots
			for i := 1; i < slots; i++ {
				values[i][1].Sub(values[i][1], values[slots-i][0])
			}
		}

		s := &scale.Value

		for i := range values {
			values[i][0].Quo(values[i][0], s)
			values[i][1].Quo(values[i][1], s)
		}

	default:
		return fmt.Errorf("cannot polyToComplexNoCRT: values.(Type) must be []complex128 or []*bignum.Complex but is %T", values)
	}

	return
}

// polyToFloatCRT decodes a multiple-level CRT poly on a real valued [FloatSlice].
func (ecd *Encoder) polyToFloatCRT(p ring.Poly, values FloatSlice, scale rlwe.Scale, logSlots int, r *ring.Ring) (err error) {

	var slots int
	switch values := values.(type) {
	case []float64:
		slots = utils.Min(len(p.Coeffs[0]), len(values))
	case []complex128:
		slots = utils.Min(len(p.Coeffs[0]), len(values))
	case []*big.Float:
		slots = utils.Min(len(p.Coeffs[0]), len(values))
	case []*bignum.Complex:
		slots = utils.Min(len(p.Coeffs[0]), len(values))
	}

	buffRef := ecd.BuffBigIntPool.Get()
	defer ecd.BuffBigIntPool.Put(buffRef)
	bigintCoeffs := *buffRef

	ecd.parameters.RingQ().PolyToBigint(p, 1, bigintCoeffs)

	Q := r.ModulusAtLevel[r.Level()]

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
	}

	switch values := values.(type) {

	case []float64:
		sf64 := scale.Float64()
		for i := 0; i < slots; i++ {
			values[i] = scaleDown(bigintCoeffs[i], sf64)
		}
	case []complex128:
		sf64 := scale.Float64()
		for i := 0; i < slots; i++ {
			values[i] = complex(scaleDown(bigintCoeffs[i], sf64), 0)
		}
	case []*big.Float:
		s := &scale.Value
		for i := 0; i < slots; i++ {

			if values[i] == nil {
				values[i] = new(big.Float)
			}

			values[i].SetInt(bigintCoeffs[i])
			values[i].Quo(values[i], s)
		}
	case []*bignum.Complex:
		s := &scale.Value
		for i := 0; i < slots; i++ {

			if values[i] == nil {
				values[i] = &bignum.Complex{
					new(big.Float),
					new(big.Float),
				}
			} else {
				if values[i][0] == nil {
					values[i][0] = new(big.Float)
				}
			}

			values[i][0].SetInt(bigintCoeffs[i])
			values[i][0].Quo(values[i][0], s)
		}
	default:
		return fmt.Errorf("cannot polyToComplexNoCRT: values.(Type) must be []complex128, []*bignum.Complex, []float64 or []*big.Float but is %T", values)
	}

	return
}

// polyToFloatNoCRT decodes a single-level CRT poly on a real valued [FloatSlice].
func (ecd *Encoder) polyToFloatNoCRT(coeffs []uint64, values FloatSlice, scale rlwe.Scale, logSlots int, r *ring.Ring) (err error) {

	Q := r.SubRings[0].Modulus

	var slots int
	switch values := values.(type) {
	case []float64:
		slots = utils.Min(len(coeffs), len(values))
	case []complex128:
		slots = utils.Min(len(coeffs), len(values))
	case []*big.Float:
		slots = utils.Min(len(coeffs), len(values))
	case []*bignum.Complex:
		slots = utils.Min(len(coeffs), len(values))
	}

	switch values := values.(type) {

	case []float64:

		sf64 := scale.Float64()

		for i := 0; i < slots; i++ {
			if coeffs[i] >= Q>>1 {
				values[i] = -float64(Q-coeffs[i]) / sf64
			} else {
				values[i] = float64(coeffs[i]) / sf64
			}
		}

	case []complex128:

		sf64 := scale.Float64()

		for i := 0; i < slots; i++ {
			if coeffs[i] >= Q>>1 {
				values[i] = complex(-float64(Q-coeffs[i])/sf64, 0)
			} else {
				values[i] = complex(float64(coeffs[i])/sf64, 0)
			}
		}

	case []*big.Float:

		s := &scale.Value

		for i := 0; i < slots; i++ {

			if values[i] == nil {
				values[i] = new(big.Float)
			}

			if coeffs[i] >= Q>>1 {
				values[i].SetInt64(-int64(Q - coeffs[i]))
			} else {
				values[i].SetInt64(int64(coeffs[i]))
			}

			values[i].Quo(values[i], s)
		}

	case []*bignum.Complex:

		s := &scale.Value

		for i := 0; i < slots; i++ {

			if values[i] == nil {
				values[i] = &bignum.Complex{
					new(big.Float),
					nil,
				}
			} else {
				if values[i][0] == nil {
					values[i][0] = new(big.Float)
				}
			}

			if coeffs[i] >= Q>>1 {
				values[i][0].SetInt64(-int64(Q - coeffs[i]))
			} else {
				values[i][0].SetInt64(int64(coeffs[i]))
			}

			values[i][0].Quo(values[i][0], s)
		}

	default:
		return fmt.Errorf("cannot polyToComplexNoCRT: values.(Type) must be []complex128, []*bignum.Complex, []float64 or []*big.Float but is %T", values)
	}

	return
}

// ShallowCopy returns a lightweight copy of the target object
// that can be used concurrently with the original object.
func (ecd Encoder) ShallowCopy() *Encoder {

	return &Encoder{
		prec:            ecd.prec,
		parameters:      ecd.parameters,
		m:               ecd.m,
		rotGroup:        ecd.rotGroup,
		roots:           ecd.roots,
		BuffPolyPool:    ecd.BuffPolyPool,
		BuffBigIntPool:  ecd.BuffBigIntPool,
		BuffComplexPool: ecd.BuffComplexPool,
	}
}
