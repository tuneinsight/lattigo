package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/big"
)

// Encoder is an interface implenting the encoding algorithms.
type Encoder interface {
	Encode(plaintext *Plaintext, values []complex128, slots uint64)
	EncodeNew(values []complex128, slots uint64) (plaintext *Plaintext)
	EncodeNTT(plaintext *Plaintext, values []complex128, slots uint64)
	EncodeNTTNew(values []complex128, slots uint64) (plaintext *Plaintext)
	Decode(plaintext *Plaintext, slots uint64) (res []complex128)
	EncodeCoeffs(values []float64, plaintext *Plaintext)
	DecodeCoeffs(plaintext *Plaintext) (res []float64)
}

// encoder is a struct storing the necessary parameters to encode a slice of complex number on a Plaintext.
type encoder struct {
	params       *Parameters
	ckksContext  *Context
	values       []complex128
	valuesfloat  []float64
	bigintCoeffs []*big.Int
	qHalf        *big.Int
	polypool     *ring.Poly
	m            uint64
	roots        []complex128
	rotGroup     []uint64
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a Plaintext.
func NewEncoder(params *Parameters) Encoder {

	if !params.isValid {
		panic("cannot newEncoder: parameters are invalid (check if the generation was done properly)")
	}

	m := uint64(2 << params.LogN)

	rotGroup := make([]uint64, m>>1)
	fivePows := uint64(1)
	for i := uint64(0); i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= GaloisGen
		fivePows &= (m - 1)
	}

	var angle float64
	roots := make([]complex128, m+1)
	for i := uint64(0); i < m; i++ {
		angle = 2 * 3.141592653589793 * float64(i) / float64(m)
		roots[i] = complex(math.Cos(angle), math.Sin(angle))
	}
	roots[m] = roots[0]

	ckksContext := newContext(params)

	return &encoder{
		params:       params.Copy(),
		ckksContext:  ckksContext,
		values:       make([]complex128, m>>2),
		valuesfloat:  make([]float64, m>>1),
		bigintCoeffs: make([]*big.Int, m>>1),
		qHalf:        ring.NewUint(0),
		polypool:     ckksContext.contextQ.NewPoly(),
		m:            m,
		rotGroup:     rotGroup,
		roots:        roots,
	}
}

func (encoder *encoder) EncodeNew(values []complex128, slots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel, encoder.params.Scale)
	encoder.Encode(plaintext, values, slots)
	return
}

// Encode takes a slice of complex128 values of size at most N/2 (the number of slots) and encodes it in the receiver Plaintext.
func (encoder *encoder) Encode(plaintext *Plaintext, values []complex128, slots uint64) {

	if uint64(len(values)) > encoder.ckksContext.maxSlots || uint64(len(values)) > slots {
		panic("cannot Encode: too many values for the given number of slots")
	}

	if slots == 0 && slots&(slots-1) != 0 {
		panic("cannot Encode: slots must be a power of two between 1 and N/2")
	}

	for i := range values{
		encoder.values[i] = values[i]
	}

	encoder.invfft(encoder.values, slots)

	gap := encoder.ckksContext.maxSlots / slots

	for i, jdx, idx := uint64(0), encoder.ckksContext.maxSlots, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx] = real(encoder.values[i])
		encoder.valuesfloat[jdx] = imag(encoder.values[i])
	}

	scaleUpVecExact(encoder.valuesfloat, plaintext.scale, encoder.ckksContext.contextQ.Modulus[:plaintext.Level()+1], plaintext.value.Coeffs)

	for i := uint64(0); i < encoder.ckksContext.maxSlots; i++ {
		encoder.values[i] = 0
	}

	for i := uint64(0); i < encoder.ckksContext.n; i++ {
		encoder.valuesfloat[i] = 0
	}

	plaintext.isNTT = false
}

func (encoder *encoder) EncodeNTTNew(values []complex128, slots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel, encoder.params.Scale)
	encoder.EncodeNTT(plaintext, values, slots)
	return
}

func (encoder *encoder) EncodeNTT(plaintext *Plaintext, values []complex128, slots uint64) {

	encoder.Encode(plaintext, values, slots)

	encoder.ckksContext.contextQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)

	plaintext.isNTT = true
}

// EncodeCoefficients takes as input a polynomial a0 + a1x + a2x^2 + ... + an-1x^n-1 with float coefficient
// and returns a scaled integer plaintext polynomial in NTT.
func (encoder *encoder) EncodeCoeffs(values []float64, plaintext *Plaintext) {

	if uint64(len(values)) > encoder.params.N {
		panic("cannot EncodeCoeffs : too many values (maximum is N)")
	}

	scaleUpVecExact(values, plaintext.scale, encoder.ckksContext.contextQ.Modulus[:plaintext.Level()+1], plaintext.value.Coeffs)

	plaintext.isNTT = false
}

// EncodeCoefficients takes as input a polynomial a0 + a1x + a2x^2 + ... + an-1x^n-1 with float coefficient
// and returns a scaled integer plaintext polynomial in NTT.
func (encoder *encoder) EncodeCoeffsNTT(values []float64, plaintext *Plaintext) {
	encoder.EncodeCoeffs(values, plaintext)
	encoder.ckksContext.contextQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// DecodeCoeffs takes as input a plaintext and returns the scaled down coefficient of the plaintext in flaot64
func (encoder *encoder) DecodeCoeffs(plaintext *Plaintext) (res []float64) {

	if plaintext.isNTT {
		encoder.ckksContext.contextQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	} else {
		encoder.ckksContext.contextQ.CopyLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	}

	res = make([]float64, encoder.params.N)

	// We have more than one moduli and need the CRT reconstruction
	if plaintext.Level() > 0 {

		encoder.ckksContext.contextQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

		Q := encoder.ckksContext.bigintChain[plaintext.Level()]

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

		Q := encoder.ckksContext.contextQ.Modulus[0]
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

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoder) Decode(plaintext *Plaintext, slots uint64) (res []complex128) {

	if plaintext.isNTT {
		encoder.ckksContext.contextQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	} else {
		encoder.ckksContext.contextQ.CopyLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	}

	maxSlots := encoder.ckksContext.maxSlots
	gap := maxSlots / slots

	// We have more than one moduli and need the CRT reconstruction
	if plaintext.Level() > 0 {

		encoder.ckksContext.contextQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

		Q := encoder.ckksContext.bigintChain[plaintext.Level()]

		encoder.qHalf.Set(Q)
		encoder.qHalf.Rsh(encoder.qHalf, 1)

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

			encoder.values[i] = complex(scaleDown(encoder.bigintCoeffs[idx], plaintext.scale), scaleDown(encoder.bigintCoeffs[idx+maxSlots], plaintext.scale))
		}
		// We can directly get the coefficients
	} else {

		Q := encoder.ckksContext.contextQ.Modulus[0]
		coeffs := encoder.polypool.Coeffs[0]

		var real, imag float64
		for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

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

			encoder.values[i] = complex(real, imag) / complex(plaintext.scale, 0)
		}
	}

	encoder.fft(encoder.values, slots)

	res = make([]complex128, slots)

	for i := range res {
		res[i] = encoder.values[i]
	}

	for i := uint64(0); i < encoder.ckksContext.maxSlots; i++ {
		encoder.values[i] = 0
	}

	return
}

func (encoder *encoder) invfftlazy(values []complex128, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v complex128

	for len := N; len >= 1; len >>= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (lenq - (encoder.rotGroup[j] % lenq)) * gap
				u = values[i+j] + values[i+j+lenh]
				v = values[i+j] - values[i+j+lenh]
				v *= encoder.roots[idx]
				values[i+j] = u
				values[i+j+lenh] = v

			}
		}
	}

	sliceBitReverseInPlaceComplex128(values, N)
}

func (encoder *encoder) invfft(values []complex128, N uint64) {

	encoder.invfftlazy(values, N)

	for i := uint64(0); i < N; i++ {
		values[i] /= complex(float64(N), 0)
	}
}

func (encoder *encoder) fft(values []complex128, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v complex128

	sliceBitReverseInPlaceComplex128(values, N)

	for len := uint64(2); len <= N; len <<= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (encoder.rotGroup[j] % lenq) * gap
				u = values[i+j]
				v = values[i+j+lenh]
				v *= encoder.roots[idx]
				values[i+j] = u + v
				values[i+j+lenh] = u - v
			}
		}
	}
}
