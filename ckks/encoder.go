package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math"
)

// Encoder is a struct storing the necessary parameters to encode a slice of complex number on a plaintext.
type Encoder struct {
	ckkscontext   *CkksContext
	values        []complex128
	valuesfloat   []float64
	bigint_coeffs []*ring.Int
	q_half        *ring.Int
	polypool      *ring.Poly
	m             uint64
	roots         []complex128
	rotGroup      []uint64
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a plaintext.
func (ckkscontext *CkksContext) NewEncoder() (encoder *Encoder) {
	encoder = new(Encoder)
	encoder.ckkscontext = ckkscontext
	encoder.values = make([]complex128, ckkscontext.maxSlots)
	encoder.valuesfloat = make([]float64, ckkscontext.n)
	encoder.bigint_coeffs = make([]*ring.Int, ckkscontext.n)
	encoder.q_half = ring.NewUint(0)
	encoder.polypool = ckkscontext.contextQ.NewPoly()

	encoder.m = ckkscontext.n << 1

	encoder.rotGroup = make([]uint64, ckkscontext.n)
	fivePows := uint64(1)
	for i := uint64(0); i < ckkscontext.maxSlots; i++ {
		encoder.rotGroup[i] = fivePows
		fivePows *= 5
		fivePows &= (encoder.m - 1)
	}

	var angle float64
	encoder.roots = make([]complex128, encoder.m+1)
	for i := uint64(0); i < encoder.m; i++ {
		angle = 2 * 3.141592653589793 * float64(i) / float64(encoder.m)
		encoder.roots[i] = complex(math.Cos(angle), math.Sin(angle))
	}
	encoder.roots[encoder.m] = encoder.roots[0]

	return
}

// EncodeFloat takes a slice of complex128 values of size at most N/2 (the number of slots) and encodes it on the receiver plaintext.
func (encoder *Encoder) Encode(plaintext *Plaintext, values []complex128, slots uint64) (err error) {

	if uint64(len(values)) > encoder.ckkscontext.maxSlots || uint64(len(values)) > slots {
		return errors.New("cannot encode -> to many values for the given number of slots")
	}

	plaintext.slots = slots

	for i := uint64(0); i < slots; i++ {
		encoder.values[i] = values[i]
	}

	encoder.invfft(encoder.values, slots)

	gap := encoder.ckkscontext.maxSlots / slots

	for i, jdx, idx := uint64(0), encoder.ckkscontext.maxSlots, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx] = real(encoder.values[i])
		encoder.valuesfloat[jdx] = imag(encoder.values[i])
	}

	scaleUpVecExact(encoder.valuesfloat, plaintext.scale, encoder.ckkscontext.moduli[:plaintext.Level()+1], plaintext.value.Coeffs)

	encoder.ckkscontext.contextQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)

	for i := uint64(0); i < encoder.ckkscontext.maxSlots; i++ {
		encoder.values[i] = 0
	}

	for i := uint64(0); i < encoder.ckkscontext.n; i++ {
		encoder.valuesfloat[i] = 0
	}

	return nil
}

// DecodeFloat decodes the plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *Encoder) Decode(plaintext *Plaintext, slots uint64) (res []complex128) {

	encoder.ckkscontext.contextQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	encoder.ckkscontext.contextQ.PolyToBigint(encoder.polypool, encoder.bigint_coeffs)

	encoder.q_half.SetBigInt(plaintext.currentModulus)
	encoder.q_half.Rsh(encoder.q_half, 1)

	gap := encoder.ckkscontext.maxSlots / slots

	var sign int

	for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

		// Centers the value arounds the current modulus
		encoder.bigint_coeffs[idx].Mod(encoder.bigint_coeffs[idx], plaintext.currentModulus)
		sign = encoder.bigint_coeffs[idx].Compare(encoder.q_half)
		if sign == 1 || sign == 0 {
			encoder.bigint_coeffs[idx].Sub(encoder.bigint_coeffs[idx], plaintext.currentModulus)
		}

		// Centers the value arounds the current modulus
		encoder.bigint_coeffs[idx+encoder.ckkscontext.maxSlots].Mod(encoder.bigint_coeffs[idx+encoder.ckkscontext.maxSlots], plaintext.currentModulus)
		sign = encoder.bigint_coeffs[idx+encoder.ckkscontext.maxSlots].Compare(encoder.q_half)
		if sign == 1 || sign == 0 {
			encoder.bigint_coeffs[idx+encoder.ckkscontext.maxSlots].Sub(encoder.bigint_coeffs[idx+encoder.ckkscontext.maxSlots], plaintext.currentModulus)
		}

		encoder.values[i] = complex(scaleDown(encoder.bigint_coeffs[idx], plaintext.scale), scaleDown(encoder.bigint_coeffs[idx+encoder.ckkscontext.maxSlots], plaintext.scale))
	}

	encoder.fft(encoder.values, slots)

	res = make([]complex128, slots)

	for i := range res {
		res[i] = encoder.values[i]

	}

	for i := uint64(0); i < encoder.ckkscontext.maxSlots; i++ {
		encoder.values[i] = 0
	}

	return
}

func (encoder *Encoder) invfftlazy(values []complex128, N uint64) {

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

	sliceBitReverse64(values, N)
}

func (encoder *Encoder) invfft(values []complex128, N uint64) {

	encoder.invfftlazy(values, N)

	for i := uint64(0); i < N; i++ {
		values[i] /= complex(float64(N), 0)
	}
}

func (encoder *Encoder) fft(values []complex128, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v complex128

	sliceBitReverse64(values, N)

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
