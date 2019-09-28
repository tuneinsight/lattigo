package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/bits"
	"math/cmplx"
)

// Encoder is a struct storing the necessary parameters to encode a slice of complex number on a plaintext.
type Encoder struct {
	ckkscontext   *CkksContext
	values        []complex128
	bigint_coeffs []*ring.Int
	q_half        *ring.Int
	polypool      *ring.Poly

	// Encoding and Decoding params
	slots       uint64
	indexMatrix []uint64
	gap         uint64
	roots       []complex128
	inv_roots   []complex128
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a plaintext.
func (ckkscontext *CkksContext) NewEncoder() (encoder *Encoder) {
	encoder = new(Encoder)
	encoder.ckkscontext = ckkscontext
	encoder.values = make([]complex128, ckkscontext.n)
	encoder.bigint_coeffs = make([]*ring.Int, ckkscontext.n)
	encoder.q_half = ring.NewUint(0)
	encoder.polypool = ckkscontext.contextLevel[ckkscontext.levels-1].NewPoly()

	encoder.slots = ckkscontext.slots

	var pos, index1, index2, m, mask uint64

	m = ckkscontext.n << 1

	mask = m - 1

	encoder.gap = 1 //gap-1 is the gap between each slot, here 1 means no gap.

	encoder.indexMatrix = make([]uint64, ckkscontext.n)

	pos = 1

	for i := uint64(0); i < ckkscontext.slots; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		encoder.indexMatrix[i] = bitReverse64(index1, ckkscontext.logN)
		encoder.indexMatrix[i|encoder.slots] = bitReverse64(index2, ckkscontext.logN)

		pos *= ckkscontext.gen
		pos &= mask
	}

	encoder.roots = make([]complex128, ckkscontext.n)
	encoder.inv_roots = make([]complex128, ckkscontext.n)

	angle := 6.283185307179586 / float64(m)
	psi := complex(math.Cos(angle), math.Sin(angle))
	psiInv := complex(1, 0) / psi

	encoder.roots[0] = 1
	encoder.inv_roots[0] = 1

	for j := uint64(1); j < ckkscontext.n; j++ {

		indexReversePrev := bitReverse64(j-1, ckkscontext.logN)
		indexReverseNext := bitReverse64(j, ckkscontext.logN)

		encoder.roots[indexReverseNext] = encoder.roots[indexReversePrev] * psi
		encoder.inv_roots[indexReverseNext] = encoder.inv_roots[indexReversePrev] * psiInv
	}

	return
}

// EncodeFloat takes a slice of float64 values of size at most N/2 (the number of slots) and encodes it on the receiver plaintext.
func (encoder *Encoder) EncodeFloat(plaintext *Plaintext, coeffs []float64) (err error) {

	if len(coeffs) > (len(encoder.indexMatrix)>>1)/int(encoder.gap) {
		return errors.New("error : invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		return errors.New("error : invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	preprocessFloat(coeffs, encoder)
	encodeFromComplex(plaintext, encoder)

	return nil
}

func preprocessFloat(coeffs []float64, encoder *Encoder) {

	for i := 0; i < len(coeffs); i++ {

		encoder.values[encoder.indexMatrix[i*int(encoder.gap)]] = complex(coeffs[i], 0)
		encoder.values[encoder.indexMatrix[i*int(encoder.gap)+int(encoder.slots)]] = complex(coeffs[i], 0)

		for j := 0; j < int(encoder.gap)-1; j++ {
			encoder.values[encoder.indexMatrix[i*int(encoder.gap)+j+1]] = complex(0, 0)
			encoder.values[encoder.indexMatrix[i*int(encoder.gap)+int(encoder.slots)+j+1]] = complex(0, 0)
		}
	}

	for i := len(coeffs) * int(encoder.gap); i < int(encoder.slots); i++ {
		encoder.values[encoder.indexMatrix[i]] = complex(0, 0)
		encoder.values[encoder.indexMatrix[i+int(encoder.slots)]] = complex(0, 0)
	}

	return
}

// EncodeFloat takes a slice of complex128 values of size at most N/2 (the number of slots) and encodes it on the receiver plaintext.
func (encoder *Encoder) EncodeComplex(plaintext *Plaintext, coeffs []complex128) (err error) {

	if len(coeffs) > (len(encoder.indexMatrix)>>1)/int(encoder.gap) {
		return errors.New("error : invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		return errors.New("error : invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	preprocessCmplx(coeffs, encoder)
	encodeFromComplex(plaintext, encoder)

	return nil
}

// DecodeFloat decodes the plaintext values to a slice of float64 values of size at most N/2.
func (encoder *Encoder) DecodeFloat(plaintext *Plaintext) (res []float64) {

	decodeToComplex(plaintext, encoder)

	res = make([]float64, int(encoder.slots)/int(encoder.gap))

	for i := range res {
		res[i] = real(encoder.values[encoder.indexMatrix[i*int(encoder.gap)]])
	}

	return
}

// DecodeFloat decodes the plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *Encoder) DecodeComplex(plaintext *Plaintext) (res []complex128) {

	decodeToComplex(plaintext, encoder)

	res = make([]complex128, int(encoder.slots)/int(encoder.gap))

	for i := range res {
		res[i] = encoder.values[encoder.indexMatrix[i*int(encoder.gap)]]
	}

	return
}

func encodeFromComplex(plaintext *Plaintext, encoder *Encoder) {

	invfft(encoder.values, encoder.inv_roots)

	for i, qi := range encoder.ckkscontext.moduli {

		for j := uint64(0); j < encoder.ckkscontext.n; j++ {

			if real(encoder.values[j]) != 0 {
				plaintext.value.Coeffs[i][j] = scaleUp(real(encoder.values[j]), plaintext.scale, qi)
			} else {
				plaintext.value.Coeffs[i][j] = 0
			}
		}
	}

	encoder.ckkscontext.contextLevel[plaintext.Level()].NTT(plaintext.value, plaintext.value)
}

func decodeToComplex(plaintext *Plaintext, encoder *Encoder) {

	encoder.ckkscontext.contextLevel[plaintext.Level()].InvNTT(plaintext.value, encoder.polypool)
	encoder.ckkscontext.contextLevel[plaintext.Level()].PolyToBigint(encoder.polypool, encoder.bigint_coeffs)

	encoder.q_half.SetBigInt(plaintext.currentModulus)
	encoder.q_half.Rsh(encoder.q_half, 1)

	var sign int

	for i := range encoder.bigint_coeffs {

		// Centers the value arounds the current modulus
		encoder.bigint_coeffs[i].Mod(encoder.bigint_coeffs[i], plaintext.currentModulus)
		sign = encoder.bigint_coeffs[i].Compare(encoder.q_half)
		if sign == 1 || sign == 0 {
			encoder.bigint_coeffs[i].Sub(encoder.bigint_coeffs[i], plaintext.currentModulus)
		}

		encoder.values[i] = complex(scaleDown(encoder.bigint_coeffs[i], plaintext.scale), 0)
	}

	fft(encoder.values, encoder.roots)

	return
}

func preprocessCmplx(coeffs []complex128, encoder *Encoder) {

	for i := 0; i < len(coeffs); i++ {

		encoder.values[encoder.indexMatrix[i*int(encoder.gap)]] = coeffs[i]
		encoder.values[encoder.indexMatrix[i*int(encoder.gap)+int(encoder.slots)]] = cmplx.Conj(coeffs[i])

		for j := 0; j < int(encoder.gap)-1; j++ {
			encoder.values[encoder.indexMatrix[i*int(encoder.gap)+j+1]] = complex(0, 0)
			encoder.values[encoder.indexMatrix[i*int(encoder.gap)+int(encoder.slots)+j+1]] = complex(0, 0)
		}
	}

	for i := len(coeffs) * int(encoder.gap); i < int(encoder.slots); i++ {
		encoder.values[encoder.indexMatrix[i]] = complex(0, 0)
		encoder.values[encoder.indexMatrix[i+int(encoder.slots)]] = complex(0, 0)
	}

	return
}

func invfft(values, inv_roots []complex128) {

	var logN, N, mm, k_start, k_end, h, t uint64
	var u, v, psi complex128

	logN = uint64(bits.Len64(uint64(len(values))) - 1)
	N = 1 << logN

	t = 1

	for i := uint64(0); i < logN; i++ {

		mm = 1 << (logN - i)
		k_start = 0
		h = mm >> 1

		for j := uint64(0); j < h; j++ {

			k_end = k_start + t
			psi = inv_roots[h+j]

			for k := k_start; k < k_end; k++ {

				u = values[k]
				v = values[k+t]
				values[k] = u + v
				values[k+t] = (u - v) * psi
			}

			k_start += (t << 1)
		}

		t <<= 1
	}

	for i := uint64(0); i < N; i++ {
		values[i] /= complex(float64(N), 0)
	}
}

func fft(values, roots []complex128) {

	var logN, t, mm, j1, j2 uint64
	var psi, u, v complex128

	t = uint64(len(values))

	logN = uint64(bits.Len64(t) - 1)

	for i := uint64(0); i < logN; i++ {
		mm = 1 << i
		t >>= 1

		for j := uint64(0); j < mm; j++ {
			j1 = 2 * j * t
			j2 = j1 + t - 1
			psi = roots[mm+j]

			for k := j1; k < j2+1; k++ {

				u = values[k]
				v = values[k+t] * psi
				values[k] = u + v
				values[k+t] = u - v
			}
		}
	}
}

func arrayBitReverse(values []complex128, logN uint64) {

	for i := uint64(0); i < 1<<logN; i++ {
		j := bitReverse64(i, logN)
		if i < j {
			values[i], values[j] = values[j], values[i]
		}
	}
}
