package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
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
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a plaintext.
func (ckkscontext *CkksContext) NewEncoder() (encoder *Encoder) {
	encoder = new(Encoder)
	encoder.ckkscontext = ckkscontext
	encoder.values = make([]complex128, ckkscontext.n)
	encoder.bigint_coeffs = make([]*ring.Int, ckkscontext.n)
	encoder.q_half = ring.NewUint(0)
	encoder.polypool = ckkscontext.contextLevel[ckkscontext.levels-1].NewPoly()
	return
}

// EncodeFloat takes a slice of float64 values of size at most N/2 (the number of slots) and encodes it on the receiver plaintext.
func (encoder *Encoder) EncodeFloat(plaintext *Plaintext, coeffs []float64) (err error) {

	if len(coeffs) > (len(encoder.ckkscontext.indexMatrix)>>1)/int(encoder.ckkscontext.gap) {
		return errors.New("error : invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value[0].Coeffs[0]) != len(encoder.ckkscontext.indexMatrix) {
		return errors.New("error : invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	preprocessFloat(coeffs, encoder)
	encodeFromComplex(plaintext, encoder)

	return nil
}

func preprocessFloat(coeffs []float64, encoder *Encoder) {

	for i := 0; i < len(coeffs); i++ {

		encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)]] = complex(coeffs[i], 0)
		encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)+int(encoder.ckkscontext.slots)]] = complex(coeffs[i], 0)

		for j := 0; j < int(encoder.ckkscontext.gap)-1; j++ {
			encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)+j+1]] = complex(0, 0)
			encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)+int(encoder.ckkscontext.slots)+j+1]] = complex(0, 0)
		}
	}

	for i := len(coeffs) * int(encoder.ckkscontext.gap); i < int(encoder.ckkscontext.slots); i++ {
		encoder.values[encoder.ckkscontext.indexMatrix[i]] = complex(0, 0)
		encoder.values[encoder.ckkscontext.indexMatrix[i+int(encoder.ckkscontext.slots)]] = complex(0, 0)
	}

	return
}

// EncodeFloat takes a slice of complex128 values of size at most N/2 (the number of slots) and encodes it on the receiver plaintext.
func (encoder *Encoder) EncodeComplex(plaintext *Plaintext, coeffs []complex128) (err error) {

	if len(coeffs) > (len(encoder.ckkscontext.indexMatrix)>>1)/int(encoder.ckkscontext.gap) {
		return errors.New("error : invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value[0].Coeffs[0]) != len(encoder.ckkscontext.indexMatrix) {
		return errors.New("error : invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	preprocessCmplx(coeffs, encoder)
	encodeFromComplex(plaintext, encoder)

	return nil
}

// DecodeFloat decodes the plaintext values to a slice of float64 values of size at most N/2.
func (encoder *Encoder) DecodeFloat(plaintext *Plaintext) (res []float64) {

	decodeToComplex(plaintext, encoder)

	res = make([]float64, int(encoder.ckkscontext.slots)/int(encoder.ckkscontext.gap))

	for i := range res {
		res[i] = real(encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)]])
	}

	return
}

// DecodeFloat decodes the plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *Encoder) DecodeComplex(plaintext *Plaintext) (res []complex128) {

	decodeToComplex(plaintext, encoder)

	res = make([]complex128, int(encoder.ckkscontext.slots)/int(encoder.ckkscontext.gap))

	for i := range res {
		res[i] = encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)]]
	}

	return
}

func encodeFromComplex(plaintext *Plaintext, encoder *Encoder) {

	invfft(encoder.values, encoder.ckkscontext.inv_roots)

	for i, qi := range encoder.ckkscontext.moduli {

		for j := uint64(0); j < encoder.ckkscontext.n; j++ {

			if real(encoder.values[j]) != 0 {
				plaintext.value[0].Coeffs[i][j] = scaleUp(real(encoder.values[j]), plaintext.scale, qi)
			} else {
				plaintext.value[0].Coeffs[i][j] = 0
			}
		}
	}

	encoder.ckkscontext.contextLevel[plaintext.Level()].NTT(plaintext.value[0], plaintext.value[0])
}

func decodeToComplex(plaintext *Plaintext, encoder *Encoder) {

	encoder.ckkscontext.contextLevel[plaintext.Level()].InvNTT(plaintext.value[0], encoder.polypool)
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

	fft(encoder.values, encoder.ckkscontext.roots)

	return
}

func preprocessCmplx(coeffs []complex128, encoder *Encoder) {

	for i := 0; i < len(coeffs); i++ {

		encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)]] = coeffs[i]
		encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)+int(encoder.ckkscontext.slots)]] = cmplx.Conj(coeffs[i])

		for j := 0; j < int(encoder.ckkscontext.gap)-1; j++ {
			encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)+j+1]] = complex(0, 0)
			encoder.values[encoder.ckkscontext.indexMatrix[i*int(encoder.ckkscontext.gap)+int(encoder.ckkscontext.slots)+j+1]] = complex(0, 0)
		}
	}

	for i := len(coeffs) * int(encoder.ckkscontext.gap); i < int(encoder.ckkscontext.slots); i++ {
		encoder.values[encoder.ckkscontext.indexMatrix[i]] = complex(0, 0)
		encoder.values[encoder.ckkscontext.indexMatrix[i+int(encoder.ckkscontext.slots)]] = complex(0, 0)
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
