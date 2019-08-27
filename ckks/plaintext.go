package ckks

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	"math/bits"
	"math/cmplx"
)

// The plaintext is a ring of N coefficients with two contexts.
// The first context is defined by the BFV parameters. The second
// context defines a NTT around its modulus it it permits it.

type Plaintext BigPoly

// NewPlaintext creates a new plaintext of level level and scale scale.
func (ckkscontext *CkksContext) NewPlaintext(level uint64, scale uint64) *Plaintext {
	plaintext := new(Plaintext)
	plaintext.value = []*ring.Poly{ckkscontext.contextLevel[level].NewPoly()}
	plaintext.scale = scale
	plaintext.currentModulus = ring.Copy(ckkscontext.contextLevel[level].ModulusBigint)
	plaintext.isNTT = true
	return plaintext
}

// Value returns the value (polynomial) of the plaintext.
func (P *Plaintext) Value() []*ring.Poly {
	return P.value
}

// SetValue sets the value (polynomial) of the plaintext to the provided value.
func (P *Plaintext) SetValue(value []*ring.Poly) {
	P.value = value
}

// Resize does nothing on a plaintext since it is always of degree 0.
func (P *Plaintext) Resize(ckkscontext *CkksContext, degree uint64) {

}

// CurrentModulus returns the current modulus of the plaintext.
// This variable is only used during the decoding.
func (P *Plaintext) CurrentModulus() *ring.Int {
	return P.currentModulus
}

// SetCurrentModulus sets the current modulus to the provided values.
// This variable is only used during the decoding.
func (P *Plaintext) SetCurrentModulus(modulus *ring.Int) {
	P.currentModulus = ring.Copy(modulus)
}

// Degree returns the degree of the plaintext,
// this value should always be zero.
func (P *Plaintext) Degree() uint64 {
	return uint64(len(P.value) - 1)
}

// Level returns the current level of the plaintext.
func (P *Plaintext) Level() uint64 {
	return uint64(len(P.value[0].Coeffs) - 1)
}

// Scale returns the current scale of the plaintext (in log2).
func (P *Plaintext) Scale() uint64 {
	return P.scale
}

// SetScale sets the scale of the plaintext to the provided value (in log2).
func (P *Plaintext) SetScale(scale uint64) {
	P.scale = scale
}

// IsNTT returns true or false depending on if the plaintext is in the NTT domain or not.
func (P *Plaintext) IsNTT() bool {
	return P.isNTT
}

// SetIsNTT sets the isNTT value of the plaintext to the provided value.
func (P *Plaintext) SetIsNTT(isNTT bool) {
	P.isNTT = isNTT
}

// NTT applies the NTT transform to a plaintext and returns the result on the receiver element.
// Can only be used if the plaintext is not already in the NTT domain.
func (P *Plaintext) NTT(ckkscontext *CkksContext, ct0 CkksElement) {

	if P.isNTT != true {
		for i := range ct0.Value() {
			ckkscontext.contextLevel[P.Level()].NTT(P.value[i], ct0.Value()[i])
		}
		ct0.SetIsNTT(true)
	}
}

// InvNTT applies the inverse NTT transform to a plaintext and returns the result on the receiver element.
// Can only be used it the plaintext is in the NTT domain
func (P *Plaintext) InvNTT(ckkscontext *CkksContext,ct0 CkksElement) {

	if P.isNTT != false {
		for i := range ct0.Value() {
			ckkscontext.contextLevel[P.Level()].InvNTT(P.value[i], ct0.Value()[i])
		}
		ct0.SetIsNTT(false)
	}
}

// CopyNew creates a new plaintext with the same value and same parameters.
func (P *Plaintext) CopyNew() CkksElement {
	PCopy := new(Plaintext)
	PCopy.value = make([]*ring.Poly, 1)
	PCopy.value[0] = P.value[0].CopyNew()
	P.CopyParams(PCopy)
	return PCopy
}

// Copy copies the value and parameters of the reference plaintext ot the receiver plaintext.
func (P *Plaintext) Copy(PCopy CkksElement) error {
	P.value[0].Copy(PCopy.Value()[0])
	P.CopyParams(PCopy)
	return nil
}

// CopyParams copies the parameters of the reference plaintext to the receiver plaintext.
func (P *Plaintext) CopyParams(ckkselement CkksElement) {
	ckkselement.SetCurrentModulus(P.CurrentModulus())
	ckkselement.SetScale(P.Scale())
	ckkselement.SetIsNTT(P.IsNTT())
}

// EncodeFloat encode a float64 slice of at most N/2 values.
func (plaintext *Plaintext) EncodeFloat(ckkscontext *CkksContext, coeffs []float64) error {

	if len(coeffs) > (len(ckkscontext.indexMatrix)>>1)/int(ckkscontext.gap) {
		return errors.New("error : invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value[0].Coeffs[0]) != len(ckkscontext.indexMatrix) {
		return errors.New("error : invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	values := make([]complex128, len(coeffs)*int(ckkscontext.gap))

	for i := range coeffs {

		values[i*int(ckkscontext.gap)] = complex(coeffs[i], 0)

		for j := 0; j < int(ckkscontext.gap)-1; j++ {
			values[i*int(ckkscontext.gap)+j+1] = complex(0, 0)
		}

	}

	encodeFromComplex(values, plaintext, ckkscontext)

	return nil
}

// EncodeFloat encode a complex128 slice of at most N/2 values.
func (plaintext *Plaintext) EncodeComplex(ckkscontext *CkksContext, coeffs []complex128) error {

	if len(coeffs) > (len(ckkscontext.indexMatrix)>>1)/int(ckkscontext.gap) {
		return errors.New("error : invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value[0].Coeffs[0]) != len(ckkscontext.indexMatrix) {
		return errors.New("error : invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	values := make([]complex128, len(coeffs)*int(ckkscontext.gap))

	for i := range coeffs {
		values[i*int(ckkscontext.gap)] = coeffs[i]

		for j := 0; j < int(ckkscontext.gap)-1; j++ {
			values[i*int(ckkscontext.gap)+j+1] = complex(0, 0)
		}

	}

	encodeFromComplex(values, plaintext, ckkscontext)

	return nil
}

// DecodeFloat decodes the plaintext to a slice of float64 values of size at most N/2.
func (plaintext *Plaintext) DecodeFloat(ckkscontext *CkksContext) (res []float64) {

	values := decodeToComplex(plaintext.value[0], plaintext.currentModulus, ckkscontext.contextLevel[plaintext.Level()], ckkscontext.roots, plaintext.scale)

	res = make([]float64, int(ckkscontext.slots)/int(ckkscontext.gap))

	for i := range res {
		res[i] = real(values[ckkscontext.indexMatrix[i*int(ckkscontext.gap)]])
	}

	return
}

// DecodeFloat decodes the plaintext to a slice of complex128 values of size at most N/2.
func (plaintext *Plaintext) DecodeComplex(ckkscontext *CkksContext) (res []complex128) {

	values := decodeToComplex(plaintext.value[0], plaintext.currentModulus, ckkscontext.contextLevel[plaintext.Level()], ckkscontext.roots, plaintext.scale)

	res = make([]complex128, int(ckkscontext.slots)/int(ckkscontext.gap))

	for i := range res {
		res[i] = values[ckkscontext.indexMatrix[i*int(ckkscontext.gap)]]
	}

	return
}

func encodeFromComplex(coeffs []complex128, plaintext *Plaintext, ckkscontext *CkksContext) {

	values := make([]complex128, ckkscontext.n)

	for i := 0; i < len(coeffs); i++ {
		values[ckkscontext.indexMatrix[i]] = coeffs[i]
		values[ckkscontext.indexMatrix[i+int(ckkscontext.slots)]] = cmplx.Conj(coeffs[i])
	}

	invfft(values, ckkscontext.inv_roots)

	for i, qi := range ckkscontext.modulie {

		for j := uint64(0); j < ckkscontext.n; j++ {

			tmp := real(values[j]) / float64(ckkscontext.n)

			if tmp != 0 {
				plaintext.value[0].Coeffs[i][j] = scaleUp(tmp, plaintext.scale, qi)
			} else {
				plaintext.value[0].Coeffs[i][j] = 0
			}
		}
	}

	ckkscontext.contextLevel[plaintext.Level()].NTT(plaintext.value[0], plaintext.value[0])
}

func decodeToComplex(pol *ring.Poly, Q *ring.Int, context *ring.Context, roots []complex128, scale uint64) (values []complex128) {

	tmp := context.NewPoly()

	context.InvNTT(pol, tmp)

	bigint_coeffs := context.PolyToBigint(tmp)

	Q_half := new(ring.Int)
	Q_half.SetBigInt(Q)
	Q_half.Rsh(Q_half, 1)

	var sign int
	values = make([]complex128, context.N)
	for i := range bigint_coeffs {

		// Centers the value arounds the current modulus
		bigint_coeffs[i].Mod(bigint_coeffs[i], Q)
		sign = bigint_coeffs[i].Compare(Q_half)
		if sign == 1 || sign == 0 {
			bigint_coeffs[i].Sub(bigint_coeffs[i], Q)
		}

		values[i] = complex(scaleDown(bigint_coeffs[i], scale), 0)
	}

	fft(values, roots)

	return
}

func invfft(values, inv_roots []complex128) {

	var logN, mm, k_start, k_end, h, t uint64
	var u, v, psi complex128

	logN = uint64(bits.Len64(uint64(len(values))) - 1)

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
