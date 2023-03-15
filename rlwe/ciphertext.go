package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Ciphertext is a generic type for RLWE ciphertexts.
type Ciphertext struct {
	MetaData
	Value []*ring.Poly
}

// NewCiphertext returns a new Ciphertext with zero values and an associated
// MetaData set to the Parameters default value.
func NewCiphertext(params Parameters, degree, level int) (ct *Ciphertext) {
	ct = new(Ciphertext)
	ct.Value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		ct.Value[i] = ring.NewPoly(params.N(), level)
	}
	ct.MetaData = MetaData{Scale: params.defaultScale, IsNTT: params.defaultNTTFlag}
	return
}

// NewCiphertextAtLevelFromPoly constructs a new Ciphertext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Ciphertext will share its backing array of coefficients.
// Returned Ciphertext's MetaData is empty.
func NewCiphertextAtLevelFromPoly(level int, poly []*ring.Poly) (ct *Ciphertext) {
	Value := make([]*ring.Poly, len(poly))
	for i := range Value {
		Value[i] = new(ring.Poly)
		Value[i].Coeffs = poly[i].Coeffs[:level+1]
		Value[i].Buff = poly[i].Buff[:poly[i].N()*(level+1)]
	}
	return &Ciphertext{Value: Value}
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level.
func NewCiphertextRandom(prng utils.PRNG, params Parameters, degree, level int) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(params, degree, level)
	PopulateElementRandom(prng, params, ciphertext)
	return
}

// Degree returns the degree of the target Ciphertext.
func (ct *Ciphertext) Degree() int {
	return len(ct.Value) - 1
}

// Level returns the level of the target Ciphertext.
func (ct *Ciphertext) Level() int {
	return len(ct.Value[0].Coeffs) - 1
}

// GetScale gets the scale of the target ciphertext.
func (ct *Ciphertext) GetScale() Scale {
	return ct.Scale
}

// SetScale sets the scale of the target ciphertext.
func (ct *Ciphertext) SetScale(scale Scale) {
	ct.Scale = scale
}

// Resize resizes the degree of the target element.
// Sets the NTT flag of the added poly equal to the NTT flag
// to the poly at degree zero.
func (ct *Ciphertext) Resize(degree, level int) {

	if ct.Level() != level {
		for i := range ct.Value {
			ct.Value[i].Resize(level)
		}
	}

	if ct.Degree() > degree {
		ct.Value = ct.Value[:degree+1]
	} else if ct.Degree() < degree {
		for ct.Degree() < degree {
			ct.Value = append(ct.Value, []*ring.Poly{ring.NewPoly(ct.Value[0].N(), level)}...)
		}
	}
}

// SwitchCiphertextRingDegreeNTT changes the ring degree of ctIn to the one of ctOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of ctOut is larger than the one of ctIn, then the ringQ of ctOut
// must be provided (otherwise, a nil pointer).
// The ctIn must be in the NTT domain and ctOut will be in the NTT domain.
func SwitchCiphertextRingDegreeNTT(ctIn *Ciphertext, ringQLargeDim *ring.Ring, ctOut *Ciphertext) {

	NIn, NOut := len(ctIn.Value[0].Coeffs[0]), len(ctOut.Value[0].Coeffs[0])

	if NIn > NOut {

		gap := NIn / NOut
		buff := make([]uint64, NIn)
		for i := range ctOut.Value {
			for j := range ctOut.Value[i].Coeffs {

				tmpIn, tmpOut := ctIn.Value[i].Coeffs[j], ctOut.Value[i].Coeffs[j]

				ringQLargeDim.SubRings[j].INTT(tmpIn, buff)

				for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+1, w1+gap {
					tmpOut[w0] = buff[w1]
				}

				s := ringQLargeDim.SubRings[j]

				switch ringQLargeDim.Type() {
				case ring.Standard:
					ring.NTTStandard(tmpOut, tmpOut, NOut, s.Modulus, s.MRedConstant, s.BRedConstant, s.RootsForward)
				case ring.ConjugateInvariant:
					ring.NTTConjugateInvariant(tmpOut, tmpOut, NOut, s.Modulus, s.MRedConstant, s.BRedConstant, s.RootsForward)
				}
			}
		}

	} else {
		for i := range ctOut.Value {
			ring.MapSmallDimensionToLargerDimensionNTT(ctIn.Value[i], ctOut.Value[i])
		}
	}

	ctOut.MetaData = ctIn.MetaData
}

// SwitchCiphertextRingDegree changes the ring degree of ctIn to the one of ctOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of ctOut is larger than the one of ctIn, then the ringQ of ctIn
// must be provided (otherwise, a nil pointer).
func SwitchCiphertextRingDegree(ctIn *Ciphertext, ctOut *Ciphertext) {

	NIn, NOut := len(ctIn.Value[0].Coeffs[0]), len(ctOut.Value[0].Coeffs[0])

	gapIn, gapOut := NOut/NIn, 1
	if NIn > NOut {
		gapIn, gapOut = 1, NIn/NOut
	}

	for i := range ctOut.Value {
		for j := range ctOut.Value[i].Coeffs {
			tmp0, tmp1 := ctOut.Value[i].Coeffs[j], ctIn.Value[i].Coeffs[j]
			for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+gapIn, w1+gapOut {
				tmp0[w0] = tmp1[w1]
			}
		}
	}

	ctOut.MetaData = ctIn.MetaData
}

// CopyNew creates a new element as a copy of the target element.
func (ct *Ciphertext) CopyNew() *Ciphertext {

	ctxCopy := new(Ciphertext)

	ctxCopy.Value = make([]*ring.Poly, ct.Degree()+1)
	for i := range ct.Value {
		ctxCopy.Value[i] = ct.Value[i].CopyNew()
	}

	ctxCopy.MetaData = ct.MetaData

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (ct *Ciphertext) Copy(ctxCopy *Ciphertext) {

	if ct != ctxCopy {
		for i := range ctxCopy.Value {
			ct.Value[i].Copy(ctxCopy.Value[i])
		}

		ct.MetaData = ctxCopy.MetaData
	}
}

// El returns a pointer to this Element
func (ct *Ciphertext) El() *Ciphertext {
	return ct
}

// GetSmallestLargest returns the provided element that has the smallest degree as a first
// returned value and the largest degree as second return value. If the degree match, the
// order is the same as for the input.
func GetSmallestLargest(el0, el1 *Ciphertext) (smallest, largest *Ciphertext, sameDegree bool) {
	switch {
	case el0.Degree() > el1.Degree():
		return el1, el0, false
	case el0.Degree() < el1.Degree():
		return el0, el1, false
	}
	return el0, el1, true
}

// PopulateElementRandom creates a new rlwe.Element with random coefficients.
func PopulateElementRandom(prng utils.PRNG, params Parameters, ct *Ciphertext) {
	sampler := ring.NewUniformSampler(prng, params.RingQ()).AtLevel(ct.Level())
	for i := range ct.Value {
		sampler.Read(ct.Value[i])
	}
}

// MarshalBinarySize returns the length in bytes of the target Ciphertext.
func (ct *Ciphertext) MarshalBinarySize() (dataLen int) {

	// 1 byte : Degree
	dataLen++

	for _, ct := range ct.Value {
		dataLen += ct.MarshalBinarySize64()
	}

	dataLen += ct.MetaData.MarshalBinarySize()

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in bytes is 4 + 8* N * numberModuliQ * (degree + 1).
func (ct *Ciphertext) MarshalBinary() (data []byte, err error) {
	data = make([]byte, ct.MarshalBinarySize())
	_, err = ct.MarshalBinaryInPlace(data)
	return
}

// MarshalBinaryInPlace encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (ct *Ciphertext) MarshalBinaryInPlace(data []byte) (ptr int, err error) {

	if len(data) < ct.MarshalBinarySize() {
		return 0, fmt.Errorf("cannot write: len(data) is too small")
	}

	if ptr, err = ct.MetaData.MarshalBinaryInPlace(data); err != nil {
		return
	}

	data[ptr] = uint8(ct.Degree() + 1)
	ptr++

	var inc int
	for _, pol := range ct.Value {

		if inc, err = pol.Encode64(data[ptr:]); err != nil {
			return
		}

		ptr += inc
	}

	return
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
func (ct *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	_, err = ct.UnmarshalBinaryInPlace(data)
	return
}

// UnmarshalBinaryInPlace decodes a slice of bytes generated by MarshalBinary or
// MarshalBinaryInPlace on the object and returns the number of bytes read.
func (ct *Ciphertext) UnmarshalBinaryInPlace(data []byte) (ptr int, err error) {

	if ptr, err = ct.MetaData.UnmarshalBinaryInPlace(data); err != nil {
		return
	}

	if degree := int(data[ptr]); ct.Value == nil {
		ct.Value = make([]*ring.Poly, degree)
	} else {
		if len(ct.Value) > degree {
			ct.Value = ct.Value[:degree]
		} else {
			ct.Value = append(ct.Value, make([]*ring.Poly, degree-len(ct.Value))...)
		}
	}
	ptr++

	var inc int
	for i := range ct.Value {

		if ct.Value[i] == nil {
			ct.Value[i] = new(ring.Poly)
		}

		if inc, err = ct.Value[i].Decode64(data[ptr:]); err != nil {
			return
		}

		ptr += inc
	}

	return
}
