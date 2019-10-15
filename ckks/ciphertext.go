package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math/bits"
)

// Ciphertext is a BigPoly of degree > 0.
type Ciphertext struct {
	*ckksElement
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
func (ckkscontext *Context) NewCiphertext(degree uint64, level uint64, scale uint64) *Ciphertext {
	ciphertext := &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextLevel[level].NewPoly()
	}

	ciphertext.scale = scale
	ciphertext.currentModulus = ring.Copy(ckkscontext.contextLevel[level].ModulusBigint)
	ciphertext.isNTT = true

	return ciphertext
}

// NewRandomCiphertext generates a new uniformely distributed ciphertext of degree, level and scale.
func (ckkscontext *Context) NewRandomCiphertext(degree, level, scale uint64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextLevel[level].NewUniformPoly()
	}

	ciphertext.scale = scale
	ciphertext.currentModulus = ring.Copy(ckkscontext.contextLevel[level].ModulusBigint)
	ciphertext.isNTT = true

	return ciphertext
}

// MarshalBinary encodes a ciphertext on a byte slice. The total size
// in byte is 4 + 2 * 8 * N * (level + 1).
func (ciphertext *Ciphertext) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(ciphertext.Value()[0].Coeffs[0]))
	level := uint64(len(ciphertext.Value()[0].Coeffs))
	degree := ciphertext.Degree()
	scale := ciphertext.Scale()

	if level > 0xFF {
		return nil, errors.New("error : ciphertext level overflow uint8")
	}

	if degree > 0xFF {
		return nil, errors.New("error : ciphertext degree overflow uint8")
	}

	if scale > 0xFF {
		return nil, errors.New("error : ciphertext scale overflow uint8")
	}

	data := make([]byte, 4+((N*level*(degree+1))<<3))

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(degree)
	data[3] = uint8(scale)

	pointer := uint64(4)

	for i := uint64(0); i < degree+1; i++ {
		if pointer, err = ring.WriteCoeffsTo(pointer, N, level, ciphertext.Value()[i].Coeffs, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled ciphertext on the target ciphertext.
// The target ciphertext must be of the appropriate format and size, it can be created with the
// method NewCiphertext(uint64, uin64t, uint64).
func (ciphertext *Ciphertext) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	degree := uint64(data[2])

	pointer := uint64(4)

	if ciphertext.Degree() != degree {
		return errors.New("error : invalid ciphertext encoding (unexpected degree)")
	}

	if ciphertext.Level() != level-1 {
		return errors.New("error : invalid ciphertext encoding (unexpected level)")
	}

	if ((uint64(len(data)) - pointer) >> 3) != N*level*(degree+1) {
		return errors.New("error : invalid ciphertext encoding (unexpected data length)")
	}

	ciphertext.SetScale(uint64(data[3]))

	for x := uint64(0); x < degree+1; x++ {
		pointer, _ = ring.DecodeCoeffs(pointer, N, level, ciphertext.Value()[x].Coeffs, data)
	}

	return nil
}
