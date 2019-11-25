package ckks

import (
	"encoding/binary"
	"github.com/ldsec/lattigo/ring"
	"math"
)

// Ciphertext is a BigPoly of degree > 0.
type Ciphertext struct {
	*ckksElement
}

func NewCiphertext() (ciphertext *Ciphertext) {
	return &Ciphertext{&ckksElement{}}
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
func (ckkscontext *Context) NewCiphertext(degree uint64, level uint64, scale float64) *Ciphertext {
	ciphertext := &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextQ.NewPolyLvl(level)
	}

	ciphertext.scale = scale
	ciphertext.isNTT = true

	return ciphertext
}

// NewRandoMCiphertext generates a new uniformely distributed ciphertext of degree, level and scale.
func (ckkscontext *Context) NewRandomCiphertext(degree, level uint64, scale float64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextQ.NewUniformPolyLvl(level)
	}

	ciphertext.scale = scale
	ciphertext.isNTT = true

	return ciphertext
}

// MarshalBinary encodes a ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(ciphertext.Degree() + 1)

	binary.LittleEndian.PutUint64(data[1:9], math.Float64bits(ciphertext.Scale()))

	if ciphertext.isNTT {
		data[10] = 1
	}

	var pointer, inc uint64

	pointer = 11

	for _, el := range ciphertext.value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled ciphertext on the target ciphertext.
// The target ciphertext must be of the appropriate format and size, it can be created with the
// methode NewCiphertext(uint64).
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {

	ciphertext.value = make([]*ring.Poly, uint8(data[0]))

	ciphertext.scale = math.Float64frombits(binary.LittleEndian.Uint64(data[1:9]))

	if uint8(data[10]) == 1 {
		ciphertext.isNTT = true
	}

	var pointer, inc uint64
	pointer = 11

	for i := range ciphertext.value {

		ciphertext.value[i] = new(ring.Poly)

		if inc, err = ciphertext.value[i].DecodePolyNew(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	return nil
}

func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen uint64) {
	if WithMetaData {
		dataLen += 11
	}

	for _, el := range ciphertext.value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}
