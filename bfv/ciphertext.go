package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Ciphertext is a *ring.Poly array representing a polynomial of degree > 0 where coefficients are in R_Q.
type Ciphertext struct {
	*bfvElement
}

func NewCiphertext() (ciphertext *Ciphertext) {
	return &Ciphertext{&bfvElement{}}
}

// NewCiphertext creates a new empty ciphertext of degree degree.
func (context *Context) NewCiphertext(degree uint64) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = context.contextQ.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// NewRandomCiphertext creates a new ciphertext with uniform coefficients.
func (context *Context) NewRandomCiphertext(degree uint64) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = context.contextQ.NewUniformPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// MarshalBinary encodes a ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(len(ciphertext.value))
	if ciphertext.isNTT {
		data[1] = 1
	}

	var pointer, inc uint64

	pointer = 2

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

	if uint8(data[1]) == 1 {
		ciphertext.isNTT = true
	}

	var pointer, inc uint64
	pointer = 2

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
		dataLen += 2
	}

	for _, el := range ciphertext.value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}
