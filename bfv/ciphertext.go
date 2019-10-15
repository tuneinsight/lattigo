package bfv

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math/bits"
)

// Ciphertext is a *ring.Poly array representing a polynomial of degree > 0 where coefficients are in R_Q.
type Ciphertext struct {
	*bfvElement
}

// NewCiphertext creates a new empty ciphertext of degree degree.
func (bfvcontext *BfvContext) NewCiphertext(degree uint64) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// NewCiphertextBig creates a new empty ciphertext of degree degree in the extended ciphertext context (Q + P).
func (bfvcontext *BfvContext) NewCiphertextBig(degree uint64) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQP.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// NewRandomCiphertext creates a new ciphertext with uniform coefficients.
func (bfvcontext *BfvContext) NewRandomCiphertext(degree uint64) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewUniformPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// MarshalBinary encodes a ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() ([]byte, error) {

	if ciphertext.IsNTT() {
		return nil, errors.New("cannot marshal ciphertext -> ciphertext is in the NTT domain")
	}

	var err error

	N := uint64(len(ciphertext.Value()[0].Coeffs[0]))
	level := uint64(len(ciphertext.Value()[0].Coeffs))
	degree := ciphertext.Degree()

	if level > 0xFF {
		return nil, errors.New("cannot marshal ciphertext -> ciphertext numberModuli overflow uint8")
	}

	if degree > 0xFF {
		return nil, errors.New("cannot marshal ciphertext -> ciphertext degree overflow uint8")
	}

	data := make([]byte, 4+((N*level*(degree+1))<<3))

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(degree)

	pointer := uint64(3)

	for i := uint64(0); i < degree+1; i++ {
		if pointer, err = ring.WriteCoeffsTo(pointer, N, level, ciphertext.Value()[i].Coeffs, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled ciphertext on the target ciphertext.
// The target ciphertext must be of the appropriate format and size, it can be created with the
// methode NewCiphertext(uint64).
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) error {
	N := uint64(1 << data[0])
	level := uint64(data[1])
	degree := uint64(data[2])

	//Do this to allocate the memory for the underlying bfvElement.
	ciphertext.bfvElement = &bfvElement{
		value: make([]*ring.Poly, degree+1),
		isNTT: false,
	}
	//allocate memory for the coeffs.
	for i := uint64(0); i < degree+1; i++ {
		coeffs := make([][]uint64, level)

		for j := uint64(0); j < level; j++ {
			coeffs[j] = make([]uint64, N)
		}
		ciphertext.value[i] = &ring.Poly{coeffs}

	}

	pointer := uint64(3)

	if ciphertext.Degree() != degree {
		return errors.New("cannot unmarshal ciphertext -> invalid ciphertext encoding (unexpected degree)")
	}

	if uint64(len(ciphertext.Value()[0].Coeffs)) != level {
		return errors.New("cannot unmarshal ciphertext -> invalid ciphertext encoding (unexpected number of moduli)")
	}

	if ((uint64(len(data)) - pointer) >> 3) != N*level*(degree+1) {
		return errors.New("cannot unmarshal ciphertext -> invalid ciphertext encoding (unexpected data length)")
	}

	for x := uint64(0); x < degree+1; x++ {
		pointer, _ = ring.DecodeCoeffs(pointer, N, level, ciphertext.Value()[x].Coeffs, data)
	}

	return nil
}
