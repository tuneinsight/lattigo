package bgv

import (
	"encoding/binary"
	"errors"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Ciphertext is *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*rlwe.Ciphertext
	Scale uint64
}

// NewCiphertext creates a new Ciphertext parameterized by degree, level and scale.
func NewCiphertext(params Parameters, degree, level int, scale uint64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{Ciphertext: rlwe.NewCiphertext(params.Parameters, degree, level)}
	for _, pol := range ciphertext.Value {
		pol.IsNTT = true
	}
	ciphertext.Scale = scale
	return ciphertext
}

// ScalingFactor returns the scaling factor of the target ciphertext.
func (ct *Ciphertext) ScalingFactor() uint64 {
	return ct.Scale
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params Parameters, degree, level int, scale uint64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{rlwe.NewCiphertextRandom(prng, params.Parameters, degree, level), scale}
	for i := range ciphertext.Value {
		ciphertext.Value[i].IsNTT = true
	}
	return
}

// NewCiphertextAtLevelFromPoly construct a new Ciphetext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Ciphertext will share its backing array of coefficient.
func NewCiphertextAtLevelFromPoly(level int, poly [2]*ring.Poly) *Ciphertext {
	ct := rlwe.NewCiphertextAtLevelFromPoly(level, poly)
	ct.Value[0].IsNTT, ct.Value[1].IsNTT = true, true
	return &Ciphertext{Ciphertext: ct, Scale: 0}
}

// Copy copies the given ciphertext ctp into the receiver ciphertext.
func (ct *Ciphertext) Copy(ctp *Ciphertext) {
	ct.Ciphertext.Copy(ctp.Ciphertext)
	ct.Scale = ctp.Scale
}

// CopyNew makes a deep copy of the receiver ciphertext and returns it.
func (ct *Ciphertext) CopyNew() (ctc *Ciphertext) {
	ctc = &Ciphertext{Ciphertext: ct.Ciphertext.CopyNew(), Scale: ct.Scale}
	return
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ct *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	// MetaData is :
	// 8 byte : Scale
	if WithMetaData {
		dataLen += 8
	}

	dataLen += ct.Ciphertext.GetDataLen(WithMetaData)

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ct *Ciphertext) MarshalBinary() (data []byte, err error) {

	dataScale := make([]byte, 8)

	binary.LittleEndian.PutUint64(dataScale, ct.Scale)

	var dataCt []byte
	if dataCt, err = ct.Ciphertext.MarshalBinary(); err != nil {
		return nil, err
	}

	return append(dataScale, dataCt...), nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
func (ct *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if len(data) < 10 { // cf. ct.GetDataLen()
		return errors.New("too small bytearray")
	}

	ct.Scale = binary.LittleEndian.Uint64(data[0:8])
	ct.Ciphertext = new(rlwe.Ciphertext)
	return ct.Ciphertext.UnmarshalBinary(data[8:])
}
