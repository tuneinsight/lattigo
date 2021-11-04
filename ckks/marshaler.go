package ckks

import (
	"encoding/binary"
	"errors"
	"github.com/ldsec/lattigo/v2/rlwe"
	"math"
)

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	// MetaData is :
	// 8 byte : Scale
	if WithMetaData {
		dataLen += 8
	}

	dataLen += ciphertext.Ciphertext.GetDataLen(WithMetaData)

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	dataScale := make([]byte, 8)

	binary.LittleEndian.PutUint64(dataScale, math.Float64bits(ciphertext.Scale))

	var dataCt []byte
	if dataCt, err = ciphertext.Ciphertext.MarshalBinary(); err != nil {
		return nil, err
	}

	return append(dataScale, dataCt...), nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if len(data) < 10 { // cf. ciphertext.GetDataLen()
		return errors.New("too small bytearray")
	}

	ciphertext.Scale = math.Float64frombits(binary.LittleEndian.Uint64(data[0:8]))
	ciphertext.Ciphertext = new(rlwe.Ciphertext)
	return ciphertext.Ciphertext.UnmarshalBinary(data[8:])
}
