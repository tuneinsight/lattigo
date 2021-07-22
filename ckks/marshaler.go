package ckks

import (
	"encoding/binary"
	"errors"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"math"
)

// MarshalBinary encode the target EncodingMatricesParameters on a slice of bytes.
func (mParams *EncodingMatricesParameters) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 8)
	data[0] = uint8(mParams.LinearTransformType)
	data[1] = uint8(mParams.LevelStart)
	if mParams.BitReversed {
		data[2] = uint8(1)
	}

	binary.BigEndian.PutUint32(data[3:7], math.Float32bits(float32(mParams.BSGSRatio)))
	data[7] = uint8(len(mParams.ScalingFactor))

	for _, d := range mParams.ScalingFactor {
		data = append(data, uint8(len(d)))
		for j := range d {
			tmp := make([]byte, 8)
			binary.BigEndian.PutUint64(tmp, math.Float64bits(float64(d[j])))
			data = append(data, tmp...)
		}
	}
	return
}

// UnmarshalBinary decodes a slice of bytes on the target EncodingMatricesParameters.
func (mParams *EncodingMatricesParameters) UnmarshalBinary(data []byte) error {

	mParams.LinearTransformType = LinearTransformType(int(data[0]))
	mParams.LevelStart = int(data[1])
	if data[2] == 1 {
		mParams.BitReversed = true
	}
	mParams.BSGSRatio = float64(math.Float32frombits(binary.BigEndian.Uint32(data[3:7])))

	mParams.ScalingFactor = make([][]float64, data[7])
	pt := 8
	for i := range mParams.ScalingFactor {
		tmp := make([]float64, data[pt])
		pt++
		for j := range tmp {
			tmp[j] = math.Float64frombits(binary.BigEndian.Uint64(data[pt : pt+8]))
			pt += 8
		}
		mParams.ScalingFactor[i] = tmp
	}

	return nil
}

// MarshalBinary encode the target EvalModParameters on a slice of bytes.
func (evmParams *EvalModParameters) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 35)
	binary.BigEndian.PutUint64(data[:8], evmParams.Q)
	data[8] = uint8(evmParams.LevelStart)
	binary.BigEndian.PutUint64(data[9:17], math.Float64bits(evmParams.ScalingFactor))
	data[17] = uint8(evmParams.SineType)
	binary.BigEndian.PutUint64(data[18:26], math.Float64bits(evmParams.MessageRatio))
	binary.BigEndian.PutUint32(data[26:30], uint32(evmParams.K))
	binary.BigEndian.PutUint16(data[30:32], uint16(evmParams.SineDeg))
	data[33] = uint8(evmParams.DoubleAngle)
	data[34] = uint8(evmParams.ArcSineDeg)
	return
}

// UnmarshalBinary decodes a slice of bytes on the target EvalModParameters.
func (evmParams *EvalModParameters) UnmarshalBinary(data []byte) (err error) {
	evmParams.Q = binary.BigEndian.Uint64(data[:8])
	evmParams.LevelStart = int(data[8])
	evmParams.ScalingFactor = math.Float64frombits(binary.BigEndian.Uint64(data[9:17]))
	evmParams.SineType = SineType(int(data[17]))
	evmParams.MessageRatio = math.Float64frombits(binary.BigEndian.Uint64(data[18:26]))
	evmParams.K = int(binary.BigEndian.Uint32(data[26:30]))
	evmParams.SineDeg = int(binary.BigEndian.Uint16(data[30:32]))
	evmParams.DoubleAngle = int(data[33])
	evmParams.ArcSineDeg = int(data[34])
	return
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	// MetaData is :
	// 1 byte : Degree
	// 9 byte : Scale
	// 1 byte : isNTT
	if WithMetaData {
		dataLen += 10
	}

	for _, el := range ciphertext.Value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(ciphertext.Degree() + 1)

	binary.LittleEndian.PutUint64(data[1:9], math.Float64bits(ciphertext.Scale))

	var pointer, inc int

	pointer = 10

	for _, el := range ciphertext.Value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if len(data) < 10 { // cf. ciphertext.GetDataLen()
		return errors.New("too small bytearray")
	}

	ciphertext.Ciphertext = new(rlwe.Ciphertext)

	ciphertext.Value = make([]*ring.Poly, uint8(data[0]))

	ciphertext.Scale = math.Float64frombits(binary.LittleEndian.Uint64(data[1:9]))

	var pointer, inc int
	pointer = 10

	for i := range ciphertext.Value {

		ciphertext.Value[i] = new(ring.Poly)

		if inc, err = ciphertext.Value[i].DecodePolyNew(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	if pointer != len(data) {
		return errors.New("remaining unparsed data")
	}

	return nil
}
