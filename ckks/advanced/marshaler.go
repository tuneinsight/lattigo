package advanced

import (
	"encoding/binary"
	"math"
)

func (mParams *EncodingMatrixLiteral) MarshalBinarySize() (dataLen int) {
	dataLen++                      // LinearTransform Type
	dataLen++                      // LogN
	dataLen++                      // LogSlots
	dataLen++                      // LevelStart
	dataLen++                      // #Levels
	dataLen += len(mParams.Levels) // Levels
	dataLen++                      // RepackImag2Real
	dataLen += 8                   // Scaling
	dataLen++                      // BitReversed
	dataLen++                      // LogBSGSRatio
	return
}

// MarshalBinary encode the target EncodingMatrixParameters on a slice of bytes.
func (mParams *EncodingMatrixLiteral) MarshalBinary() (data []byte, err error) {

	data = make([]byte, mParams.MarshalBinarySize())

	var ptr int

	data[ptr] = uint8(mParams.LinearTransformType)
	ptr++

	data[ptr] = uint8(mParams.LogN)
	ptr++

	data[ptr] = uint8(mParams.LogSlots)
	ptr++

	data[ptr] = uint8(mParams.LevelStart)
	ptr++

	data[ptr] = uint8(len(mParams.Levels))
	ptr++

	for _, d := range mParams.Levels {
		data[ptr] = uint8(d)
		ptr++
	}

	if mParams.RepackImag2Real {
		data[ptr] = uint8(1)
	}
	ptr++

	binary.BigEndian.PutUint64(data[ptr:], math.Float64bits(mParams.Scaling))
	ptr += 8

	if mParams.BitReversed {
		data[ptr] = uint8(1)
	}
	ptr++

	data[ptr] = uint8(mParams.LogBSGSRatio)

	return
}

// UnmarshalBinary decodes a slice of bytes on the target EncodingMatrixParameters.
func (mParams *EncodingMatrixLiteral) UnmarshalBinary(data []byte) error {

	var ptr int

	mParams.LinearTransformType = LinearTransformType(int(data[ptr]))
	ptr++

	mParams.LogN = int(data[ptr])
	ptr++

	mParams.LogSlots = int(data[ptr])
	ptr++

	mParams.LevelStart = int(data[ptr])
	ptr++

	mParams.Levels = make([]int, data[ptr])
	ptr++

	for i := range mParams.Levels {
		mParams.Levels[i] = int(data[ptr])
		ptr++
	}

	mParams.RepackImag2Real = data[ptr] == uint8(1)
	ptr++

	mParams.Scaling = math.Float64frombits(binary.BigEndian.Uint64(data[ptr:]))
	ptr += 8

	mParams.BitReversed = data[ptr] == uint8(1)
	ptr++

	mParams.LogBSGSRatio = int(data[ptr])

	return nil
}

// MarshalBinarySize returns the size in bytes of the target EvalModLiteral.
func (evm *EvalModLiteral) MarshalBinarySize() (dataLen int) {
	dataLen++    // LevelStart
	dataLen++    // LogScale
	dataLen++    // SineType
	dataLen++    // LogMessageRatio
	dataLen += 4 // K
	dataLen += 2 // SineDeg
	dataLen++    // DoubleAngle
	dataLen++    // ArcSineDeg
	return
}

// MarshalBinary encode the target EvalModLiteral on a slice of bytes.
func (evm *EvalModLiteral) MarshalBinary() (data []byte, err error) {
	data = make([]byte, evm.MarshalBinarySize())
	data[0] = uint8(evm.LevelStart)
	data[1] = uint8(evm.LogScale)
	data[2] = uint8(evm.SineType)
	data[3] = uint8(evm.LogMessageRatio)
	binary.BigEndian.PutUint32(data[4:8], uint32(evm.K))
	binary.BigEndian.PutUint16(data[8:10], uint16(evm.SineDeg))
	data[10] = uint8(evm.DoubleAngle)
	data[11] = uint8(evm.ArcSineDeg)
	return
}

// UnmarshalBinary decodes a slice of bytes on the target EvalModLiteral.
func (evm *EvalModLiteral) UnmarshalBinary(data []byte) (err error) {
	evm.LevelStart = int(data[0])
	evm.LogScale = int(data[1])
	evm.SineType = SineType(int(data[2]))
	evm.LogMessageRatio = int(data[3])
	evm.K = int(binary.BigEndian.Uint32(data[4:8]))
	evm.SineDeg = int(binary.BigEndian.Uint16(data[8:10]))
	evm.DoubleAngle = int(data[10])
	evm.ArcSineDeg = int(data[11])
	return
}
