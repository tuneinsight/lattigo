package advanced

import (
	"encoding/binary"
	"math"
)

// MarshalBinary encode the target EncodingMatrixParameters on a slice of bytes.
func (mParams *EncodingMatrixLiteral) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 9)
	data[0] = uint8(mParams.LinearTransformType)
	data[1] = uint8(mParams.LevelStart)
	if mParams.BitReversed {
		data[2] = uint8(1)
	}

	binary.BigEndian.PutUint32(data[3:7], math.Float32bits(float32(mParams.BSGSRatio)))
	data[7] = uint8(len(mParams.ScalingFactor))

	if mParams.RepackImag2Real {
		data[8] = uint8(1)
	}

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

// UnmarshalBinary decodes a slice of bytes on the target EncodingMatrixParameters.
func (mParams *EncodingMatrixLiteral) UnmarshalBinary(data []byte) error {

	mParams.LinearTransformType = LinearTransformType(int(data[0]))
	mParams.LevelStart = int(data[1])
	if data[2] == 1 {
		mParams.BitReversed = true
	}
	mParams.BSGSRatio = float64(math.Float32frombits(binary.BigEndian.Uint32(data[3:7])))
	mParams.ScalingFactor = make([][]float64, data[7])
	mParams.RepackImag2Real = data[8] == 1

	pt := 9
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
func (evmParams *EvalModLiteral) MarshalBinary() (data []byte, err error) {
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
func (evmParams *EvalModLiteral) UnmarshalBinary(data []byte) (err error) {
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
