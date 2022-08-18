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
func (evm *EvalModLiteral) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 26)
	data[0] = uint8(evm.LevelStart)
	binary.BigEndian.PutUint64(data[1:9], math.Float64bits(evm.ScalingFactor))
	data[9] = uint8(evm.SineType)
	binary.BigEndian.PutUint64(data[10:18], math.Float64bits(evm.MessageRatio))
	binary.BigEndian.PutUint32(data[18:22], uint32(evm.K))
	binary.BigEndian.PutUint16(data[22:24], uint16(evm.SineDeg))
	data[24] = uint8(evm.DoubleAngle)
	data[25] = uint8(evm.ArcSineDeg)
	return
}

// UnmarshalBinary decodes a slice of bytes on the target EvalModParameters.
func (evm *EvalModLiteral) UnmarshalBinary(data []byte) (err error) {
	evm.LevelStart = int(data[0])
	evm.ScalingFactor = math.Float64frombits(binary.BigEndian.Uint64(data[1:9]))
	evm.SineType = SineType(int(data[9]))
	evm.MessageRatio = math.Float64frombits(binary.BigEndian.Uint64(data[10:18]))
	evm.K = int(binary.BigEndian.Uint32(data[18:22]))
	evm.SineDeg = int(binary.BigEndian.Uint16(data[22:24]))
	evm.DoubleAngle = int(data[24])
	evm.ArcSineDeg = int(data[25])
	return
}
