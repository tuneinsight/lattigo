package bgv

import (
	"encoding/binary"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

type Scale struct {
	T          uint64
	BredParams []uint64
	Value      uint64
}

func NewScale(T uint64) *Scale {
	return &Scale{T: T, BredParams: ring.BRedParams(T), Value: 0}
}

func CheckScaleType(scale rlwe.Scale) *Scale {
	switch scale := scale.(type) {
	case *Scale:
		return scale
	default:
		panic("NewCiphertext: invalid scale.(type), must be ckks.Scale")
	}
}

func (s *Scale) Set(scale interface{}) {
	switch scale := scale.(type) {
	case uint64:
		s.Value = scale
	case *Scale:
		s.Set(scale.Value)
	default:
		panic("ckks.Scale.Set: invalid input type, must be uint64")
	}
}

func (s *Scale) Div(scale interface{}) {
	switch scale := scale.(type) {
	case uint64:
		s.Value = ring.BRed(s.Value, ring.ModExp(scale, s.T-2, s.T), s.T, s.BredParams)
	case *Scale:
		s.Div(scale.Value)
	default:
		panic("ckks.Scale.Div: invalid input type, must be uint64")
	}
}

func (s *Scale) Mul(scale interface{}) {
	switch scale := scale.(type) {
	case uint64:
		s.Value = ring.BRed(s.Value, scale, s.T, s.BredParams)
	case *Scale:
		s.Mul(scale.Value)
	default:
		panic("ckks.Scale.Mul: invalid input type, must be uint64")
	}
}

func (s *Scale) Max(scale interface{}) (max rlwe.Scale) {
	panic("bgv.Scale.Max: method is not supported for this scheme")
}

func (s *Scale) Min(scale interface{}) (min rlwe.Scale) {
	panic("bgv.Scale.Min: method is not supported for this scheme")
}

// Compare returns
//	-1 if scale < target
//   0 if scale = target
//   1 if scale > target
func (s *Scale) Compare(scale interface{}) (cmp int) {
	panic("bgv.Scale.Compare: method is not supported for this scheme")
}

func (s *Scale) Equal(scale interface{}) (cmp bool) {
	switch scale := scale.(type) {
	case uint64:
		cmp = s.Value == scale
	case *Scale:
		return s.Equal(scale.Value)
	case rlwe.Scale:
		return s.Equal(scale.(*Scale).Value)
	default:
		panic("ckks.Scale.Compare: invalid input type, must be uint64")
	}
	return
}

func (s *Scale) Get() interface{} {
	return s.Value
}

func (s *Scale) CopyNew() rlwe.Scale {
	return &Scale{Value: s.Value}
}

func (s *Scale) GetDataLen() int {
	return 8
}

func (s *Scale) Encode(data []byte) {
	binary.LittleEndian.PutUint64(data, s.Value)
}

func (s *Scale) Decode(data []byte) {
	s.Value = binary.LittleEndian.Uint64(data)
}
