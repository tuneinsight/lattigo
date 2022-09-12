package ckks

import (
	"encoding/binary"
	"math"

	"github.com/tuneinsight/lattigo/v3/rlwe"
)

type Scale struct {
	Value float64
}

func NewScale(scale float64) *Scale {
	return &Scale{Value: scale}
}

func CheckScaleType(scale rlwe.Scale) *Scale {
	switch scale := scale.(type) {
	case *Scale:
		return scale
	default:
		panic("CheckScaleType: invalid scale.(type), must be *ckks.Scale")
	}
}

func (s *Scale) Set(scale rlwe.Scale) {
	switch scale := scale.(type) {
	case *Scale:
		s.Value = scale.Value
	default:
		panic("ckks.Scale.Set: invalid input type, must be ckks.Scale")
	}
}

func (s *Scale) Div(scale interface{}) {
	switch scale := scale.(type) {
	case float64:
		s.Value /= scale
	case *Scale:
		s.Div(scale.Value)
	default:
		panic("ckks.Scale.Div: invalid input type, must be float64 or *ckks.Scale")
	}
}

func (s *Scale) Mul(scale interface{}) {
	switch scale := scale.(type) {
	case float64:
		s.Value *= scale
	case *Scale:
		s.Mul(scale.Value)
	default:
		panic("ckks.Scale.Mul: invalid input type, must be float64 or *ckks.Scale")
	}
}

func (s *Scale) Max(scale interface{}) (max rlwe.Scale) {
	switch scale := scale.(type) {
	case float64:
		if s.Value <= scale {
			return &Scale{scale}
		}
		return &Scale{s.Value}
	case *Scale:
		return s.Max(scale.Value)
	case rlwe.Scale:
		return s.Max(scale.(*Scale).Value)
	default:
		panic("ckks.Scale.Compare: invalid input type, must be float64 or *ckks.Scale")
	}
}

func (s *Scale) Min(scale interface{}) (min rlwe.Scale) {
	switch scale := scale.(type) {
	case float64:
		if s.Value >= scale {
			return &Scale{scale}
		}
		return &Scale{s.Value}
	case *Scale:
		return s.Min(scale.Value)
	case rlwe.Scale:
		return s.Min(scale.(*Scale).Value)
	default:
		panic("ckks.Scale.Compare: invalid input type, must be float64 or *ckks.Scale")
	}
}

// Compare returns
//	-1 if scale < target
//   0 if scale = target
//   1 if scale > target
func (s *Scale) Compare(scale interface{}) (cmp int) {
	switch scale := scale.(type) {
	case float64:
		if s.Value == scale {
			cmp = 0
		} else if s.Value > scale {
			return 1
		} else {
			return -1
		}
	case *Scale:
		s.Compare(scale.Value)
	default:
		panic("ckks.Scale.Compare: invalid input type, must be float64 or *ckks.Scale")
	}
	return
}

func (s *Scale) Equal(scale interface{}) (cmp bool) {
	switch scale := scale.(type) {
	case float64:
		cmp = s.Value == scale
	case *Scale:
		return s.Equal(scale.Value)
	case rlwe.Scale:
		return s.Equal(scale.(*Scale).Value)
	default:
		panic("ckks.Scale.Compare: invalid input type, must be float64 or *ckks.Scale")
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
	binary.LittleEndian.PutUint64(data, math.Float64bits(s.Value))
}

func (s *Scale) Decode(data []byte) {
	s.Value = math.Float64frombits(binary.LittleEndian.Uint64(data))
}
