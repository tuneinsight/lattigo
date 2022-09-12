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

func NewScale(params Parameters, scale uint64) *Scale {
	return &Scale{T: params.T(), BredParams: ring.BRedParams(params.T()), Value: scale}
}

func CheckScaleType(scale rlwe.Scale) *Scale {
	switch scale := scale.(type) {
	case *Scale:
		return scale
	default:
		panic("CheckScaleType: invalid scale.(type), must be bgv.Scale")
	}
}

func (s *Scale) Set(scale rlwe.Scale) {
	switch scale := scale.(type) {
	case *Scale:
		s.T = scale.T
		s.BredParams = scale.BredParams[:]
		s.Value = scale.Value
	default:
		panic("bgv.Scale.Set: invalid input type, must be uint64")
	}
}

func (s *Scale) Div(scale interface{}) {
	switch scale := scale.(type) {
	case uint64:
		s.Value = ring.BRed(s.Value, ring.ModExp(scale, s.T-2, s.T), s.T, s.BredParams)
	case *Scale:
		s.Div(scale.Value)
	default:
		panic("bgv.Scale.Div: invalid input type, must be uint64")
	}
}

func (s *Scale) Mul(scale interface{}) {
	switch scale := scale.(type) {
	case uint64:

		s.Value = ring.BRed(s.Value, scale, s.T, s.BredParams)
	case *Scale:
		s.Mul(scale.Value)
	default:
		panic("bgv.Scale.Mul: invalid input type, must be uint64")
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
	case int:
		return s.Equal(uint64(scale))
	case uint64:
		cmp = s.Value == scale
	case *Scale:
		return s.Equal(scale.Value)
	case rlwe.Scale:
		return s.Equal(scale.(*Scale).Value)
	default:
		panic("bgv.Scale.Compare: invalid input type, must be either int, uint64, *bgv.Scale or rlwe.Scale")
	}
	return
}

func (s *Scale) Get() interface{} {
	return s.Value
}

func (s *Scale) CopyNew() rlwe.Scale {
	return &Scale{T: s.T, BredParams: s.BredParams[:], Value: s.Value}
}

func (s *Scale) GetDataLen() int {
	return 32
}

func (s *Scale) Encode(data []byte) {
	binary.LittleEndian.PutUint64(data[0:], s.T)
	binary.LittleEndian.PutUint64(data[8:], s.BredParams[0])
	binary.LittleEndian.PutUint64(data[16:], s.BredParams[1])
	binary.LittleEndian.PutUint64(data[24:], s.Value)
}

func (s *Scale) Decode(data []byte) {
	s.T = binary.LittleEndian.Uint64(data[0:])
	s.BredParams = []uint64{binary.LittleEndian.Uint64(data[8:]), binary.LittleEndian.Uint64(data[16:])}
	s.Value = binary.LittleEndian.Uint64(data[24:])
}
