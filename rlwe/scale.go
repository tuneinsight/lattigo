package rlwe

import (
	"fmt"
	"math/big"
	"reflect"
)

const (
	// ScalePrecision is the default precision of the scale.
	ScalePrecision = uint(128)
)

// Scale is a struct used to track the scaling factor
// of Plaintext and Ciphertext structs.
// The scale is managed as an 128-bit precision real.
type Scale struct {
	Value big.Float
}

// NewScale instantiates a new Scale.
// Accepted types are int, int64, uint64, float64, *big.Int, *big.Float and *Scale.
// If the input type is not an accepted type, returns an error.
func NewScale(s1 interface{}) (s2 Scale) {
	return Scale{Value: *scaleToBigFloat(s1)}
}

// Float64 returns the underlying scale as a float64 value.
func (s Scale) Float64() float64 {
	f64, _ := s.Value.Float64()
	return f64
}

// Uint64 returns the underlying scale as an uint64 value.
func (s Scale) Uint64() uint64 {
	u64, _ := s.Value.Uint64()
	return u64
}

// Mul multiplies the target s with s2, returning the result in
// a new Scale struct. If mod is specified, performs the multiplication
// modulo mod.
func (s Scale) Mul(s1 Scale, mod ...uint64) Scale {

	res := new(big.Float)

	si, _ := s.Value.Int(nil)
	s1i, _ := s1.Value.Int(nil)
	s1i.Mul(si, s1i)

	if len(mod) > 0 && mod[0] > 0 {
		T := big.NewInt(0).SetUint64(mod[0])
		s1i.Mod(s1i, T)
	}

	res.SetPrec(ScalePrecision)
	res.SetInt(s1i)

	return Scale{Value: *res}
}

// Div multiplies the target s with s1^-1, returning the result in
// a new Scale struct. If mod is specified, performs the multiplication
// modulo t with the multiplicative inverse of s1. Otherwise, performs
// the quotient operation.
func (s Scale) Div(s1 Scale, mod ...uint64) Scale {

	res := new(big.Float)

	if len(mod) > 0 && mod[0] > 0 {
		s1i, _ := s.Value.Int(nil)
		s2i, _ := s1.Value.Int(nil)

		T := big.NewInt(0).SetUint64(mod[0])
		s2i.ModInverse(s2i, T)

		s1i.Mul(s1i, s2i)
		s1i.Mod(s1i, T)

		res.SetPrec(ScalePrecision)
		res.SetInt(s1i)
	} else {
		res.Quo(&s.Value, &s1.Value)
	}

	return Scale{Value: *res}
}

// Cmp compares the target scale with s1.
// Returns 0 if the scales are equal, 1 if
// the target scale is greater and -1 if
// the target scale is smaller.
func (s Scale) Cmp(s1 Scale) (cmp int) {
	return s.Value.Cmp(&s1.Value)
}

// Max returns the a new scale which is the maximum
// between the target scale and s1.
func (s Scale) Max(s1 Scale) (max Scale) {

	if s.Cmp(s1) < 0 {
		return NewScale(s1)
	}

	return NewScale(s)
}

// Min returns the a new scale which is the minimum
// between the target scale and s1.
func (s Scale) Min(s1 Scale) (max Scale) {

	if s.Cmp(s1) > 0 {
		return NewScale(s1)
	}

	return NewScale(s)
}

// MarshalBinarySize returns the size in bytes required to
// encode the target scale.
func (s Scale) MarshalBinarySize() int {
	return 40
}

// Encode encode the target scale on the input slice of bytes.
// If the slice of bytes given as input is smaller than the
// value of .MarshalBinarySize(), the method will return an error.
func (s Scale) Encode(data []byte) (err error) {
	var sBytes []byte
	if sBytes, err = s.Value.MarshalText(); err != nil {
		return
	}

	b := make([]byte, s.MarshalBinarySize())

	if len(data) < len(b) {
		return fmt.Errorf("len(data) < %d", len(b))
	}

	b[0] = uint8(len(sBytes))
	copy(b[1:], sBytes)
	copy(data, b)
	return
}

// Decode decodes the input slice of bytes on the target scale.
// If the input slice of bytes is smaller than .GetDataLen(),
// the method will return an error.
func (s *Scale) Decode(data []byte) (err error) {

	if dLen := s.MarshalBinarySize(); len(data) < dLen {
		return fmt.Errorf("len(data) < %d", dLen)
	}

	bLen := data[0]

	v := new(big.Float)

	if data[1] != 0x30 || bLen > 1 { // 0x30 indicates an empty big.Float
		if err = v.UnmarshalText(data[1 : bLen+1]); err != nil {
			return
		}

		v.SetPrec(ScalePrecision)
	}

	s.Value = *v

	return
}

func scaleToBigFloat(scale interface{}) (s *big.Float) {

	switch scale := scale.(type) {
	case float64:
		if scale < 0 {
			panic(fmt.Errorf("scale cannot be  negative but is %f", scale))
		}
		s = new(big.Float).SetPrec(ScalePrecision)
		s.SetFloat64(scale)
		return
	case *big.Float:
		if scale.Cmp(new(big.Float).SetFloat64(0)) < 0 {
			panic(fmt.Errorf("scale cannot be negative, but is %f", scale))
		}
		s = new(big.Float).SetPrec(ScalePrecision)
		s.Set(scale)
		return
	case big.Float:
		if scale.Cmp(new(big.Float).SetFloat64(0)) < 0 {
			panic(fmt.Errorf("scale cannot be negative, but is %f", &scale))
		}
		s = new(big.Float).SetPrec(ScalePrecision)
		s.Set(&scale)
		return
	case *big.Int:
		if scale.Cmp(new(big.Int).SetInt64(0)) < 0 {
			panic(fmt.Errorf("scale cannot be negative, but is %f", scale))
		}
		s = new(big.Float).SetPrec(ScalePrecision)
		s.SetInt(scale)
		return
	case big.Int:
		if scale.Cmp(new(big.Int).SetInt64(0)) < 0 {
			panic(fmt.Errorf("scale cannot be negative, but is %f", &scale))
		}
		s = new(big.Float).SetPrec(ScalePrecision)
		s.SetInt(&scale)
		return
	case int:
		return scaleToBigFloat(float64(scale))
	case int64:
		return scaleToBigFloat(float64(scale))
	case uint64:
		return scaleToBigFloat(float64(scale))
	case Scale:
		return scaleToBigFloat(scale.Value)
	default:
		panic(fmt.Errorf("invalid scale.(type): must be int, int64, uint64, float64, *big.Int, *big.Float or *Scale but is %s", reflect.TypeOf(scale)))
	}
}
