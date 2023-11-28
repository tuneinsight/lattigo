package rlwe

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo-enterprise/v5/utils/bignum"
)

const (
	// ScalePrecision is the default precision of the scale.
	ScalePrecision = uint(128)
)

var ScalePrecisionLog10 = int(math.Ceil(float64(ScalePrecision) / math.Log2(10)))

// Scale is a struct used to track the scaling factor
// of Plaintext and Ciphertext structs.
// The scale is managed as an 128-bit precision real and can
// be either a floating point value or a mod T
// prime integer, which is determined at instantiation.
type Scale struct {
	Value big.Float //`json:",omitempty"`
	Mod   *big.Int  //`json:",omitempty"`
}

// NewScale instantiates a new floating point Scale.
// Accepted types for s are int, int64, uint64, float64, *big.Int, *big.Float and *Scale.
// If the input type is not an accepted type, returns an error.
func NewScale(s interface{}) Scale {
	v := scaleToBigFloat(s)
	return Scale{Value: *v}
}

// NewScaleModT instantiates a new integer mod T Scale.
// Accepted types for s are int, int64, uint64, float64, *big.Int, *big.Float and *Scale.
// If the input type is not an accepted type, returns an error.
func NewScaleModT(s interface{}, mod uint64) Scale {
	scale := NewScale(s)
	if mod != 0 {
		scale.Mod = big.NewInt(0).SetUint64(mod)
	}
	return scale
}

// BigInt returns the scale as a big.Int, truncating the rational part and rounding ot the nearest integer.
// The rounding assumes that the scale is a positive value.
func (s Scale) BigInt() (sInt *big.Int) {
	sInt = new(big.Int)
	new(big.Float).SetPrec(s.Value.Prec()).Add(&s.Value, new(big.Float).SetFloat64(0.5)).Int(sInt)
	return
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

// Mul multiplies the target s with s1, returning the result in
// a new Scale struct. If mod is specified, performs the multiplication
// modulo mod.
func (s Scale) Mul(s1 Scale) Scale {

	res := new(big.Float)

	if s.Mod != nil {
		si, _ := s.Value.Int(nil)
		s1i, _ := s1.Value.Int(nil)
		s1i.Mul(si, s1i)
		s1i.Mod(s1i, s.Mod)
		res.SetPrec(ScalePrecision)
		res.SetInt(s1i)
	} else {
		res.Mul(&s.Value, &s1.Value)
	}

	return Scale{Value: *res, Mod: s.Mod}
}

// Div multiplies the target s with s1^-1, returning the result in
// a new Scale struct. If mod is specified, performs the multiplication
// modulo t with the multiplicative inverse of s1. Otherwise, performs
// the quotient operation.
func (s Scale) Div(s1 Scale) Scale {

	res := new(big.Float)

	if s.Mod != nil {
		s1i, _ := s.Value.Int(nil)
		s2i, _ := s1.Value.Int(nil)

		s2i.ModInverse(s2i, s.Mod)

		s1i.Mul(s1i, s2i)
		s1i.Mod(s1i, s.Mod)

		res.SetPrec(ScalePrecision)
		res.SetInt(s1i)
	} else {
		res.Quo(&s.Value, &s1.Value)
	}

	return Scale{Value: *res, Mod: s.Mod}
}

// Cmp compares the target scale with s1.
// Returns 0 if the scales are equal, 1 if
// the target scale is greater and -1 if
// the target scale is smaller.
func (s Scale) Cmp(s1 Scale) (cmp int) {
	return s.Value.Cmp(&s1.Value)
}

// Equal returns true if a == b.
func (s Scale) Equal(s1 Scale) bool {
	return s.Cmp(s1) == 0
}

// InDelta returns true if abs(a-b) <= 2^{-log2Delta}
func (s Scale) InDelta(s1 Scale, log2Delta float64) bool {
	return s.Log2Delta(s1) >= log2Delta
}

// Log2Delta returns -log2(abs(a-b)/max(a, b))
func (s Scale) Log2Delta(s1 Scale) float64 {
	d := new(big.Float).Sub(&s.Value, &s1.Value)
	d.Abs(d)
	max := s.Max(s1)
	d.Quo(d, &max.Value)
	d.Quo(bignum.Log(d), bignum.Log2(s.Value.Prec()))
	d.Neg(d)
	f64, _ := d.Float64()
	return f64
}

// Max returns the a new scale which is the maximum
// between the target scale and s1.
func (s Scale) Max(s1 Scale) (max Scale) {

	if s.Cmp(s1) < 0 {
		return s1
	}

	return s
}

// Min returns the a new scale which is the minimum
// between the target scale and s1.
func (s Scale) Min(s1 Scale) (max Scale) {

	if s.Cmp(s1) > 0 {
		return s1
	}

	return s
}

// BinarySize returns the serialized size of the object in bytes.
// Each value is encoded with .Text('x', ceil(ScalePrecision / log2(10))).
func (s Scale) BinarySize() int {
	return 21 + (ScalePrecisionLog10+8)<<1 // 21 for JSON formatting and 2*(8 + ScalePrecisionLog10)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (s Scale) MarshalBinary() (p []byte, err error) {
	return s.MarshalJSON()
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (s Scale) UnmarshalBinary(p []byte) (err error) {
	return s.UnmarshalJSON(p)
}

// MarshalJSON encodes the object into a binary form on a newly allocated slice of bytes.
func (s Scale) MarshalJSON() (p []byte, err error) {

	var mod string

	if s.Mod != nil {
		mod = new(big.Float).SetPrec(ScalePrecision).SetInt(s.Mod).Text('x', ScalePrecisionLog10)
	} else {

		var m string
		for i := 0; i < ScalePrecisionLog10; i++ {
			m += "0"
		}

		mod = "0x0." + m + "p+00"
	}

	aux := &struct {
		Value string
		Mod   string
	}{
		Value: s.Value.Text('x', ScalePrecisionLog10),
		Mod:   mod,
	}

	p, err = json.Marshal(aux)

	return
}

func (s *Scale) UnmarshalJSON(p []byte) (err error) {

	aux := &struct {
		Value string
		Mod   string
	}{}

	if err = json.Unmarshal(p, aux); err != nil {
		return
	}

	s.Value.SetPrec(ScalePrecision)

	f, _ := new(big.Float).SetString(aux.Value)
	s.Value.Set(f)

	mod, bool := new(big.Float).SetString(aux.Mod)

	if mod.Cmp(new(big.Float)) != 0 {

		if s.Mod == nil {
			s.Mod = new(big.Int)
		}

		if !bool {
			return fmt.Errorf("Scale: UnmarshalJSON: s.Mod != exact")
		}

		mod.Int(s.Mod)
	}

	return
}

func scaleToBigFloat(scale interface{}) (s *big.Float) {

	switch scale := scale.(type) {
	case float64:
		if scale < 0 || math.IsNaN(scale) || math.IsInf(scale, 0) {
			panic(fmt.Errorf("scale cannot be negative, NaN or Inf, but is %f", scale))
		}

		s = new(big.Float).SetPrec(ScalePrecision)
		s.SetFloat64(scale)
		return
	case *big.Float:
		if scale.Cmp(new(big.Float).SetFloat64(0)) < 0 || scale.IsInf() {
			panic(fmt.Errorf("scale cannot be negative, but is %f", scale))
		}
		s = new(big.Float).SetPrec(ScalePrecision)
		s.Set(scale)
		return
	case big.Float:
		if scale.Cmp(new(big.Float).SetFloat64(0)) < 0 || scale.IsInf() {
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
		return scaleToBigFloat(new(big.Int).SetInt64(int64(scale)))
	case int64:
		return scaleToBigFloat(new(big.Int).SetInt64(scale))
	case uint64:
		return scaleToBigFloat(new(big.Int).SetUint64(scale))
	case Scale:
		return scaleToBigFloat(scale.Value)
	default:
		panic(fmt.Errorf("invalid scale.(type): must be int, int64, uint64, float64, *big.Int, *big.Float or *Scale but is %T", scale))
	}
}
