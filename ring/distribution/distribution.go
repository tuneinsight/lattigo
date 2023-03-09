// Package distribution implements definition for sampling distributions.
package distribution

import (
	"encoding/binary"
	"encoding/json"
	"fmt"
	"math"
)

type Type uint8

const (
	uniform Type = iota + 1
	ternary
	discreteGaussian
)

var typeToString = [5]string{"Undefined", "Uniform", "Ternary", "DiscreteGaussian"}

var typeFromString = map[string]Type{
	"Undefined":        0,
	"Uniform":          uniform,
	"Ternary":          ternary,
	"DiscreteGaussian": discreteGaussian,
}

func (t Type) String() string {
	if int(t) >= len(typeToString) {
		return "Unknown"
	}
	return typeToString[int(t)]
}

// Distribution is a interface for distributions
type Distribution interface {
	Type() Type
	StandardDeviation(LogN int, LogQP float64) StandardDeviation // TODO: properly define
	Equals(Distribution) bool
	CopyNew() Distribution

	MarshalBinarySize() int
	Encode(data []byte) (ptr int, err error)
	Decode(data []byte) (ptr int, err error)
}

func NewFromMap(distDef map[string]interface{}) (Distribution, error) {
	distTypeVal, specified := distDef["Type"]
	if !specified {
		return nil, fmt.Errorf("map specifies no distribution type")
	}
	distTypeStr, isString := distTypeVal.(string)
	if !isString {
		return nil, fmt.Errorf("value for key Type of map should be of type string")
	}
	distType, exists := typeFromString[distTypeStr]
	if !exists {
		return nil, fmt.Errorf("distribution type %s does not exist", distTypeStr)
	}
	switch distType {
	case uniform:
		return NewUniform(distDef)
	case ternary:
		return NewTernary(distDef)
	case discreteGaussian:
		return NewDiscreteGaussian(distDef)
	default:
		return nil, fmt.Errorf("invalid distribution type")
	}
}

func Encode(X Distribution, data []byte) (ptr int, err error) {
	if len(data) == 1+X.MarshalBinarySize() {
		return 0, fmt.Errorf("buffer is too small for encoding distribution (size %d instead of %d)", len(data), 1+X.MarshalBinarySize())
	}
	data[0] = byte(X.Type())
	ptr, err = X.Encode(data[1:])

	return ptr + 1, err
}

func Decode(data []byte) (ptr int, X Distribution, err error) {
	if len(data) == 0 {
		return 0, nil, fmt.Errorf("data should have length >= 1")
	}
	switch Type(data[0]) {
	case uniform:
		X = &Uniform{}
	case ternary:
		X = &Ternary{}
	case discreteGaussian:
		X = &DiscreteGaussian{}
	default:
		return 0, nil, fmt.Errorf("invalid distribution type: %s", Type(data[0]))
	}

	ptr, err = X.Decode(data[1:])

	return ptr + 1, X, err
}

// StandardDeviation is a float64 type storing
// a value representing a standard deviation
type StandardDeviation float64

// DiscreteGaussian is a discrete Gaussian distribution
// with a given standard deviation and a bound
// in number of standard deviations.
type DiscreteGaussian struct {
	Sigma StandardDeviation
	Bound int
}

func NewDiscreteGaussian(distDef map[string]interface{}) (d *DiscreteGaussian, err error) {
	sigma, errSigma := getFloatFromMap(distDef, "Sigma")
	if errSigma != nil {
		return nil, err
	}
	bound, errBound := getIntFromMap(distDef, "Bound")
	if errBound != nil {
		return nil, err
	}
	return &DiscreteGaussian{Sigma: StandardDeviation(sigma), Bound: bound}, nil
}

func (d *DiscreteGaussian) Type() Type {
	return discreteGaussian
}

func (d *DiscreteGaussian) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return StandardDeviation(d.Sigma)
}

func (d *DiscreteGaussian) Equals(other Distribution) bool {
	if other == d {
		return true
	}
	if otherGaus, isGaus := other.(*DiscreteGaussian); isGaus {
		return *d == *otherGaus
	}
	return false
}

func (d *DiscreteGaussian) MarshalJSON() ([]byte, error) {
	return json.Marshal(map[string]interface{}{
		"Type":  discreteGaussian.String(),
		"Sigma": d.Sigma,
		"Bound": d.Bound,
	})
}

// NoiseBound returns floor(StandardDeviation * Bound)
func (d *DiscreteGaussian) NoiseBound() uint64 {
	return uint64(float64(d.Sigma) * float64(d.Bound)) // TODO: is bound really given as a factor of sigma ?
}

func (d *DiscreteGaussian) CopyNew() Distribution {
	return &DiscreteGaussian{d.Sigma, d.Bound}
}

func (d *DiscreteGaussian) MarshalBinarySize() int {
	return 16
}

func (d *DiscreteGaussian) Encode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("data stream is too small: should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}

	binary.LittleEndian.PutUint64(data, math.Float64bits(float64(d.Sigma)))
	binary.LittleEndian.PutUint64(data[8:], uint64(d.Bound))

	return 16, nil
}

func (d *DiscreteGaussian) Decode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("data length should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	d.Sigma = StandardDeviation(math.Float64frombits(binary.LittleEndian.Uint64(data[0:])))
	d.Bound = int(binary.LittleEndian.Uint64(data[8:]))
	return 16, nil
}

// Ternary is a distribution with coefficient uniformly distributed
// in [-1, 0, 1] with probability [(1-P)/2, P, (1-P)/2].
type Ternary struct {
	P float64
	H int
}

func NewTernary(distDef map[string]interface{}) (*Ternary, error) {
	_, hasP := distDef["P"]
	_, hasH := distDef["H"]
	var p float64
	var h int
	var err error
	switch {
	case !hasH && hasP:
		p, err = getFloatFromMap(distDef, "P")
	case hasH && !hasP:
		h, err = getIntFromMap(distDef, "H")
	default:
		err = fmt.Errorf("exactly one of the field P or H need to be set")
	}
	if err != nil {
		return nil, err
	}
	return &Ternary{P: p, H: h}, nil
}

func (d *Ternary) Type() Type {
	return ternary
}

func (d *Ternary) Equals(other Distribution) bool {
	if other == d {
		return true
	}
	if otherTern, isTern := other.(*Ternary); isTern {
		return *d == *otherTern
	}
	return false
}

func (d *Ternary) MarshalJSON() ([]byte, error) {
	return json.Marshal(map[string]interface{}{
		"Type": ternary.String(),
		"P":    d.P,
	})
}

func (d *Ternary) CopyNew() Distribution {
	return &Ternary{d.P, d.H}
}

func (d *Ternary) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return StandardDeviation(math.Sqrt(1 - d.P))
}

func (d *Ternary) MarshalBinarySize() int {
	return 16
}

func (d *Ternary) Encode(data []byte) (ptr int, err error) { // TODO: seems not tested for H
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("data stream is too small: should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	binary.LittleEndian.PutUint64(data, math.Float64bits(d.P))
	binary.LittleEndian.PutUint64(data[8:], uint64(d.H))
	return 16, nil
}

func (d *Ternary) Decode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("invalid data stream: length should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	d.P = math.Float64frombits(binary.LittleEndian.Uint64(data))
	d.H = int(binary.LittleEndian.Uint64(data[8:]))
	return 16, nil

}

// Uniform is a distribution with coefficients uniformly distributed in the given ring.
type Uniform struct{}

func NewUniform(_ map[string]interface{}) (*Uniform, error) {
	return &Uniform{}, nil
}

func (d *Uniform) Type() Type {
	return uniform
}

func (d *Uniform) Equals(other Distribution) bool {
	if other == d {
		return true
	}
	if otherUni, isUni := other.(*Uniform); isUni {
		return *d == *otherUni
	}
	return false
}

func (d *Uniform) MarshalJSON() ([]byte, error) {
	return json.Marshal(map[string]interface{}{
		"Type": uniform.String(),
	})
}

// func (d *Uniform) NewSampler(prng utils.PRNG, baseRing *Ring, montgomery bool) Sampler {
// 	return NewSampler(prng, baseRing, d, montgomery)
// }

func (d *Uniform) CopyNew() Distribution {
	return &Uniform{}
}

func (d *Uniform) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return StandardDeviation(math.Exp2(LogQP) / math.Sqrt(12.0))
}

func (d *Uniform) MarshalBinarySize() int {
	return 0
}

func (d *Uniform) Encode(data []byte) (ptr int, err error) {
	return 0, nil
}

func (d *Uniform) Decode(data []byte) (ptr int, err error) {
	return
}

func getFloatFromMap(distDef map[string]interface{}, key string) (float64, error) {
	val, hasVal := distDef[key]
	if !hasVal {
		return 0, fmt.Errorf("map specifies no value for %s", key)
	}
	f, isFloat := val.(float64)
	if !isFloat {
		return 0, fmt.Errorf("value for key %s in map should be of type float", key)
	}
	return f, nil
}

func getIntFromMap(distDef map[string]interface{}, key string) (int, error) {
	val, hasVal := distDef[key]
	if !hasVal {
		return 0, fmt.Errorf("map specifies no value for %s", key)
	}
	f, isNumeric := val.(float64)
	if !isNumeric && f == float64(int(f)) {
		return 0, fmt.Errorf("value for key %s in map should be an integer", key)
	}
	return int(f), nil
}
