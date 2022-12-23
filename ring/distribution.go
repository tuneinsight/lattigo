package ring

import (
	"encoding/binary"
	"encoding/json"
	"fmt"
	"math"
)

type DistributionType uint8

const (
	Uniform DistributionType = iota + 1
	Ternary
	Gaussian
)

var distributionTypeToString = [5]string{"Undefined", "Uniform", "Ternary", "Gaussian"}

var distributionTypeFromString = map[string]DistributionType{
	"Undefined": 0, "Uniform": Uniform, "Ternary": Ternary, "Gaussian": Gaussian,
}

func (t DistributionType) String() string {
	if int(t) >= len(distributionTypeToString) {
		return "Unknown"
	}
	return distributionTypeToString[int(t)]
}

// Distribution is a interface for distributions
type Distribution interface {
	Type() DistributionType
	StandardDeviation(LogN int, LogQP float64) StandardDeviation // TODO: properly define
	Equals(Distribution) bool
	CopyNew() Distribution

	MarshalBinarySize() int
	Encode(data []byte) (ptr int, err error)
	Decode(data []byte) (ptr int, err error)
}

func NewDistributionFromMap(distDef map[string]interface{}) (Distribution, error) {
	distTypeVal, specified := distDef["Type"]
	if !specified {
		return nil, fmt.Errorf("map specifies no distribution type")
	}
	distTypeStr, isString := distTypeVal.(string)
	if !isString {
		return nil, fmt.Errorf("value for key Type of map should be of type string")
	}
	distType, exists := distributionTypeFromString[distTypeStr]
	if !exists {
		return nil, fmt.Errorf("distribution type %s does not exist", distTypeStr)
	}
	switch distType {
	case Uniform:
		return NewUniformDistributionFromMap(distDef)
	case Ternary:
		return NewTernaryUniformDistribution(distDef)
	case Gaussian:
		return NewDiscreteGaussianDistribution(distDef)
	default:
		return nil, fmt.Errorf("invalid distribution type")
	}
}

func EncodeDistribution(X Distribution, data []byte) (ptr int, err error) {
	if len(data) == 1+X.MarshalBinarySize() {
		return 0, fmt.Errorf("buffer is too small for encoding distribution (size %d instead of %d)", len(data), 1+X.MarshalBinarySize())
	}
	data[0] = byte(X.Type())
	ptr, err = X.Encode(data[1:])

	return ptr + 1, err
}

func DecodeDistribution(data []byte) (ptr int, X Distribution, err error) {
	if len(data) == 0 {
		return 0, nil, fmt.Errorf("data should have length >= 1")
	}
	switch DistributionType(data[0]) {
	case Uniform:
		X = &UniformDistribution{}
	case Ternary:
		X = &TernaryDistribution{}
	case Gaussian:
		X = &DiscreteGaussianDistribution{}
	default:
		return 0, nil, fmt.Errorf("invalid distribution type: %s", DistributionType(data[0]))
	}

	ptr, err = X.Decode(data[1:])

	return ptr + 1, X, err
}

// StandardDeviation is a float64 type storing
// a value representing a standard deviation
type StandardDeviation float64

// DiscreteGaussianDistribution is a discrete Gaussian distribution
// with a given standard deviation and a bound
// in number of standard deviations.
type DiscreteGaussianDistribution struct {
	Sigma StandardDeviation
	Bound int
}

func NewDiscreteGaussianDistribution(distDef map[string]interface{}) (d *DiscreteGaussianDistribution, err error) {
	sigma, errSigma := getFloatFromMap(distDef, "Sigma")
	if errSigma != nil {
		return nil, err
	}
	bound, errBound := getIntFromMap(distDef, "Bound")
	if errBound != nil {
		return nil, err
	}
	return &DiscreteGaussianDistribution{Sigma: StandardDeviation(sigma), Bound: bound}, nil
}

func (d *DiscreteGaussianDistribution) Type() DistributionType {
	return Gaussian
}

func (d *DiscreteGaussianDistribution) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return StandardDeviation(d.Sigma)
}

func (d *DiscreteGaussianDistribution) Equals(other Distribution) bool {
	if other == d {
		return true
	}
	if otherGaus, isGaus := other.(*DiscreteGaussianDistribution); isGaus {
		return *d == *otherGaus
	}
	return false
}

func (d *DiscreteGaussianDistribution) MarshalJSON() ([]byte, error) {
	return json.Marshal(map[string]interface{}{
		"Type":  Gaussian.String(),
		"Sigma": d.Sigma,
		"Bound": d.Bound,
	})
}

// NoiseBound returns floor(StandardDeviation * Bound)
func (d *DiscreteGaussianDistribution) NoiseBound() uint64 {
	return uint64(float64(d.Sigma) * float64(d.Bound)) // TODO: is bound really given as a factor of sigma ?
}

func (d *DiscreteGaussianDistribution) CopyNew() Distribution {
	return &DiscreteGaussianDistribution{d.Sigma, d.Bound}
}

func (d *DiscreteGaussianDistribution) MarshalBinarySize() int {
	return 16
}

func (d *DiscreteGaussianDistribution) Encode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("data stream is too small: should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}

	binary.LittleEndian.PutUint64(data, math.Float64bits(float64(d.Sigma)))
	binary.LittleEndian.PutUint64(data[8:], uint64(d.Bound))

	return 16, nil
}

func (d *DiscreteGaussianDistribution) Decode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("data length should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	d.Sigma = StandardDeviation(math.Float64frombits(binary.LittleEndian.Uint64(data[0:])))
	d.Bound = int(binary.LittleEndian.Uint64(data[8:]))
	return 16, nil
}

// TernaryDistribution is a distribution with coefficient uniformly distributed
// in [-1, 0, 1] with probability [(1-P)/2, P, (1-P)/2].
type TernaryDistribution struct {
	P float64
	H int
}

func NewTernaryUniformDistribution(distDef map[string]interface{}) (*TernaryDistribution, error) {
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
	return &TernaryDistribution{P: p, H: h}, nil
}

func (d *TernaryDistribution) Type() DistributionType {
	return Ternary
}

func (d *TernaryDistribution) Equals(other Distribution) bool {
	if other == d {
		return true
	}
	if otherTern, isTern := other.(*TernaryDistribution); isTern {
		return *d == *otherTern
	}
	return false
}

func (d *TernaryDistribution) MarshalJSON() ([]byte, error) {
	return json.Marshal(map[string]interface{}{
		"Type": Ternary.String(),
		"P":    d.P,
	})
}

func (d *TernaryDistribution) CopyNew() Distribution {
	return &TernaryDistribution{d.P, d.H}
}

func (d *TernaryDistribution) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return StandardDeviation(math.Sqrt(1 - d.P))
}

func (d *TernaryDistribution) MarshalBinarySize() int {
	return 16
}

func (d *TernaryDistribution) Encode(data []byte) (ptr int, err error) { // TODO: seems not tested for H
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("data stream is too small: should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	binary.LittleEndian.PutUint64(data, math.Float64bits(d.P))
	binary.LittleEndian.PutUint64(data[8:], uint64(d.H))
	return 16, nil
}

func (d *TernaryDistribution) Decode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("invalid data stream: length should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	d.P = math.Float64frombits(binary.LittleEndian.Uint64(data))
	d.H = int(binary.LittleEndian.Uint64(data[8:]))
	return 16, nil

}

// UniformDistribution is a distribution with coefficients uniformly distributed in the given ring.
type UniformDistribution struct{}

func NewUniformDistributionFromMap(_ map[string]interface{}) (*UniformDistribution, error) {
	return &UniformDistribution{}, nil
}

func (d *UniformDistribution) Type() DistributionType {
	return Uniform
}

func (d *UniformDistribution) Equals(other Distribution) bool {
	if other == d {
		return true
	}
	if otherUni, isUni := other.(*UniformDistribution); isUni {
		return *d == *otherUni
	}
	return false
}

func (d *UniformDistribution) MarshalJSON() ([]byte, error) {
	return json.Marshal(map[string]interface{}{
		"Type": Uniform.String(),
	})
}

// func (d *Uniform) NewSampler(prng utils.PRNG, baseRing *Ring, montgomery bool) Sampler {
// 	return NewSampler(prng, baseRing, d, montgomery)
// }

func (d *UniformDistribution) CopyNew() Distribution {
	return &UniformDistribution{}
}

func (d *UniformDistribution) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return StandardDeviation(math.Exp2(LogQP) / math.Sqrt(12.0))
}

func (d *UniformDistribution) MarshalBinarySize() int {
	return 0
}

func (d *UniformDistribution) Encode(data []byte) (ptr int, err error) {
	return 0, nil
}

func (d *UniformDistribution) Decode(data []byte) (ptr int, err error) {
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
