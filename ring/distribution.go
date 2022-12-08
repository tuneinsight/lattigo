package ring

import (
	"encoding/binary"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// Distribution is a interface for distributions
type Distribution interface {
	CopyNew() Distribution
	StandardDeviation(LogN int, LogQP float64) StandardDeviation
	MarshalBinarySize() int
	Encode(data []byte) (ptr int, err error)
	Decode(data []byte) (ptr int, err error)
	MarshalBinary() (data []byte, err error)
	UnmarshalBinary(data []byte) (err error)
	NewSampler(prng utils.PRNG, baseRing *Ring, montgomery bool) Sampler
}

func EncodeDistribution(X Distribution, data []byte) (ptr int, err error) {
	switch X.(type) {
	case *DiscreteGaussian:
		data[0] = 0
	case *UniformTernary:
		data[0] = 1
	case *SparseTernary:
		data[0] = 2
	case *Uniform:
		data[0] = 3
	}

	ptr, err = X.Encode(data[1:])

	return ptr + 1, err
}

func MarshalDistribution(X Distribution) (data []byte, err error) {

	b := utils.NewBuffer([]byte{})

	switch X.(type) {
	case *DiscreteGaussian:
		b.WriteUint8(0)
	case *UniformTernary:
		b.WriteUint8(1)
	case *SparseTernary:
		b.WriteUint8(2)
	case *Uniform:
		b.WriteUint8(3)
	}

	var Xdata []byte
	if Xdata, err = X.MarshalBinary(); err != nil {
		return
	}

	b.Write(Xdata)

	return b.Bytes(), nil
}

func DecodeDistribution(data []byte) (ptr int, X Distribution, err error) {
	switch data[0] {
	case 0:
		X = &DiscreteGaussian{}
	case 1:
		X = &UniformTernary{}
	case 2:
		X = &SparseTernary{}
	case 3:
		X = &Uniform{}
	}

	ptr, err = X.Decode(data[1:])

	return ptr + 1, X, err
}

func UnmarshalDistribution(data []byte) (X Distribution, err error) {

	switch data[0] {
	case 0:
		X = &DiscreteGaussian{}
	case 1:
		X = &UniformTernary{}
	case 2:
		X = &SparseTernary{}
	case 3:
		X = &Uniform{}
	}

	return X, X.UnmarshalBinary(data[1:])
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

func (d *DiscreteGaussian) NewSampler(prng utils.PRNG, baseRing *Ring, montgomery bool) Sampler {
	return NewSampler(prng, baseRing, d, montgomery)
}

// NoiseBound returns floor(StandardDeviation * Bound)
func (d *DiscreteGaussian) NoiseBound() uint64 {
	return uint64(float64(d.Sigma) * float64(d.Bound))
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

	binary.LittleEndian.PutUint64(data[0:], math.Float64bits(float64(d.Sigma)))
	binary.LittleEndian.PutUint64(data[8:], uint64(d.Bound))

	return 16, nil
}

func (d *DiscreteGaussian) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 16)
	_, err = d.Encode(data)
	return
}

func (d *DiscreteGaussian) Decode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("invalid data stream: length should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	d.Sigma = StandardDeviation(math.Float64frombits(binary.LittleEndian.Uint64(data[0:])))
	d.Bound = int(binary.LittleEndian.Uint64(data[8:]))
	return 16, nil
}

func (d *DiscreteGaussian) UnmarshalBinary(data []byte) (err error) {
	var ptr int
	ptr, err = d.Decode(data)

	if len(data) > ptr {
		return fmt.Errorf("remaining unparsed data")
	}

	return
}

func (d *DiscreteGaussian) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return d.Sigma
}

// UniformTernary is a distribution with coefficient uniformly distributed
// in [-1, 0, 1] with probability [(1-P)/2, P, (1-P)/2].
type UniformTernary struct {
	P float64
}

func (d *UniformTernary) NewSampler(prng utils.PRNG, baseRing *Ring, montgomery bool) Sampler {
	return NewSampler(prng, baseRing, d, montgomery)
}

func (d *UniformTernary) CopyNew() Distribution {
	return &UniformTernary{d.P}
}

func (d *UniformTernary) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return StandardDeviation(math.Sqrt(1 - d.P))
}

func (d *UniformTernary) MarshalBinarySize() int {
	return 8
}

func (d *UniformTernary) Encode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("data stream is too small: should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	binary.LittleEndian.PutUint64(data, math.Float64bits(d.P))
	return 8, nil
}

func (d *UniformTernary) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 8)
	_, err = d.Encode(data)
	return
}

func (d *UniformTernary) Decode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("invalid data stream: length should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}

	d.P = math.Float64frombits(binary.LittleEndian.Uint64(data))
	return 8, nil

}

func (d *UniformTernary) UnmarshalBinary(data []byte) (err error) {
	var ptr int
	ptr, err = d.Decode(data)

	if len(data) > ptr {
		return fmt.Errorf("remaining unparsed data")
	}

	return
}

// SparseTernary is a distribution with exactly `HammingWeight`coefficients  uniformly distributed in [-1, 1]
type SparseTernary struct {
	HammingWeight int
}

func (d *SparseTernary) NewSampler(prng utils.PRNG, baseRing *Ring, montgomery bool) Sampler {
	return NewSampler(prng, baseRing, d, montgomery)
}

func (d *SparseTernary) CopyNew() Distribution {
	return &SparseTernary{d.HammingWeight}
}

func (d *SparseTernary) StandardDeviation(LogN int, LogQP float64) StandardDeviation {
	return StandardDeviation(math.Sqrt(float64(d.HammingWeight) / math.Exp2(float64(LogN))))
}

func (d *SparseTernary) MarshalBinarySize() int {
	return 4
}

func (d *SparseTernary) Encode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("data stream is too small: should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	binary.LittleEndian.PutUint32(data, uint32(d.HammingWeight))
	return 4, nil
}

func (d *SparseTernary) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 4)
	_, err = d.Encode(data)
	return
}

func (d *SparseTernary) Decode(data []byte) (ptr int, err error) {
	if len(data) < d.MarshalBinarySize() {
		return ptr, fmt.Errorf("invalid data stream: length should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}
	d.HammingWeight = int(binary.LittleEndian.Uint32(data))
	return 4, nil
}

func (d *SparseTernary) UnmarshalBinary(data []byte) (err error) {
	if len(data) < d.MarshalBinarySize() {
		return fmt.Errorf("invalid data stream: length should be at least %d but is %d", d.MarshalBinarySize(), len(data))
	}

	var ptr int
	ptr, err = d.Decode(data)

	if len(data) > ptr {
		return fmt.Errorf("remaining unparsed data")
	}

	return
}

// Uniform is a distribution with coefficients uniformly distributed in the given ring.
type Uniform struct{}

func (d *Uniform) NewSampler(prng utils.PRNG, baseRing *Ring, montgomery bool) Sampler {
	return NewSampler(prng, baseRing, d, montgomery)
}

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
	return
}

func (d *Uniform) MarshalBinary() (data []byte, err error) {
	return
}

func (d *Uniform) Decode(data []byte) (ptr int, err error) {
	return
}

func (d *Uniform) UnmarshalBinary(data []byte) (err error) {
	return
}
