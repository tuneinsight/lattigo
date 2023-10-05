package structs

import (
	"testing"

	"github.com/google/go-cmp/cmp"
	"github.com/stretchr/testify/require"
	"golang.org/x/exp/constraints"
)

func TestStructs(t *testing.T) {
	t.Run("Vector/W64/Serialization&Equatable", func(t *testing.T) {
		testVector[uint64](t)
	})

	t.Run("Vector/W32/Serialization&Equatable", func(t *testing.T) {
		testVector[uint32](t)
	})

	t.Run("Vector/W16/Serialization&Equatable", func(t *testing.T) {
		testVector[uint16](t)
	})

	t.Run("Vector/W8/Serialization&Equatable", func(t *testing.T) {
		testVector[uint8](t)
	})

	t.Run("Matrix/W64/Serialization&Equatable", func(t *testing.T) {
		testMatrix[float64](t)
	})

	t.Run("Matrix/W32/Serialization&Equatable", func(t *testing.T) {
		testMatrix[float64](t)
	})

	t.Run("Matrix/W16/Serialization&Equatable", func(t *testing.T) {
		testMatrix[float64](t)
	})

	t.Run("Matrix/W8/Serialization&Equatable", func(t *testing.T) {
		testMatrix[float64](t)
	})
}

func testVector[T constraints.Float | constraints.Integer](t *testing.T) {
	v := Vector[T](make([]T, 64))
	for i := range v {
		v[i] = T(i)
	}
	data, err := v.MarshalBinary()
	require.NoError(t, err)
	vNew := Vector[T]{}
	require.NoError(t, vNew.UnmarshalBinary(data))
	require.True(t, cmp.Equal(v, vNew)) // also tests Equatable
}

func testMatrix[T constraints.Float | constraints.Integer](t *testing.T) {
	v := Matrix[T](make([][]T, 64))
	for i := range v {
		vi := make([]T, 64)
		for j := range vi {
			vi[j] = T(i & j)
		}

		v[i] = vi
	}

	data, err := v.MarshalBinary()
	require.NoError(t, err)
	vNew := Matrix[T]{}
	require.NoError(t, vNew.UnmarshalBinary(data))
	require.True(t, cmp.Equal(v, vNew)) // also tests Equatable
}
