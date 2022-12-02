package utils

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestAllDistinct(t *testing.T) {
	require.True(t, AllDistinct([]uint64{}))
	require.True(t, AllDistinct([]uint64{1}))
	require.True(t, AllDistinct([]uint64{1, 2, 3}))
	require.False(t, AllDistinct([]uint64{1, 1}))
	require.False(t, AllDistinct([]uint64{1, 2, 3, 4, 5, 5}))
}

func TestRotateUint64(t *testing.T) {
	s := []uint64{0, 1, 2, 3, 4, 5, 6, 7}
	sout := make([]uint64, len(s))

	RotateUint64SliceAllocFree(s, 3, sout)
	require.Equal(t, []uint64{3, 4, 5, 6, 7, 0, 1, 2}, sout)
	require.Equal(t, []uint64{0, 1, 2, 3, 4, 5, 6, 7}, s, "should not modify input slice")

	RotateUint64SliceAllocFree(s, 0, sout)
	require.Equal(t, []uint64{0, 1, 2, 3, 4, 5, 6, 7}, sout)

	RotateUint64SliceAllocFree(s, -2, sout)
	require.Equal(t, []uint64{6, 7, 0, 1, 2, 3, 4, 5}, sout)

	RotateUint64SliceAllocFree(s, 9, sout)
	require.Equal(t, []uint64{1, 2, 3, 4, 5, 6, 7, 0}, sout)

	RotateUint64SliceAllocFree(s, -11, sout)
	require.Equal(t, []uint64{5, 6, 7, 0, 1, 2, 3, 4}, sout)

	RotateUint64SliceAllocFree(s, 0, s)
	require.Equal(t, []uint64{0, 1, 2, 3, 4, 5, 6, 7}, s)

	RotateUint64SliceAllocFree(s, 1, s)
	require.Equal(t, []uint64{1, 2, 3, 4, 5, 6, 7, 0}, s)

	RotateUint64SliceAllocFree(s, -2, s)
	require.Equal(t, []uint64{7, 0, 1, 2, 3, 4, 5, 6}, s)
}

func TestGetFactors(t *testing.T) {

	m := new(big.Int).SetUint64(35184372088631)

	t.Run("ECM", func(t *testing.T) {

		factor := GetFactorECM(m)

		if factor.Cmp(new(big.Int).SetUint64(6292343)) != 0 && factor.Cmp(new(big.Int).SetUint64(5591617)) != 0 {
			t.Fail()
		}
	})

	t.Run("PollardRho", func(t *testing.T) {
		factor := GetFactorPollardRho(m)

		if factor.Cmp(new(big.Int).SetUint64(6292343)) != 0 && factor.Cmp(new(big.Int).SetUint64(5591617)) != 0 {
			t.Fail()
		}
	})
}
