package utils

import (
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

	RotateSliceAllocFree(s, 3, sout)
	require.Equal(t, []uint64{3, 4, 5, 6, 7, 0, 1, 2}, sout)
	require.Equal(t, []uint64{0, 1, 2, 3, 4, 5, 6, 7}, s, "should not modify input slice")

	RotateSliceAllocFree(s, 0, sout)
	require.Equal(t, []uint64{0, 1, 2, 3, 4, 5, 6, 7}, sout)

	RotateSliceAllocFree(s, -2, sout)
	require.Equal(t, []uint64{6, 7, 0, 1, 2, 3, 4, 5}, sout)

	RotateSliceAllocFree(s, 9, sout)
	require.Equal(t, []uint64{1, 2, 3, 4, 5, 6, 7, 0}, sout)

	RotateSliceAllocFree(s, -11, sout)
	require.Equal(t, []uint64{5, 6, 7, 0, 1, 2, 3, 4}, sout)

	RotateSliceAllocFree(s, 0, s)
	require.Equal(t, []uint64{0, 1, 2, 3, 4, 5, 6, 7}, s)

	RotateSliceAllocFree(s, 1, s)
	require.Equal(t, []uint64{1, 2, 3, 4, 5, 6, 7, 0}, s)

	RotateSliceAllocFree(s, -2, s)
	require.Equal(t, []uint64{7, 0, 1, 2, 3, 4, 5, 6}, s)
}
