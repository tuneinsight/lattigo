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
