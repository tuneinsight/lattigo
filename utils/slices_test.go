package utils

import (
	"sort"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestGetSortedKeys(t *testing.T) {
	m := map[int]int{1: 1, 3: 3, 2: 2}
	require.Equal(t, []int{1, 2, 3}, GetSortedKeys(m))
	m = map[int]int{-1: 1, -3: 3, -2: 2}
	require.Equal(t, []int{-3, -2, -1}, GetSortedKeys(m))
}

func TestGetDistincts(t *testing.T) {
	actual := GetDistincts([]int{1, 2})
	expected := []int{1, 2}
	sort.Ints(expected)
	sort.Ints(actual)
	require.Equal(t, expected, actual)

	actual = GetDistincts([]int{1, 2, 3, 1, 2, 3})
	expected = []int{1, 2, 3}
	sort.Ints(expected)
	sort.Ints(actual)
	require.Equal(t, expected, actual)

	actual = GetDistincts([]int{-1, 1, 1, 1})
	expected = []int{-1, 1}
	sort.Ints(expected)
	sort.Ints(actual)
	require.Equal(t, expected, actual)
}

func TestRotateSlice(t *testing.T) {
	actual := RotateSlice([]int{1, 2, 3, 4, 5}, 2)
	expected := []int{3, 4, 5, 1, 2}
	require.Equal(t, expected, actual)

	actual = RotateSlice([]int{1, 2, 3, 4, 5}, -2)
	expected = []int{4, 5, 1, 2, 3}
	require.Equal(t, expected, actual)

	actual = RotateSlice([]int{1, 2, 3, 4, 5}, 0)
	expected = []int{1, 2, 3, 4, 5}
	require.Equal(t, expected, actual)
}

func TestRotateSliceInPlace(t *testing.T) {
	slice := []int{1, 2, 3, 4, 5}
	RotateSliceInPlace(slice, 2)
	expected := []int{3, 4, 5, 1, 2}
	require.Equal(t, expected, slice)

	slice = []int{1, 2, 3, 4, 5}
	RotateSliceInPlace(slice, -2)
	expected = []int{4, 5, 1, 2, 3}
	require.Equal(t, expected, slice)

	slice = []int{1, 2, 3, 4, 5}
	RotateSliceInPlace(slice, 0)
	expected = []int{1, 2, 3, 4, 5}
	require.Equal(t, expected, slice)
}
