package utils

// PointyInt creates a new int variable and returns its pointer.
func PointyInt(x int) *int {
	return &x
}

// PointyUint64 creates a new uint64 variable and returns its pointer.
func PointyUint64(x uint64) *uint64 {
	return &x
}
