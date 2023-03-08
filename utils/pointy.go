package utils

// PointyInt creates a new int variable and returns its pointer.
func PointyInt(x int) *int {
	return &x
}
