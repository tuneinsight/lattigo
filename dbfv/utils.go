package dbfv

// equalslice compares two slices of uint64 values, and return true if they are equal, else false.
func equalslice(a, b []uint64) bool {

	if len(a) != len(b) {
		return false
	}

	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
