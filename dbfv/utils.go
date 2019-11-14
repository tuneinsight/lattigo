package dbfv

// equalslice compares two slices of uint64 values, and return true if they are equal, else false.
func equalslice(a, b []uint64) bool {

	if len(a) != len(b) {
		return false
	}

	v := true
	for i := range a {
		v = v && (a[i] == b[i])
	}
	return v
}
