package ring

// RNSScalar represents a scalar value in the Ring (i.e., a degree-0 polynomial) in RNS form.
type RNSScalar []uint64

// NewRNSScalar creates a new Scalar value.
func (r *Ring) NewRNSScalar() RNSScalar {
	return make(RNSScalar, r.NbModuli())
}

// NewRNSScalarFromUInt64 creates a new Scalar initialized with value v.
func (r *Ring) NewRNSScalarFromUInt64(v uint64) RNSScalar {
	s := make(RNSScalar, r.NbModuli())
	for i, table := range r.Tables {
		s[i] = v % table.Modulus
	}
	return s
}

// SubRNSScalar subtracts s2 to s1 and stores the result in sout.
func (r *Ring) SubRNSScalar(s1, s2, sout RNSScalar) {
	for i, table := range r.Tables {
		if s2[i] > s1[i] {
			sout[i] = s1[i] + table.Modulus - s2[i]
		} else {
			sout[i] = s1[i] - s2[i]
		}
	}
}

// MulRNSScalar multiplies s1 and s2 and stores the result in sout.
// Multiplication is oprated with Montgomery.
func (r *Ring) MulRNSScalar(s1, s2, sout RNSScalar) {
	for i, table := range r.Tables {
		sout[i] = MRedConstant(s1[i], s2[i], table.Modulus, table.MRedParams)
	}
}

// Inverse computes the modular inverse of a scalar a expressed in a CRT decomposition.
// The inversion is done in-place and assumes that a is in Montgomery form.
func (r *Ring) Inverse(a RNSScalar) {
	for i, table := range r.Tables {
		a[i] = ModexpMontgomery(a[i], int(table.Modulus-2), table.Modulus, table.MRedParams, table.BRedParams)
	}
}
