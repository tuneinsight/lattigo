package ring

type Dimensions struct {
	Rows, Cols int
}

// EvalPolyModP evaluates y = sum poly[i] * x^{i} mod p.
func EvalPolyModP(x uint64, poly []uint64, p uint64) (y uint64) {
	brc := GenBRedConstant(p)
	y = poly[len(poly)-1]
	for i := len(poly) - 2; i >= 0; i-- {
		y = BRed(y, x, p, brc)
		y = CRed(y+poly[i], p)
	}

	return
}

// Min returns the minimum between to int
func Min(x, y int) int {
	if x > y {
		return y
	}

	return x
}

// ModExp performs the modular exponentiation x^e mod p,
// x and p are required to be at most 64 bits to avoid an overflow.
func ModExp(x, e, p uint64) (result uint64) {
	brc := GenBRedConstant(p)
	result = 1
	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result = BRed(result, x, p, brc)
		}
		x = BRed(x, x, p, brc)
	}
	return result
}

// ModExpPow2 performs the modular exponentiation x^e mod p, where p is a power of two,
// x and p are required to be at most 64 bits to avoid an overflow.
func ModExpPow2(x, e, p uint64) (result uint64) {

	result = 1
	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result *= x
		}
		x *= x
	}
	return result & (p - 1)
}

// ModexpMontgomery performs the modular exponentiation x^e mod p,
// where x is in Montgomery form, and returns x^e in Montgomery form.
func ModexpMontgomery(x uint64, e int, q, mredconstant uint64, bredconstant [2]uint64) (result uint64) {

	result = MForm(1, q, bredconstant)

	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result = MRed(result, x, q, mredconstant)
		}
		x = MRed(x, x, q, mredconstant)
	}
	return result
}
