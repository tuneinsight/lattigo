package ring

import (
	"math/bits"
)

// EvalPolyModP evaluates y = sum poly[i] * x^{i} mod p.
func EvalPolyModP(x uint64, poly []uint64, p uint64) (y uint64) {
	bredParams := BRedParams(p)
	y = poly[len(poly)-1]
	for i := len(poly) - 2; i >= 0; i-- {
		y = BRed(y, x, p, bredParams)
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

// PowerOf2 returns (x*2^n)%q where x is in Montgomery form
func PowerOf2(x uint64, n int, q, qInv uint64) (r uint64) {
	ahi, alo := x>>(64-n), x<<n
	R := alo * qInv
	H, _ := bits.Mul64(R, q)
	r = ahi - H + q
	if r >= q {
		r -= q
	}
	return
}

// ModExp performs the modular exponentiation x^e mod p,
// x and p are required to be at most 64 bits to avoid an overflow.
func ModExp(x, e, p uint64) (result uint64) {
	params := BRedParams(p)
	result = 1
	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result = BRed(result, x, p, params)
		}
		x = BRed(x, x, p, params)
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
func ModexpMontgomery(x uint64, e int, q, qInv uint64, bredParams []uint64) (result uint64) {

	result = MForm(1, q, bredParams)

	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result = MRed(result, x, q, qInv)
		}
		x = MRed(x, x, q, qInv)
	}
	return result
}

// primitiveRoot computes one primitive root (the smallest) of for the given prime q
func primitiveRoot(q uint64) (g uint64) {
	var tmp uint64

	notFoundPrimitiveRoot := true

	factors := GetFactors(q - 1) //Factor q-1, might be slow

	g = 2

	for notFoundPrimitiveRoot {
		g++
		for _, factor := range factors {
			tmp = (q - 1) / factor
			// if for any factor of q-1, g^(q-1)/factor = 1 mod q, g is not a primitive root
			if ModExp(g, tmp, q) == 1 {
				notFoundPrimitiveRoot = true
				break
			}
			notFoundPrimitiveRoot = false
		}
	}
	return
}
