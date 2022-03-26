package lwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
)

// m0 : [0, 1/4]
// m1 : [0, 1/4]
// | 0 0 -> 1 (0/8) -> 2/8
// | 0 1 -> 1 (2/8) -> 2/8
// | 1 0 -> 1 (2/8) -> 2/8
// | 1 1 -> 0 (4/8) -> 0/8
func nandGate(x float64) float64 {
	if x > -1/8.0 && x < 3/8.0 {
		return 2 / 8.0
	}

	return 0
}

// m0 : [0, 1/4]
// m1 : [0, 1/4]
// | 0 0 -> 0 (0/8) -> 0/8
// | 0 1 -> 0 (2/8) -> 0/8
// | 1 0 -> 0 (2/8) -> 0/8
// | 1 1 -> 1 (4/8) -> 2/8
func andGate(x float64) float64 {
	if x > -1/8.0 && x < 3/8.0 {
		return 0
	}

	return 1 / 4.0
}

// m0 : [0, 1/4]
// m1 : [0, 1/4]
// | 0 0 -> 0 (0/8) -> 0/8
// | 0 1 -> 1 (2/8) -> 2/8
// | 1 0 -> 1 (2/8) -> 2/8
// | 1 1 -> 0 (4/8) -> 0/8
func xorGate(x float64) float64 {
	if x > 1/8.0 && x < 3/8.0 {
		return 2 / 8.0
	}

	return 0
}

// m0 : [0, 1/4]
// m1 : [0, 1/4]
// | 0 0 -> 1 (0/8) -> 2/8
// | 0 1 -> 0 (2/8) -> 0/8
// | 1 0 -> 0 (2/8) -> 0/8
// | 1 1 -> 1 (4/8) -> 2/8
func nxorGate(x float64) float64 {
	if x > 1/8.0 && x < 3/8.0 {
		return 0
	}

	return 2 / 8.0
}

// m0 : [0, 1/4]
// m1 : [0, 1/4]
// | 0 0 -> 0 (0/8) -> 0/8
// | 0 1 -> 1 (2/8) -> 2/8
// | 1 0 -> 1 (2/8) -> 2/8
// | 1 1 -> 1 (4/8) -> 2/8
func orGate(x float64) float64 {
	if x > 1/8.0 && x < 5/8.0 {
		return 2 / 8.0
	}

	return 0
}

// m0 : [0, 1/4]
// m1 : [0, 1/4]
// | 0 0 -> 1 (0/8) -> 2/8
// | 0 1 -> 0 (2/8) -> 0/8
// | 1 0 -> 0 (2/8) -> 0/8
// | 1 1 -> 0 (4/8) -> 0/8
func norGate(x float64) float64 {
	if x > 1/8.0 && x < 5/8.0 {
		return 0
	}

	return 2 / 8.0
}

// m0 : [0, 1/4]
// m1 : [0, 1/4]
// | 0 -> 1 (0/8) -> 2/8
// | 1 -> 0 (2/8) -> 0/8
func notGate(x float64) float64 {
	if x > 1/8.0 && x < 3/8.0 {
		return 0
	}

	return 2 / 8.0
}

// InitGate generate the test rlwe plaintext for the selected gate.
func InitGate(gate func(x float64) float64, r *ring.Ring) *ring.Poly {
	return InitLUT(gate, float64(r.Modulus[0])/4.0, r, -1, 1)
}
