package ring

import (
	"testing"
)

func Test_Float128_DivMult(t *testing.T) {
	var x, y uint64
	var xFloat, yFloat Float128

	x = 0xb80b8d5351c4d81b
	y = 0xd3cd9f41f6606a7d

	xFloat = Float128SetUint64(x)
	yFloat = Float128SetUint64(y)
	xFloat = Float128Div(xFloat, yFloat)
	xFloat = Float128Mul(xFloat, yFloat)

	if x != Float128ToUint64(xFloat) {
		t.Errorf("DivMult")
	}
}

func Benchmark_Float128_Add(b *testing.B) {

	var x, y Float128
	x = Float128SetUint64(1)
	y = Float128SetUint64(3)

	for i := 0; i < b.N; i++ {
		Float128Add(x, y)
	}
}

func Benchmark_Float128_Sub(b *testing.B) {

	var x, y Float128
	x = Float128SetUint64(1)
	y = Float128SetUint64(3)

	for i := 0; i < b.N; i++ {
		Float128Sub(x, y)
	}
}

func Benchmark_Float128_Mul(b *testing.B) {

	var x, y Float128
	x = Float128SetUint64(1)
	y = Float128SetUint64(3)

	for i := 0; i < b.N; i++ {
		Float128Mul(x, y)
	}
}

func Benchmark_Float128_Div(b *testing.B) {

	var x, y Float128
	x = Float128SetUint64(1)
	y = Float128SetUint64(3)

	for i := 0; i < b.N; i++ {
		Float128Div(x, y)
	}
}
