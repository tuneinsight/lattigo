package ring

import (
	"testing"
)

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
