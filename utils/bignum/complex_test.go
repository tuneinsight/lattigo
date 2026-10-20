package bignum

import (
	"testing"

	"github.com/stretchr/testify/require"
)

func TestComplexMultiplierQuo(t *testing.T) {
	const prec = 128

	tests := []struct {
		name string
		a    complex128
		b    complex128
		want complex128
	}{
		{
			name: "real divided by complex",
			a:    complex(1, 0),
			b:    complex(1, 1),
			want: complex(0.5, -0.5),
		},
		{
			name: "real divided by complex with negative imaginary part",
			a:    complex(2, 0),
			b:    complex(1, -1),
			want: complex(1, 1),
		},
		{
			name: "complex divided by real",
			a:    complex(2, 3),
			b:    complex(2, 0),
			want: complex(1, 1.5),
		},
		{
			name: "complex divided by complex",
			a:    complex(2, 3),
			b:    complex(1, -1),
			want: complex(-0.5, 2.5),
		},
		{
			name: "real divided by real",
			a:    complex(6, 0),
			b:    complex(2, 0),
			want: complex(3, 0),
		},
	}

	eval := NewComplexMultiplier()

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			a := ToComplex(test.a, prec)
			b := ToComplex(test.b, prec)
			got := ToComplex(0, prec)
			want := ToComplex(test.want, prec)

			eval.Quo(a, b, got)

			require.Zero(t, got[0].Cmp(want[0]), "unexpected real part")
			require.Zero(t, got[1].Cmp(want[1]), "unexpected imaginary part")
		})
	}
}
