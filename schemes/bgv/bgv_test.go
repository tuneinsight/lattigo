package bgv

import (
	"encoding/json"
	"flag"
	"fmt"
	"runtime"
	"slices"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

var flagPrintNoise = flag.Bool("print-noise", false, "print the residual noise")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")

func name(op string, tc *TestContext, lvl int) string {
	return fmt.Sprintf("%s/%s/lvl=%d", op, tc, lvl)
}

func TestBGV(t *testing.T) {

	var err error

	paramsLiterals := testParams

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:] {

		for _, plaintextModulus := range testPlaintextModulus[:] {

			p.PlaintextModulus = plaintextModulus

			tc := NewTestContext(p, false)

			for _, testSet := range []func(tc *TestContext, t *testing.T){
				testParameters,
				testEncoder,
				testEvaluatorBvg,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}
}

func TestBFV(t *testing.T) {

	var err error

	paramsLiterals := testParams

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:] {

		for _, plaintextModulus := range testPlaintextModulus[:] {

			p.PlaintextModulus = plaintextModulus

			tc := NewTestContext(p, true)

			for _, testSet := range []func(tc *TestContext, t *testing.T){
				testEvaluatorBfv,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}
}

func testParameters(tc *TestContext, t *testing.T) {
	t.Run(name("Parameters/Binary", tc, 0), func(t *testing.T) {

		bytes, err := tc.Params.MarshalBinary()
		require.Nil(t, err)
		var p Parameters
		require.Nil(t, p.UnmarshalBinary(bytes))
		require.True(t, tc.Params.Equal(&p))

	})

	t.Run(name("Parameters/JSON", tc, 0), func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(tc.Params)
		require.Nil(t, err)
		require.NotNil(t, data)

		// checks that the Parameters can be unmarshalled without error
		var paramsRec Parameters
		err = json.Unmarshal(data, &paramsRec)
		require.Nil(t, err)
		require.True(t, tc.Params.Equal(&paramsRec))

		// checks that the Parameters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "PlaintextModulus":65537}`, tc.Params.LogN()))
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		require.Nil(t, err)
		require.Equal(t, 2, paramsWithLogModuli.QCount())
		require.Equal(t, 1, paramsWithLogModuli.PCount())
		require.Equal(t, rlwe.DefaultXe, paramsWithLogModuli.Xe()) // Omitting Xe should result in Default being used
		require.Equal(t, rlwe.DefaultXs, paramsWithLogModuli.Xs()) // Omitting Xe should result in Default being used

		// checks that one can provide custom parameters for the secret-key and error distributions
		dataWithCustomSecrets := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "PlaintextModulus":65537, "Xs": {"Type": "Ternary", "H": 192}, "Xe": {"Type": "DiscreteGaussian", "Sigma": 6.6, "Bound": 39.6}}`, tc.Params.LogN()))
		var paramsWithCustomSecrets Parameters
		err = json.Unmarshal(dataWithCustomSecrets, &paramsWithCustomSecrets)
		require.Nil(t, err)
		require.Equal(t, ring.DiscreteGaussian{Sigma: 6.6, Bound: 39.6}, paramsWithCustomSecrets.Xe())
		require.Equal(t, ring.Ternary{H: 192}, paramsWithCustomSecrets.Xs())
	})
}

func testEncoder(tc *TestContext, t *testing.T) {
	testLevel := [2]int{0, tc.Params.MaxLevel()}

	for _, lvl := range testLevel {
		t.Run(name("Encoder/Uint/IsBatched=true", tc, lvl), func(t *testing.T) {
			values, plaintext, _ := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())
			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, plaintext, values, t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Encoder/Int/IsBatched=true", tc, lvl), func(t *testing.T) {

			T := tc.Params.PlaintextModulus()
			THalf := T >> 1
			poly := tc.Sampler.ReadNew()
			coeffs := make([]int64, poly.N())
			for i, c := range poly.Coeffs[0] {
				c %= T
				if c >= THalf {
					coeffs[i] = -int64(T - c)
				} else {
					coeffs[i] = int64(c)
				}
			}

			plaintext := NewPlaintext(tc.Params, lvl)
			tc.Ecd.Encode(coeffs, plaintext)
			have := make([]int64, tc.Params.MaxSlots())
			tc.Ecd.Decode(plaintext, have)
			require.True(t, slices.Equal(coeffs, have))
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Encoder/Uint/IsBatched=false", tc, lvl), func(t *testing.T) {
			T := tc.Params.PlaintextModulus()
			poly := tc.Sampler.ReadNew()
			coeffs := make([]uint64, poly.N())
			for i, c := range poly.Coeffs[0] {
				coeffs[i] = c % T
			}

			plaintext := NewPlaintext(tc.Params, lvl)
			plaintext.IsBatched = false
			require.NoError(t, tc.Ecd.Encode(coeffs, plaintext))
			have := make([]uint64, tc.Params.MaxSlots())
			require.NoError(t, tc.Ecd.Decode(plaintext, have))
			require.True(t, slices.Equal(coeffs, have))
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Encoder/Int/IsBatched=false", tc, lvl), func(t *testing.T) {

			T := tc.Params.PlaintextModulus()
			THalf := T >> 1
			poly := tc.Sampler.ReadNew()
			coeffs := make([]int64, poly.N())
			for i, c := range poly.Coeffs[0] {
				c %= T
				if c >= THalf {
					coeffs[i] = -int64(T - c)
				} else {
					coeffs[i] = int64(c)
				}
			}

			plaintext := NewPlaintext(tc.Params, lvl)
			plaintext.IsBatched = false
			require.NoError(t, tc.Ecd.Encode(coeffs, plaintext))
			have := make([]int64, tc.Params.MaxSlots())
			require.NoError(t, tc.Ecd.Decode(plaintext, have))
			require.True(t, slices.Equal(coeffs, have))
		})
	}
}

func testEvaluatorBvg(tc *TestContext, t *testing.T) {
	testLevel := [2]int{0, tc.Params.MaxLevel()}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Add/Ct/Ct/New", tc, lvl), func(t *testing.T) {
			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			ciphertext2, err := tc.Evl.AddNew(ciphertext0, ciphertext1)
			require.NoError(t, err)
			tc.Params.RingT().Add(p0, p1, p0)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext2, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Add/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			require.NoError(t, tc.Evl.Add(ciphertext0, ciphertext1, ciphertext0))
			tc.Params.RingT().Add(p0, p1, p0)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Add/Ct/Pt/Inplace", tc, lvl), func(t *testing.T) {
			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, plaintext, _ := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(plaintext.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			require.NoError(t, tc.Evl.Add(ciphertext0, plaintext, ciphertext0))
			tc.Params.RingT().Add(p0, p1, p0)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Add/Ct/Scalar/Inplace", tc, lvl), func(t *testing.T) {
			values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

			scalar := tc.Params.PlaintextModulus() >> 1

			p := ring.Poly{Coeffs: [][]uint64{values}}

			require.NoError(t, tc.Evl.Add(ciphertext, scalar, ciphertext))
			tc.Params.RingT().AddScalar(p, scalar, p)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Add/Ct/Vector/Inplace", tc, lvl), func(t *testing.T) {
			values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

			p := ring.Poly{Coeffs: [][]uint64{values}}

			require.NoError(t, tc.Evl.Add(ciphertext, values, ciphertext))
			tc.Params.RingT().Add(p, p, p)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Sub/Ct/Ct/New", tc, lvl), func(t *testing.T) {
			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			ciphertext0, err := tc.Evl.SubNew(ciphertext0, ciphertext1)
			require.NoError(t, err)
			tc.Params.RingT().Sub(p0, p1, p0)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Sub/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			require.NoError(t, tc.Evl.Sub(ciphertext0, ciphertext1, ciphertext0))
			tc.Params.RingT().Sub(p0, p1, p0)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Sub/Ct/Pt/Inplace", tc, lvl), func(t *testing.T) {
			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, plaintext, _ := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(plaintext.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			require.NoError(t, tc.Evl.Sub(ciphertext0, plaintext, ciphertext0))
			tc.Params.RingT().Sub(p0, p1, p0)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Sub/Ct/Scalar/Inplace", tc, lvl), func(t *testing.T) {
			values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

			scalar := tc.Params.PlaintextModulus() >> 1

			p := ring.Poly{Coeffs: [][]uint64{values}}

			require.NoError(t, tc.Evl.Sub(ciphertext, scalar, ciphertext))
			tc.Params.RingT().SubScalar(p, scalar, p)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Sub/Ct/Vector/Inplace", tc, lvl), func(t *testing.T) {
			values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

			p := ring.Poly{Coeffs: [][]uint64{values}}

			require.NoError(t, tc.Evl.Sub(ciphertext, values, ciphertext))
			tc.Params.RingT().Sub(p, p, p)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Mul/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			require.NoError(t, tc.Evl.Mul(ciphertext0, ciphertext1, ciphertext0))
			tc.Params.RingT().MulCoeffsBarrett(p0, p1, p0)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Mul/Ct/Pt/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values0, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, plaintext, _ := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext.Scale.Cmp(plaintext.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			require.NoError(t, tc.Evl.Mul(ciphertext, plaintext, ciphertext))
			tc.Params.RingT().MulCoeffsBarrett(p0, p1, p0)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Mul/Ct/Scalar/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

			scalar := tc.Params.PlaintextModulus() >> 1

			p := ring.Poly{Coeffs: [][]uint64{values}}

			require.NoError(t, tc.Evl.Mul(ciphertext, scalar, ciphertext))
			tc.Params.RingT().MulScalar(p, scalar, p)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Mul/Ct/Vector/Inplace", tc, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

			p := ring.Poly{Coeffs: [][]uint64{values}}

			require.NoError(t, tc.Evl.Mul(ciphertext, values, ciphertext))
			tc.Params.RingT().MulCoeffsBarrett(p, p, p)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/Square/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

			p := ring.Poly{Coeffs: [][]uint64{values}}

			require.NoError(t, tc.Evl.Mul(ciphertext, ciphertext, ciphertext))
			tc.Params.RingT().MulCoeffsBarrett(p, p, p)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/MulRelin/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			tc.Params.RingT().MulCoeffsBarrett(p0, p1, p0)

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

			receiver := NewCiphertext(tc.Params, 1, lvl)

			require.NoError(t, tc.Evl.MulRelin(ciphertext0, ciphertext1, receiver))
			require.NoError(t, tc.Evl.Rescale(receiver, receiver))

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, receiver, p0.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/MulThenAdd/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(2))
			values2, _, ciphertext2 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}
			p2 := ring.Poly{Coeffs: [][]uint64{values2}}

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)
			require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

			require.NoError(t, tc.Evl.MulThenAdd(ciphertext0, ciphertext1, ciphertext2))
			tc.Params.RingT().MulCoeffsBarrettThenAdd(p0, p1, p2)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext2, p2.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/MulThenAdd/Ct/Pt/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())
			values1, plaintext1, _ := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(2))
			values2, _, ciphertext2 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}
			p2 := ring.Poly{Coeffs: [][]uint64{values2}}

			require.True(t, ciphertext0.Scale.Cmp(plaintext1.Scale) != 0)
			require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

			require.NoError(t, tc.Evl.MulThenAdd(ciphertext0, plaintext1, ciphertext2))
			tc.Params.RingT().MulCoeffsBarrettThenAdd(p0, p1, p2)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext2, p2.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/MulThenAdd/Ct/Scalar/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

			scalar := tc.Params.PlaintextModulus() >> 1

			require.NoError(t, tc.Evl.MulThenAdd(ciphertext0, scalar, ciphertext1))
			tc.Params.RingT().MulScalarThenAdd(p0, scalar, p1)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext1, p1.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/MulThenAdd/Ct/Vector/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}

			scale := ciphertext1.Scale

			require.NoError(t, tc.Evl.MulThenAdd(ciphertext0, values1, ciphertext1))
			tc.Params.RingT().MulCoeffsBarrettThenAdd(p0, p1, p1)

			// Checks that output scale isn't changed
			require.True(t, scale.Equal(ciphertext1.Scale))

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext1, p1.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel {
		t.Run(name("Evaluator/MulRelinThenAdd/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
			if lvl == 0 {
				t.Skip("Skipping: Level = 0")
			}

			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())
			values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(2))
			values2, _, ciphertext2 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

			require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)
			require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

			p0 := ring.Poly{Coeffs: [][]uint64{values0}}
			p1 := ring.Poly{Coeffs: [][]uint64{values1}}
			p2 := ring.Poly{Coeffs: [][]uint64{values2}}

			require.NoError(t, tc.Evl.MulRelinThenAdd(ciphertext0, ciphertext1, ciphertext2))
			tc.Params.RingT().MulCoeffsBarrettThenAdd(p0, p1, p2)

			VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext2, p2.Coeffs[0], t)
		})
	}

	for _, lvl := range testLevel[:] {
		t.Run(name("Evaluator/Rescale", tc, lvl), func(t *testing.T) {

			ringT := tc.Params.RingT()

			values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

			printNoise := func(msg string, values []uint64, ct *rlwe.Ciphertext) {
				pt := NewPlaintext(tc.Params, ct.Level())
				pt.MetaData = ciphertext0.MetaData
				require.NoError(t, tc.Ecd.Encode(values0, pt))
				ct, err := tc.Evl.SubNew(ct, pt)
				require.NoError(t, err)
				vartmp, _, _ := rlwe.Norm(ct, tc.Dec)
				t.Logf("STD(noise) %s: %f\n", msg, vartmp)
			}

			if lvl != 0 {

				values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

				if *flagPrintNoise {
					printNoise("0x", values0, ciphertext0)
				}

				p0 := ring.Poly{Coeffs: [][]uint64{values0}}
				p1 := ring.Poly{Coeffs: [][]uint64{values1}}

				for i := 0; i < lvl; i++ {
					tc.Evl.MulRelin(ciphertext0, ciphertext1, ciphertext0)

					ringT.MulCoeffsBarrett(p0, p1, p0)

					if *flagPrintNoise {
						printNoise(fmt.Sprintf("%dx", i+1), p0.Coeffs[0], ciphertext0)
					}

				}

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)

				require.Nil(t, tc.Evl.Rescale(ciphertext0, ciphertext0))

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)

			} else {
				require.NotNil(t, tc.Evl.Rescale(ciphertext0, ciphertext0))
			}
		})
	}
}

func testEvaluatorBfv(tc *TestContext, t *testing.T) {
	testLevels := []int{0, tc.Params.MaxLevel()}

	t.Run("Evaluator", func(t *testing.T) {

		for _, lvl := range testLevels {
			t.Run(name("Mul/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
				if lvl == 0 {
					t.Skip("Skipping: Level = 0")
				}

				values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
				values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				p0 := ring.Poly{Coeffs: [][]uint64{values0}}
				p1 := ring.Poly{Coeffs: [][]uint64{values1}}

				require.NoError(t, tc.Evl.Mul(ciphertext0, ciphertext1, ciphertext0))
				tc.Params.RingT().MulCoeffsBarrett(p0, p1, p0)

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext0, p0.Coeffs[0], t)
			})
		}

		for _, lvl := range testLevels {
			t.Run(name("Mul/Ct/Pt/Inplace", tc, lvl), func(t *testing.T) {
				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
				values1, plaintext, _ := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

				require.True(t, ciphertext.Scale.Cmp(plaintext.Scale) != 0)

				p0 := ring.Poly{Coeffs: [][]uint64{values0}}
				p1 := ring.Poly{Coeffs: [][]uint64{values1}}

				require.NoError(t, tc.Evl.Mul(ciphertext, plaintext, ciphertext))
				tc.Params.RingT().MulCoeffsBarrett(p0, p1, p0)

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p0.Coeffs[0], t)
			})
		}

		for _, lvl := range testLevels {
			t.Run(name("Mul/Ct/Scalar/Inplace", tc, lvl), func(t *testing.T) {
				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

				scalar := tc.Params.PlaintextModulus() >> 1

				p := ring.Poly{Coeffs: [][]uint64{values}}

				require.NoError(t, tc.Evl.Mul(ciphertext, scalar, ciphertext))
				tc.Params.RingT().MulScalar(p, scalar, p)

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
			})
		}

		for _, lvl := range testLevels {
			t.Run(name("Square/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values, _, ciphertext := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())

				p := ring.Poly{Coeffs: [][]uint64{values}}

				require.NoError(t, tc.Evl.Mul(ciphertext, ciphertext, ciphertext))
				tc.Params.RingT().MulCoeffsBarrett(p, p, p)

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext, p.Coeffs[0], t)
			})
		}

		for _, lvl := range testLevels {
			t.Run(name("MulRelin/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
				values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

				p0 := ring.Poly{Coeffs: [][]uint64{values0}}
				p1 := ring.Poly{Coeffs: [][]uint64{values1}}

				tc.Params.RingT().MulCoeffsBarrett(p0, p1, p0)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				receiver := NewCiphertext(tc.Params, 1, lvl)

				require.NoError(t, tc.Evl.MulRelin(ciphertext0, ciphertext1, receiver))
				require.NoError(t, tc.Evl.Rescale(receiver, receiver))

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, receiver, p0.Coeffs[0], t)
			})
		}

		for _, lvl := range testLevels {
			t.Run(name("MulThenAdd/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())
				values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(2))
				values2, _, ciphertext2 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

				p0 := ring.Poly{Coeffs: [][]uint64{values0}}
				p1 := ring.Poly{Coeffs: [][]uint64{values1}}
				p2 := ring.Poly{Coeffs: [][]uint64{values2}}

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)
				require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

				require.NoError(t, tc.Evl.MulThenAdd(ciphertext0, ciphertext1, ciphertext2))
				tc.Params.RingT().MulCoeffsBarrettThenAdd(p0, p1, p2)

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext2, p2.Coeffs[0], t)
			})
		}

		for _, lvl := range testLevels {
			t.Run(name("MulThenAdd/Ct/Pt/Inplace", tc, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())
				values1, plaintext1, _ := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(2))
				values2, _, ciphertext2 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

				p0 := ring.Poly{Coeffs: [][]uint64{values0}}
				p1 := ring.Poly{Coeffs: [][]uint64{values1}}
				p2 := ring.Poly{Coeffs: [][]uint64{values2}}

				require.True(t, ciphertext0.Scale.Cmp(plaintext1.Scale) != 0)
				require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

				require.NoError(t, tc.Evl.MulThenAdd(ciphertext0, plaintext1, ciphertext2))
				tc.Params.RingT().MulCoeffsBarrettThenAdd(p0, p1, p2)

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext2, p2.Coeffs[0], t)
			})
		}

		for _, lvl := range testLevels {
			t.Run(name("MulThenAdd/Ct/Scalar/Inplace", tc, lvl), func(t *testing.T) {
				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(3))
				values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

				p0 := ring.Poly{Coeffs: [][]uint64{values0}}
				p1 := ring.Poly{Coeffs: [][]uint64{values1}}

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				scalar := tc.Params.PlaintextModulus() >> 1

				require.NoError(t, tc.Evl.MulThenAdd(ciphertext0, scalar, ciphertext1))
				tc.Params.RingT().MulScalarThenAdd(p0, scalar, p1)

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext1, p1.Coeffs[0], t)
			})
		}

		for _, lvl := range testLevels {
			t.Run(name("MulRelinThenAdd/Ct/Ct/Inplace", tc, lvl), func(t *testing.T) {
				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.DefaultScale())
				values1, _, ciphertext1 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(2))
				values2, _, ciphertext2 := NewTestVector(tc.Params, tc.Ecd, tc.Enc, lvl, tc.Params.NewScale(7))

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)
				require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

				p0 := ring.Poly{Coeffs: [][]uint64{values0}}
				p1 := ring.Poly{Coeffs: [][]uint64{values1}}
				p2 := ring.Poly{Coeffs: [][]uint64{values2}}

				require.NoError(t, tc.Evl.MulRelinThenAdd(ciphertext0, ciphertext1, ciphertext2))
				tc.Params.RingT().MulCoeffsBarrettThenAdd(p0, p1, p2)

				VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, ciphertext2, p2.Coeffs[0], t)
			})
		}
	})
}

var (
	// testInsecure are insecure parameters used for the sole purpose of fast testing.
	testInsecure = ParametersLiteral{
		LogN: 10,
		Q:    []uint64{0x3fffffa8001, 0x1000090001, 0x10000c8001, 0x10000f0001, 0xffff00001},
		P:    []uint64{0x7fffffd8001},
	}

	testPlaintextModulus = []uint64{0x101, 0xffc001}

	testParams = []ParametersLiteral{testInsecure}
)
