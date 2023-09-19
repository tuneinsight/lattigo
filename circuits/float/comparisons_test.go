package float_test

import (
	"math"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/circuits/float"
	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"

	"github.com/stretchr/testify/require"
)

// CoeffsMinimaxCompositePolynomialSign20Cheby is an example of composite minimax polynomial
// for the sign function that is able to distinguish between value with a delta of up to
// 2^{-alpha=30}, tolerates a scheme error of 2^{-35} and outputs a binary value (-1, or 1)
// of up to 20x4 bits of precision.
//
// It was computed with GenMinimaxCompositePolynomialForSign(256, 30, 35, []int{15, 15, 15, 17, 31, 31, 31, 31})
// which outputs a minimax composite polynomial of precision 21.926741, which is further composed with
// CoeffsSignX4Cheby to bring it to ~80bits of precision.
var CoeffsMinimaxCompositePolynomialSignAlpha30Err35Prec20x4Cheby = [][]string{
	{"0", "0.6371462957672043333", "0", "-0.2138032460610765328", "0", "0.1300439303835664499", "0", "-0.0948842756566191044", "0", "0.0760417811618939909", "0", "-0.0647714820920817557", "0", "0.0577904411211959048", "0", "-0.5275634328386103792"},
	{"0", "0.6371463830322414578", "0", "-0.2138032749880402509", "0", "0.1300439475440832118", "0", "-0.0948842877009570762", "0", "0.0760417903036533484", "0", "-0.0647714893343788749", "0", "0.0577904470018789283", "0", "-0.5275633669027163690"},
	{"0", "0.6371474873319408921", "0", "-0.2138036410457105809", "0", "0.1300441647026617059", "0", "-0.0948844401165889295", "0", "0.0760419059884502454", "0", "-0.0647715809823254389", "0", "0.0577905214191996406", "0", "-0.5275625325136631842"},
	{"0", "0.6370469776996076431", "0", "-0.2134526779726600620", "0", "0.1294300181775238920", "0", "-0.0939692999460324791", "0", "0.0747629355709698798", "0", "-0.0630298319949635571", "0", "0.0554299627688379896", "0", "-0.0504549111784642023", "0", "0.5242368268605847996"},
	{"0", "0.6371925153898374380", "0", "-0.2127272333844484291", "0", "0.1280350175397897124", "0", "-0.0918861831051024970", "0", "0.0719237384158242601", "0", "-0.0593247422790627989", "0", "0.0506973946536399213", "0", "-0.0444605229007162961", "0", "0.0397788020190944552", "0", "-0.0361705584687241925", "0", "0.0333397971860406254", "0", "-0.0310960060432036761", "0", "0.0293126335952747929", "0", "-0.0279042579223662982", "0", "0.0268135229627401517", "0", "-0.5128179323757194002"},
	{"0", "0.6484328404896112084", "0", "-0.2164688471885406655", "0", "0.1302737771018761402", "0", "-0.0934786176742356885", "0", "0.0731553324133884104", "0", "-0.0603252338481440981", "0", "0.0515366139595849853", "0", "-0.0451803385226980999", "0", "0.0404062758116036740", "0", "-0.0367241775307736352", "0", "0.0338327393147257876", "0", "-0.0315379870551266008", "0", "0.0297110181467332488", "0", "-0.0282647625290482803", "0", "0.0271406820054187399", "0", "-0.5041440308249296747"},
	{"0", "0.8988231150519633581", "0", "-0.2996064625122592138", "0", "0.1797645789317822353", "0", "-0.1284080039344265678", "0", "0.0998837306152582349", "0", "-0.0817422066647773587", "0", "0.0691963884439569899", "0", "-0.0600136111161848355", "0", "0.0530132660795356506", "0", "-0.0475133961913746909", "0", "0.0430936248086665091", "0", "-0.0394819050695222720", "0", "0.0364958013826412785", "0", "-0.0340100990129699835", "0", "0.0319381346687564699", "0", "-0.3095637759472512887"},
	{"0", "1.2654405107323937767", "0", "-0.4015427502443620045", "0", "0.2182109348265640036", "0", "-0.1341692540177466882", "0", "0.0852282854825304735", "0", "-0.0539043807248265057", "0", "0.0332611560159092728", "0", "-0.0197419082926337129", "0", "0.0111368708758574529", "0", "-0.0058990205011466309", "0", "0.0028925861201479251", "0", "-0.0012889673944941461", "0", "0.0005081425552893727", "0", "-0.0001696330470066833", "0", "0.0000440808328172753", "0", "-0.0000071549240608255"},
	float.CoeffsSignX4Cheby, // Quadruples the output precision (up to the scheme error)
}

func TestComparisons(t *testing.T) {

	paramsLiteral := float.TestPrec90

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		paramsLiteral.RingType = ringType

		if testing.Short() {
			paramsLiteral.LogN = 10
		}

		params, err := ckks.NewParametersFromLiteral(paramsLiteral)
		require.NoError(t, err)

		var tc *ckksTestContext
		if tc, err = genCKKSTestParams(params); err != nil {
			t.Fatal(err)
		}

		enc := tc.encryptorSk
		sk := tc.sk
		ecd := tc.encoder
		dec := tc.decryptor
		kgen := tc.kgen

		btp := bootstrapper.NewSecretKeyBootstrapper(params, sk)

		var galKeys []*rlwe.GaloisKey
		if params.RingType() == ring.Standard {
			galKeys = append(galKeys, kgen.GenGaloisKeyNew(params.GaloisElementForComplexConjugation(), sk))
		}

		eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk), galKeys...))
		polyEval := float.NewPolynomialEvaluator(params, eval)

		PWFEval := float.NewMinimaxCompositePolynomialEvaluator(params, eval, polyEval, btp)

		polys := float.NewMinimaxCompositePolynomial(CoeffsMinimaxCompositePolynomialSignAlpha30Err35Prec20x4Cheby)

		CmpEval := float.NewComparisonEvaluator(PWFEval, polys)

		threshold := bignum.NewFloat(math.Exp2(-30), params.EncodingPrecision())

		t.Run(GetTestName(params, "Sign"), func(t *testing.T) {

			values, _, ct := newCKKSTestVectors(tc, enc, complex(-1, 0), complex(1, 0), t)

			var sign *rlwe.Ciphertext
			sign, err = CmpEval.Sign(ct)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(sign), have))

			want := make([]*big.Float, params.MaxSlots())

			for i := range have {

				if new(big.Float).Abs(values[i][0]).Cmp(threshold) == -1 {
					want[i] = new(big.Float).Set(values[i][0])
				} else if values[i][0].Cmp(new(big.Float)) == -1 {
					want[i] = bignum.NewFloat(-1, params.EncodingPrecision())
				} else {
					want[i] = bignum.NewFloat(1, params.EncodingPrecision())
				}
			}

			ckks.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Step"), func(t *testing.T) {

			values, _, ct := newCKKSTestVectors(tc, enc, complex(-1, 0), complex(1, 0), t)

			var step *rlwe.Ciphertext
			step, err = CmpEval.Step(ct)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(step), have))

			want := make([]*big.Float, params.MaxSlots())

			for i := range have {

				if new(big.Float).Abs(values[i][0]).Cmp(threshold) == -1 {
					want[i] = new(big.Float).Set(values[i][0])
				} else if values[i][0].Cmp(new(big.Float)) == -1 {
					want[i] = bignum.NewFloat(0, params.EncodingPrecision())
				} else {
					want[i] = bignum.NewFloat(1, params.EncodingPrecision())
				}
			}

			ckks.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Max"), func(t *testing.T) {

			values0, _, ct0 := newCKKSTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)
			values1, _, ct1 := newCKKSTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)

			var max *rlwe.Ciphertext
			max, err = CmpEval.Max(ct0, ct1)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(max), have))

			want := make([]*big.Float, params.MaxSlots())

			for i := range have {

				if values0[i][0].Cmp(values1[i][0]) == -1 {
					want[i] = values1[i][0]
				} else {
					want[i] = values0[i][0]
				}
			}

			ckks.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Min"), func(t *testing.T) {

			values0, _, ct0 := newCKKSTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)
			values1, _, ct1 := newCKKSTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)

			var max *rlwe.Ciphertext
			max, err = CmpEval.Min(ct0, ct1)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(max), have))

			want := make([]*big.Float, params.MaxSlots())

			for i := range have {

				if values0[i][0].Cmp(values1[i][0]) == 1 {
					want[i] = values1[i][0]
				} else {
					want[i] = values0[i][0]
				}
			}

			ckks.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), nil, *printPrecisionStats, t)
		})
	}
}
