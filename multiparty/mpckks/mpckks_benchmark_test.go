package mpckks

import (
	"encoding/json"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func BenchmarkMultiPartyCKKS(b *testing.B) {

	var err error

	var testParams []ckks.ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, ckks.ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			b.Fatal(err)
		}
	default:
		testParams = testParamsLiteral
	}

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		for _, paramsLiteral := range testParams {

			paramsLiteral.RingType = ringType

			var params ckks.Parameters
			if params, err = ckks.NewParametersFromLiteral(paramsLiteral); err != nil {
				b.Fatal(err)
			}
			N := 3
			var tc *testContext
			if tc, err = genTestParams(params, N); err != nil {
				b.Fatal(err)
			}

			benchRefresh(tc, b)
			benchMaskedTransform(tc, b)
		}
	}
}

func benchRefresh(tc *testContext, b *testing.B) {

	params := tc.params

	minLevel, logBound, ok := GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q())

	if ok {

		sk0Shards := tc.sk0Shards

		type Party struct {
			RefreshProtocol
			s     *rlwe.SecretKey
			share multiparty.RefreshShare
		}

		p := new(Party)
		var err error
		p.RefreshProtocol, err = NewRefreshProtocol(params, logBound, params.Xe())
		require.NoError(b, err)
		p.s = sk0Shards[0]
		p.share = p.AllocateShare(minLevel, params.MaxLevel())

		ciphertext := ckks.NewCiphertext(params, 1, minLevel)

		crp := p.SampleCRP(params.MaxLevel(), tc.crs)

		b.Run(GetTestName("Refresh/Round1/Gen", tc.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, logBound, ciphertext, crp, &p.share)
			}
		})

		b.Run(GetTestName("Refresh/Round1/Agg", tc.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(&p.share, &p.share, &p.share)
			}
		})

		b.Run(GetTestName("Refresh/Finalize", tc.NParties, params), func(b *testing.B) {
			opOut := ckks.NewCiphertext(params, 1, params.MaxLevel())
			for i := 0; i < b.N; i++ {
				p.Finalize(ciphertext, crp, p.share, opOut)
			}
		})

	} else {
		b.Log("bench skipped : not enough level to ensure correctness and 128 bit security")
	}
}

func benchMaskedTransform(tc *testContext, b *testing.B) {

	params := tc.params

	minLevel, logBound, ok := GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q())

	if ok {

		sk0Shards := tc.sk0Shards

		type Party struct {
			MaskedLinearTransformationProtocol
			s     *rlwe.SecretKey
			share multiparty.RefreshShare
		}

		ciphertext := ckks.NewCiphertext(params, 1, minLevel)

		p := new(Party)
		p.MaskedLinearTransformationProtocol, _ = NewMaskedLinearTransformationProtocol(params, params, logBound, params.Xe())
		p.s = sk0Shards[0]
		p.share = p.AllocateShare(ciphertext.Level(), params.MaxLevel())

		crp := p.SampleCRP(params.MaxLevel(), tc.crs)

		transform := &MaskedLinearTransformationFunc{
			Decode: true,
			Func: func(coeffs []*bignum.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], bignum.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], bignum.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		b.Run(GetTestName("Refresh&Transform/Round1/Gen", tc.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, p.s, logBound, ciphertext, crp, transform, &p.share)
			}
		})

		b.Run(GetTestName("Refresh&Transform/Round1/Agg", tc.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(&p.share, &p.share, &p.share)
			}
		})

		b.Run(GetTestName("Refresh&Transform/Transform", tc.NParties, params), func(b *testing.B) {
			opOut := ckks.NewCiphertext(params, 1, params.MaxLevel())
			for i := 0; i < b.N; i++ {
				p.Transform(ciphertext, transform, crp, p.share, opOut)
			}
		})

	} else {
		b.Log("bench skipped : not enough level to ensure correctness and 128 bit security")
	}
}
