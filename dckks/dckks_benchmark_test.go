package dckks

import (
	"encoding/json"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func BenchmarkDCKKS(b *testing.B) {

	var err error

	defaultParams := ckks.DefaultParams
	if testing.Short() {
		defaultParams = ckks.DefaultParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams ckks.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParams = []ckks.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	parties := 3

	for _, p := range defaultParams {

		var params ckks.Parameters
		if params, err = ckks.NewParametersFromLiteral(p); err != nil {
			b.Fatal(err)
		}

		var tc *testContext
		if tc, err = genTestParams(params, parties); err != nil {
			b.Fatal(err)
		}

		benchRefresh(tc, b)
		benchMaskedTransform(tc, b)
	}
}

func benchRefresh(tc *testContext, b *testing.B) {

	params := tc.params

	minLevel, logBound, ok := GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q())

	if ok {

		sk0Shards := tc.sk0Shards

		type Party struct {
			*RefreshProtocol
			s     *rlwe.SecretKey
			share *drlwe.RefreshShare
		}

		p := new(Party)
		p.RefreshProtocol = NewRefreshProtocol(params, logBound, 3.2)
		p.s = sk0Shards[0]
		p.share = p.AllocateShare(minLevel, params.MaxLevel())

		ciphertext := ckks.NewCiphertext(params, 1, minLevel)

		crp := p.SampleCRP(params.MaxLevel(), tc.crs)

		b.Run(testString("Refresh/Round1/Gen", tc.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, logBound, params.LogSlots(), ciphertext, crp, p.share)
			}
		})

		b.Run(testString("Refresh/Round1/Agg", tc.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Refresh/Finalize", tc.NParties, params), func(b *testing.B) {
			ctOut := ckks.NewCiphertext(params, 1, params.MaxLevel())
			for i := 0; i < b.N; i++ {
				p.Finalize(ciphertext, params.LogSlots(), crp, p.share, ctOut)
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
			*MaskedTransformProtocol
			s     *rlwe.SecretKey
			share *drlwe.RefreshShare
		}

		ciphertext := ckks.NewCiphertext(params, 1, minLevel)

		p := new(Party)
		p.MaskedTransformProtocol, _ = NewMaskedTransformProtocol(params, params, logBound, 3.2)
		p.s = sk0Shards[0]
		p.share = p.AllocateShare(ciphertext.Level(), params.MaxLevel())

		crp := p.SampleCRP(params.MaxLevel(), tc.crs)

		transform := &MaskedTransformFunc{
			Decode: true,
			Func: func(coeffs []*ring.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], ring.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], ring.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		b.Run(testString("Refresh&Transform/Round1/Gen", tc.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, p.s, logBound, params.LogSlots(), ciphertext, crp, transform, p.share)
			}
		})

		b.Run(testString("Refresh&Transform/Round1/Agg", tc.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Refresh&Transform/Transform", tc.NParties, params), func(b *testing.B) {
			ctOut := ckks.NewCiphertext(params, 1, params.MaxLevel())
			for i := 0; i < b.N; i++ {
				p.Transform(ciphertext, params.LogSlots(), transform, crp, p.share, ctOut)
			}
		})

	} else {
		b.Log("bench skipped : not enough level to ensure correctness and 128 bit security")
	}
}
