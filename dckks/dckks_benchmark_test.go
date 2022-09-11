package dckks

import (
	"encoding/json"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
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

		var testCtx *testContext
		if testCtx, err = genTestParams(params, parties); err != nil {
			b.Fatal(err)
		}

		benchKeySwitching(testCtx, b)
		benchPublicKeySwitching(testCtx, b)
		benchRefresh(testCtx, b)
		benchMaskedTransform(testCtx, b)
	}
}

func benchKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards
	params := testCtx.params

	ciphertext := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())

	type Party struct {
		*CKSProtocol
		s0    *rlwe.SecretKey
		s1    *rlwe.SecretKey
		share *drlwe.CKSShare
	}

	p := new(Party)
	p.CKSProtocol = NewCKSProtocol(params, 6.36)
	p.s0 = sk0Shards[0]
	p.s1 = sk1Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("KeySwitching/Gen", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s0, p.s1, ciphertext.Value[1], p.share)
		}
	})

	b.Run(testString("KeySwitching/Agg", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("KeySwitching/KS", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(ciphertext, p.share, ciphertext)
		}
	})
}

func benchPublicKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1
	params := testCtx.params

	ciphertext := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())

	type Party struct {
		*PCKSProtocol
		s     *rlwe.SecretKey
		share *drlwe.PCKSShare
	}

	p := new(Party)
	p.PCKSProtocol = NewPCKSProtocol(params, 6.36)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("PublicKeySwitching/Gen", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, pk1, ciphertext.Value[1], p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/Agg", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/KS", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(ciphertext, p.share, ciphertext)
		}
	})
}

func benchRefresh(testCtx *testContext, b *testing.B) {

	params := testCtx.params

	minLevel, logBound, ok := GetMinimumLevelForBootstrapping(128, params.DefaultScale().(*ckks.Scale).Value, testCtx.NParties, params.Q())

	if ok {

		sk0Shards := testCtx.sk0Shards

		type Party struct {
			*RefreshProtocol
			s     *rlwe.SecretKey
			share *RefreshShare
		}

		p := new(Party)
		p.RefreshProtocol = NewRefreshProtocol(params, logBound, 3.2)
		p.s = sk0Shards[0]
		p.share = p.AllocateShare(minLevel, params.MaxLevel())

		ciphertext := ckks.NewCiphertext(params, 1, minLevel, &ckks.Scale{Value: params.DefaultScale().(*ckks.Scale).Value})

		crp := p.SampleCRP(params.MaxLevel(), testCtx.crs)

		b.Run(testString("Refresh/Round1/Gen", testCtx.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, logBound, params.LogSlots(), ciphertext.Value[1], ciphertext.Scale().(*ckks.Scale).Value, crp, p.share)
			}
		})

		b.Run(testString("Refresh/Round1/Agg", testCtx.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Refresh/Finalize", testCtx.NParties, params), func(b *testing.B) {
			ctOut := ckks.NewCiphertext(params, 1, params.MaxLevel(), &ckks.Scale{Value: params.DefaultScale().(*ckks.Scale).Value})
			for i := 0; i < b.N; i++ {
				p.Finalize(ciphertext, params.LogSlots(), crp, p.share, ctOut)
			}
		})

	} else {
		b.Log("bench skipped : not enough level to ensure correctness and 128 bit security")
	}
}

func benchMaskedTransform(testCtx *testContext, b *testing.B) {

	params := testCtx.params

	minLevel, logBound, ok := GetMinimumLevelForBootstrapping(128, params.DefaultScale().(*ckks.Scale).Value, testCtx.NParties, params.Q())

	if ok {

		sk0Shards := testCtx.sk0Shards

		type Party struct {
			*MaskedTransformProtocol
			s     *rlwe.SecretKey
			share *MaskedTransformShare
		}

		ciphertext := ckks.NewCiphertext(params, 1, minLevel, params.DefaultScale())

		p := new(Party)
		p.MaskedTransformProtocol, _ = NewMaskedTransformProtocol(params, params, logBound, 3.2)
		p.s = sk0Shards[0]
		p.share = p.AllocateShare(ciphertext.Level(), params.MaxLevel())

		crp := p.SampleCRP(params.MaxLevel(), testCtx.crs)

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

		b.Run(testString("Refresh&Transform/Round1/Gen", testCtx.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, p.s, logBound, params.LogSlots(), ciphertext.Value[1], ciphertext.Scale().(*ckks.Scale).Value, crp, transform, p.share)
			}
		})

		b.Run(testString("Refresh&Transform/Round1/Agg", testCtx.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Refresh&Transform/Transform", testCtx.NParties, params), func(b *testing.B) {
			ctOut := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())
			for i := 0; i < b.N; i++ {
				p.Transform(ciphertext, params.LogSlots(), transform, crp, p.share, ctOut)
			}
		})

	} else {
		b.Log("bench skipped : not enough level to ensure correctness and 128 bit security")
	}
}
