package dbgv

import (
	"encoding/json"
	"testing"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func BenchmarkDBGV(b *testing.B) {

	var err error

	defaultParams := bgv.DefaultParams
	if testing.Short() {
		defaultParams = bgv.DefaultParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams bgv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParams = []bgv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {
		var params bgv.Parameters
		if params, err = bgv.NewParametersFromLiteral(p); err != nil {
			b.Fatal(err)
		}

		nParties := 3

		var tc *testContext
		if tc, err = gentestContext(nParties, params); err != nil {
			b.Fatal(err)
		}

		benchRefresh(tc, b)
	}
}

func benchRefresh(tc *testContext, b *testing.B) {

	sk0Shards := tc.sk0Shards

	minLevel := 0
	maxLevel := tc.params.MaxLevel()

	type Party struct {
		*RefreshProtocol
		s     *rlwe.SecretKey
		share *RefreshShare
	}

	p := new(Party)
	p.RefreshProtocol = NewRefreshProtocol(tc.params, 3.2)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(minLevel, maxLevel)

	ciphertext := bgv.NewCiphertext(tc.params, 1, minLevel)

	crp := p.SampleCRP(maxLevel, tc.crs)

	b.Run(testString("Refresh/Round1/Gen", tc.NParties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, ciphertext, ciphertext.Scale, crp, p.share)
		}
	})

	b.Run(testString("Refresh/Round1/Agg", tc.NParties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Refresh/Finalize", tc.NParties, tc.params), func(b *testing.B) {
		ctOut := bgv.NewCiphertext(tc.params, 1, maxLevel)
		for i := 0; i < b.N; i++ {
			p.Finalize(ciphertext, crp, p.share, ctOut)
		}
	})
}
