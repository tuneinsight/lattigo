package dbfv

import (
	"encoding/json"
	"testing"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func BenchmarkDBFV(b *testing.B) {

	var err error

	defaultParams := bfv.DefaultParams
	if testing.Short() {
		defaultParams = bfv.DefaultParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams bfv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParams = []bfv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	parties := 3

	for _, p := range defaultParams {
		var params bfv.Parameters
		if params, err = bfv.NewParametersFromLiteral(p); err != nil {
			b.Fatal(err)
		}

		var tc *testContext
		if tc, err = gentestContext(params, parties); err != nil {
			b.Fatal(err)
		}

		benchRefresh(tc, b)
	}
}

func benchRefresh(tc *testContext, b *testing.B) {

	sk0Shards := tc.sk0Shards

	type Party struct {
		*RefreshProtocol
		s     *rlwe.SecretKey
		share *drlwe.RefreshShare
	}

	p := new(Party)
	p.RefreshProtocol = NewRefreshProtocol(tc.params, 3.2)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(tc.params.MaxLevel(), tc.params.MaxLevel())

	ciphertext := bfv.NewCiphertext(tc.params, 1, tc.params.MaxLevel())

	crp := p.SampleCRP(ciphertext.Level(), tc.crs)

	b.Run(testString("Refresh/Round1/Gen", tc.NParties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, ciphertext, crp, p.share)
		}
	})

	b.Run(testString("Refresh/Round1/Agg", tc.NParties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Refresh/Finalize", tc.NParties, tc.params), func(b *testing.B) {
		ctOut := bfv.NewCiphertext(tc.params, 1, tc.params.MaxLevel())
		for i := 0; i < b.N; i++ {
			p.Finalize(ciphertext, crp, p.share, ctOut)
		}
	})
}
