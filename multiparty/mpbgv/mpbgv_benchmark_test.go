package mpbgv

import (
	"encoding/json"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
)

func BenchmarkInteger(b *testing.B) {

	var err error

	paramsLiterals := testParams

	if *flagParamString != "" {
		var jsonParams bgv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		paramsLiterals = []bgv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals {

		for _, plaintextModulus := range testPlaintextModulus[:] {

			p.PlaintextModulus = plaintextModulus

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
}

func benchRefresh(tc *testContext, b *testing.B) {

	sk0Shards := tc.sk0Shards

	minLevel := 0
	maxLevel := tc.params.MaxLevel()

	type Party struct {
		RefreshProtocol
		s     *rlwe.SecretKey
		share multiparty.RefreshShare
	}

	p := new(Party)
	var err error
	p.RefreshProtocol, err = NewRefreshProtocol(tc.params, tc.params.Xe())
	require.NoError(b, err)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(minLevel, maxLevel)

	ciphertext := bgv.NewCiphertext(tc.params, 1, minLevel)

	crp := p.SampleCRP(maxLevel, tc.crs)

	b.Run(GetTestName("Refresh/Round1/Gen", tc.params, tc.NParties), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, ciphertext, crp, &p.share)
		}
	})

	b.Run(GetTestName("Refresh/Round1/Agg", tc.params, tc.NParties), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, &p.share)
		}
	})

	b.Run(GetTestName("Refresh/Finalize", tc.params, tc.NParties), func(b *testing.B) {
		opOut := bgv.NewCiphertext(tc.params, 1, maxLevel)
		for i := 0; i < b.N; i++ {
			p.Finalize(ciphertext, crp, p.share, opOut)
		}
	})
}
