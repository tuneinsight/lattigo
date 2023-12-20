package mheint

import (
	"encoding/json"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/heint"
	"github.com/tuneinsight/lattigo/v5/mhe"
)

func BenchmarkInteger(b *testing.B) {

	var err error

	paramsLiterals := testParams

	if *flagParamString != "" {
		var jsonParams heint.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		paramsLiterals = []heint.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals {

		for _, plaintextModulus := range testPlaintextModulus[:] {

			p.PlaintextModulus = plaintextModulus

			var params heint.Parameters
			if params, err = heint.NewParametersFromLiteral(p); err != nil {
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
		share mhe.RefreshShare
	}

	p := new(Party)
	var err error
	p.RefreshProtocol, err = NewRefreshProtocol(tc.params, tc.params.Xe())
	require.NoError(b, err)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(minLevel, maxLevel)

	ciphertext := heint.NewCiphertext(tc.params, 1, minLevel)

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
		opOut := heint.NewCiphertext(tc.params, 1, maxLevel)
		for i := 0; i < b.N; i++ {
			p.Finalize(ciphertext, crp, p.share, opOut)
		}
	})
}
