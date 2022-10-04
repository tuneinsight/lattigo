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

		var testCtx *testContext
		if testCtx, err = gentestContext(params, parties); err != nil {
			b.Fatal(err)
		}

		benchKeyswitching(testCtx, b)
		benchPublicKeySwitching(testCtx, b)
		benchRefresh(testCtx, b)
	}
}

func benchKeyswitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards

	type Party struct {
		*CKSProtocol
		s0    *rlwe.SecretKey
		s1    *rlwe.SecretKey
		share *drlwe.CKSShare
	}

	ciphertext := bfv.NewCiphertext(testCtx.params, 1)

	p := new(Party)
	p.CKSProtocol = NewCKSProtocol(testCtx.params, 6.36)
	p.s0 = sk0Shards[0]
	p.s1 = sk1Shards[0]
	p.share = p.AllocateShare()

	b.Run(testString("Keyswitching/Round1/Gen/", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s0, p.s1, ciphertext.Value[1], p.share)
		}
	})

	b.Run(testString("Keyswitching/Round1/Agg/", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Keyswitching/Finalize/", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(ciphertext, p.share, ciphertext)
		}
	})
}

func benchPublicKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1

	ciphertext := bfv.NewCiphertext(testCtx.params, 1)

	type Party struct {
		*PCKSProtocol
		s     *rlwe.SecretKey
		share *drlwe.PCKSShare
	}

	p := new(Party)
	p.PCKSProtocol = NewPCKSProtocol(testCtx.params, 6.36)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare()

	b.Run(testString("PublicKeySwitching/Round1/Gen/", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, pk1, ciphertext.Value[1], p.share)

		}
	})

	b.Run(testString("PublicKeySwitching/Round1/Agg/", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/Finalize/", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(ciphertext, p.share, ciphertext)
		}
	})
}

func benchRefresh(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*RefreshProtocol
		s     *rlwe.SecretKey
		share *RefreshShare
	}

	p := new(Party)
	p.RefreshProtocol = NewRefreshProtocol(testCtx.params, 3.2)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare()

	ciphertext := bfv.NewCiphertext(testCtx.params, 1)

	crp := p.SampleCRP(ciphertext.Level(), testCtx.crs)

	b.Run(testString("Refresh/Round1/Gen/", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, ciphertext.Value[1], crp, p.share)
		}
	})

	b.Run(testString("Refresh/Round1/Agg/", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Refresh/Finalize/", testCtx.NParties, testCtx.params), func(b *testing.B) {
		ctOut := bfv.NewCiphertext(testCtx.params, 1)
		for i := 0; i < b.N; i++ {
			p.Finalize(ciphertext, crp, p.share, ctOut)
		}
	})
}
