package dckks

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/rlwe"
)

var t uint64

func Benchmark_DBFV_ThresholdProtocol(b *testing.B) {

	defaultParams := ckks.DefaultParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = ckks.DefaultParams[:2] // the short test runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = ckks.DefaultParams // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams ckks.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ckks.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	t = parties / 2

	for _, p := range defaultParams {
		params, err := ckks.NewParametersFromLiteral(p)
		var testCtx *testContext
		if testCtx, err = genTestParams(params); err != nil {
			panic(err)
		}

		benchThreshold(testCtx, b)
	}
}

func benchThreshold(testCtx *testContext, b *testing.B) {
	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*Thresholdizer
		*Combiner
		*CombinerCache
		gen *drlwe.ShareGenPoly
		s   *rlwe.SecretKey
		tsk *rlwe.SecretKey
	}

	p := new(Party)
	p.s = sk0Shards[0]

	b.Run(testString("Thresholdizer/Init/", parties, testCtx.params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Thresholdizer = NewThresholdizer(testCtx.params)
			p.gen = p.Thresholdizer.AllocateShareGenPoly()
			p.Thresholdizer.InitShareGenPoly(p.gen, p.s, t)
			p.tsk = ckks.NewSecretKey(testCtx.params)
		}
	})

	//Array of all shamir
	shamir_keys := make([]*drlwe.ThreshPublicKey, parties)
	b.Run(testString("Thresholdizer/KeyGen/", parties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := uint64(0); j < parties; j++ {
				shamir_keys[j] = p.Thresholdizer.GenKeyFromID(drlwe.PartyID{fmt.Sprintf("An arbitrary ID %d", j)})
			}
		}
	})

	b.Run(testString("Thresholdizer/Share/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			tmp_share := p.Thresholdizer.AllocateSecretShare()

			for j := uint64(0); j < parties; j++ {
				p.Thresholdizer.GenShareForParty(p.gen, shamir_keys[j], tmp_share)
			}

			for k := uint64(0); k < parties; k++ {
				p.Thresholdizer.AggregateShares(tmp_share, tmp_share, tmp_share)
			}
			p.Thresholdizer.GenThreshSecretKey(tmp_share, p.tsk)
		}
	})

	active_shamir_keys := shamir_keys[:t]
	b.Run(testString("Combiner/Init/", parties, testCtx.params)+fmt.Sprintf("threshold=%d", t)+fmt.Sprintf("precomputation=false"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner = NewCombiner(testCtx.params, t)
		}
	})

	b.Run(testString("Combiner/Init/", parties, testCtx.params)+fmt.Sprintf("threshold=%d", t)+fmt.Sprintf("precomputation=true"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.CombinerCache = NewCombinerCache(p.Combiner, active_shamir_keys[0], active_shamir_keys)
			p.CombinerCache.CacheInverses(active_shamir_keys[0], active_shamir_keys)
		}
	})
	//Nothing is cached (simulates first decryption)
	b.Run(testString("Combiner/Combine/", parties, testCtx.params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.CombinerCache.ClearCache()
			temp_sk := testCtx.dckksContext.ringQP.NewPoly()
			p.Combiner.GenFinalShare(active_shamir_keys, active_shamir_keys[0], p.tsk, p.tsk)
			p.s.Value = temp_sk.CopyNew()
		}
	})
	// Everything is cached (simulates n-th decryption)
	b.Run(testString("Combiner/CombineCached/", parties, testCtx.params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			temp_sk := testCtx.dckksContext.ringQP.NewPoly()
			p.Combiner.GenFinalShare(active_shamir_keys, active_shamir_keys[0], p.tsk, p.tsk)
			p.s.Value = temp_sk.CopyNew()
		}
	})
}
