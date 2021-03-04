package dbfv

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/rlwe"
)

var t uint64

func Benchmark_DBFV_ThresholdProtocol(b *testing.B) {

	defaultParams := bfv.DefaultParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = bfv.DefaultParams[:2] // the short test runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = bfv.DefaultParams // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams bfv.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []bfv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	t = parties / 2

	for _, p := range defaultParams {
		params, err := bfv.NewParametersFromLiteral(p)
		var testCtx *testContext
		if testCtx, err = gentestContext(params); err != nil {
			panic(err)
		}
		print("-----VARYING N------\n")
		for N := uint64(3); N < 21; N++ {
			parties = N
			t = parties / 2
			benchThreshold(testCtx, b)
		}
		print("-----VARYING t------\n")
		parties = uint64(20)
		for th := uint64(2); th < parties; th++ {
			t = th
			benchThreshold(testCtx, b)
		}
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
		sk_t *rlwe.SecretKey
	}

	p := new(Party)
	p.s = sk0Shards[0]

	b.Run(testString("Thresholdizer/Init/", parties, testCtx.params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Thresholdizer = NewThresholdizer(testCtx.params)
			p.gen = p.Thresholdizer.AllocateShareGenPoly()
			p.Thresholdizer.InitShareGenPoly(p.gen, p.s, t)
			p.tsk = bfv.NewSecretKey(testCtx.params)
			p.sk_t = bfv.NewSecretKey(testCtx.params)
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
			p.CombinerCache = NewCombinerCache(p.Combiner, active_shamir_keys[0], shamir_keys)
		}
	})
	//Nothing is cached (simulates first decryption)
	b.Run(testString("Combiner/Combine/", parties, testCtx.params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner.GenFinalShare(active_shamir_keys, active_shamir_keys[0], p.tsk, p.sk_t)
		}
	})
	p.CombinerCache.ClearCache()
	p.CombinerCache.CacheInverses(active_shamir_keys[0], active_shamir_keys)
	// Everything is cached (simulates n-th decryption)
	b.Run(testString("Combiner/CombineCached/", parties, testCtx.params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			temp_sk := testCtx.ringQP.NewPoly()
			p.CombinerCache.GenFinalShare(p.tsk, p.sk_t)
			p.s.Value = temp_sk.CopyNew()
		}
	})
}
