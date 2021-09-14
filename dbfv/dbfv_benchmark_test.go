package dbfv

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

func Benchmark_DBFV(b *testing.B) {

	defaultParams := bfv.DefaultParams
	if testing.Short() {
		defaultParams = bfv.DefaultParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams bfv.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []bfv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	parties := 3

	for _, p := range defaultParams {
		params, err := bfv.NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}
		var testCtx *testContext
		if testCtx, err = gentestContext(params, parties); err != nil {
			panic(err)
		}

		benchPublicKeyGen(testCtx, b)
		benchRelinKeyGen(testCtx, b)
		benchKeyswitching(testCtx, b)
		benchPublicKeySwitching(testCtx, b)
		benchRotKeyGen(testCtx, b)
		benchRefresh(testCtx, b)

		// Varying N
		for N := 3; N < 20; N++ {
			t := N / 2
			benchThreshold(testCtx.params, N, t, b)
		}

		// Varying t
		N := 20
		for t := 2; t < N; t++ {
			benchThreshold(testCtx.params, N, t, b)

		}

	}
}

func benchPublicKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)

	crp := crpGenerator.ReadNew()

	type Party struct {
		*CKGProtocol
		s  *rlwe.SecretKey
		s1 *drlwe.CKGShare
	}

	p := new(Party)
	p.CKGProtocol = NewCKGProtocol(testCtx.params)
	p.s = sk0Shards[0]
	p.s1 = p.AllocateShares()

	b.Run(testString("PublicKeyGen/Round1/Gen", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, crp, p.s1)
		}
	})

	b.Run(testString("PublicKeyGen/Round1/Agg", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.s1, p.s1, p.s1)
		}
	})

	pk := bfv.NewPublicKey(testCtx.params)
	b.Run(testString("PublicKeyGen/Finalize", testCtx.NParties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenPublicKey(p.s1, crp, pk)
		}
	})
}

func benchRelinKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*RKGProtocol
		ephSk  *rlwe.SecretKey
		sk     *rlwe.SecretKey
		share1 *drlwe.RKGShare
		share2 *drlwe.RKGShare

		rlk *rlwe.RelinearizationKey
	}

	p := new(Party)
	p.RKGProtocol = NewRKGProtocol(testCtx.params)
	p.sk = sk0Shards[0]
	p.ephSk, p.share1, p.share2 = p.RKGProtocol.AllocateShares()
	p.rlk = bfv.NewRelinearizationKey(testCtx.params, 2)

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)

	crp := make([]*ring.Poly, testCtx.params.Beta())

	for i := 0; i < testCtx.params.Beta(); i++ {
		crp[i] = crpGenerator.ReadNew()
	}

	b.Run(testString("RelinKeyGen/Round1/Gen", testCtx.NParties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenShareRoundOne(p.sk, crp, p.ephSk, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round1/Agg", testCtx.NParties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share1, p.share1, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round2/Gen", testCtx.NParties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenShareRoundTwo(p.ephSk, p.sk, p.share1, crp, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Round2/Agg", testCtx.NParties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share2, p.share2, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Finalize", testCtx.NParties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenRelinearizationKey(p.share1, p.share2, p.rlk)
		}
	})
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

	ciphertext := bfv.NewCiphertextRandom(testCtx.prng, testCtx.params, 1)

	p := new(Party)
	p.CKSProtocol = NewCKSProtocol(testCtx.params, 6.36)
	p.s0 = sk0Shards[0]
	p.s1 = sk1Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("Keyswitching/Round1/Gen", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s0, p.s1, ciphertext.Ciphertext, p.share)
		}
	})

	b.Run(testString("Keyswitching/Round1/Agg", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Keyswitching/Finalize", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext.Ciphertext, ciphertext.Ciphertext)
		}
	})
}

func benchPublicKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1

	ciphertext := bfv.NewCiphertextRandom(testCtx.prng, testCtx.params, 1)

	type Party struct {
		*PCKSProtocol
		s     *rlwe.SecretKey
		share *drlwe.PCKSShare
	}

	p := new(Party)
	p.PCKSProtocol = NewPCKSProtocol(testCtx.params, 6.36)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("PublicKeySwitching/Round1/Gen", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, pk1, ciphertext.Ciphertext, p.share)

		}
	})

	b.Run(testString("PublicKeySwitching/Round1/Agg", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/Finalize", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext.Ciphertext, ciphertext.Ciphertext)
		}
	})
}

func benchRotKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*RTGProtocol
		s     *rlwe.SecretKey
		share *drlwe.RTGShare
	}

	p := new(Party)
	p.RTGProtocol = NewRotKGProtocol(testCtx.params)
	p.s = sk0Shards[0]
	p.share = p.AllocateShares()

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
	crp := make([]*ring.Poly, testCtx.params.Beta())

	for i := 0; i < testCtx.params.Beta(); i++ {
		crp[i] = crpGenerator.ReadNew()
	}

	b.Run(testString("RotKeyGen/Round1/Gen", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, testCtx.params.GaloisElementForRowRotation(), crp, p.share)
		}
	})

	b.Run(testString("RotKeyGen/Round1/Agg", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share, p.share, p.share)
		}
	})

	rotKey := bfv.NewSwitchingKey(testCtx.params)
	b.Run(testString("RotKeyGen/Finalize", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenRotationKey(p.share, crp, rotKey)
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

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
	crp := crpGenerator.ReadNew()

	ciphertext := bfv.NewCiphertextRandom(testCtx.prng, testCtx.params, 1)

	b.Run(testString("Refresh/Round1/Gen", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShares(p.s, ciphertext, crp, p.share)
		}
	})

	b.Run(testString("Refresh/Round1/Agg", testCtx.NParties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Refresh/Finalize", testCtx.NParties, testCtx.params), func(b *testing.B) {
		ctOut := bfv.NewCiphertext(testCtx.params, 1)
		for i := 0; i < b.N; i++ {
			p.Finalize(ciphertext, crp, p.share, ctOut)
		}
	})
}

func benchThreshold(params bfv.Parameters, NParties, t int, b *testing.B) {

	type Party struct {
		*drlwe.Thresholdizer
		drlwe.Combiner
		*drlwe.CachedCombiner
		gen *drlwe.ShamirPolynomial
		s   *rlwe.SecretKey
		tsk *drlwe.ShamirSecretShare
		sk  *rlwe.SecretKey
	}

	p := new(Party)
	p.s = bfv.NewSecretKey(params)

	b.Run(testString("Thresholdizer/Init/", NParties, params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Thresholdizer = drlwe.NewThresholdizer(params.Parameters)
			p.Thresholdizer.GenShamirPolynomial(t, p.s)
			p.tsk = p.Thresholdizer.AllocateThresholdSecretShare()
			p.sk = bfv.NewSecretKey(params)
		}
	})

	//Array of all shamir
	shamirPks := make([]*drlwe.ShamirPublicKey, NParties)
	b.Run(testString("Thresholdizer/KeyGen/", NParties, params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := 0; j < NParties; j++ {
				shamirPks[j] = p.Thresholdizer.GenShamirPublicKey(uint64(i) + 1)
			}
		}
	})

	b.Run(testString("Thresholdizer/Share/", NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			shamirShare := p.Thresholdizer.AllocateThresholdSecretShare()

			for j := 0; j < NParties; j++ {
				p.Thresholdizer.GenShamirSecretShare(shamirPks[j], p.gen, shamirShare)
			}

			for k := 0; k < NParties; k++ {
				p.Thresholdizer.AggregateShares(shamirShare, shamirShare, shamirShare)
			}
		}
	})

	activeShamirPks := shamirPks[:t]
	b.Run(testString("Combiner/Init/", NParties, params)+fmt.Sprintf("threshold=%d", t)+fmt.Sprintf("precomputation=false"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner = drlwe.NewCombiner(params.Parameters, t)
		}
	})

	b.Run(testString("Combiner/Init/", NParties, params)+fmt.Sprintf("threshold=%d", t)+fmt.Sprintf("precomputation=true"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.CachedCombiner = drlwe.NewCachedCombiner(params.Parameters, t)
		}
	})
	//Nothing is cached (simulates first decryption)
	b.Run(testString("Combiner/Combine/", NParties, params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner.GenAdditiveShare(activeShamirPks, activeShamirPks[0], p.tsk, p.sk)
		}
	})

	p.CachedCombiner.ClearCache()
	p.CachedCombiner.Precompute(activeShamirPks, activeShamirPks[0])
	// Everything is cached (simulates n-th decryption)
	b.Run(testString("Combiner/CombineCached/", NParties, params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.CachedCombiner.GenAdditiveShare(activeShamirPks, activeShamirPks[0], p.tsk, p.sk)
		}
	})
}
