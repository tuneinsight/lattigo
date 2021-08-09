package dckks

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

func BenchmarkDCKKS(b *testing.B) {

	defaultParams := ckks.DefaultParams
	if testing.Short() {
		defaultParams = ckks.DefaultParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams ckks.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ckks.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	parties := 3

	for _, p := range defaultParams {

		params, err := ckks.NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}

		var testCtx *testContext
		if testCtx, err = genTestParams(params, parties); err != nil {
			panic(err)
		}

		// benchPublicKeyGen(testCtx, b)
		// benchRelinKeyGen(testCtx, b)
		// benchKeySwitching(testCtx, b)
		// benchPublicKeySwitching(testCtx, b)
		// benchRotKeyGen(testCtx, b)
		// benchRefresh(testCtx, b)
		// benchMaskedTransform(testCtx, b)

		// Varying N
		for N := 3; N < 6; N++ {
			t := N / 2
			benchThreshold(testCtx.params, N, t, b)
		}

		// Varying t
		N := 6
		for t := 2; t < N; t++ {
			benchThreshold(testCtx.params, N, t, b)

		}
	}
}

func benchPublicKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
	crp := crpGenerator.ReadNew()

	type Party struct {
		*CKGProtocol
		s  *rlwe.SecretKey
		s1 *drlwe.CKGShare
	}

	p := new(Party)
	p.CKGProtocol = NewCKGProtocol(params)
	p.s = sk0Shards[0]
	p.s1 = p.AllocateShares()

	b.Run(testString("PublicKeyGen/Gen/", testCtx.NParties, params), func(b *testing.B) {

		// Each party creates a new CKGProtocol instance
		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, crp, p.s1)
		}
	})

	b.Run(testString("PublicKeyGen/Agg/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.s1, p.s1, p.s1)
		}
	})

}

func benchRelinKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	type Party struct {
		*RKGProtocol
		ephSk  *rlwe.SecretKey
		sk     *rlwe.SecretKey
		share1 *drlwe.RKGShare
		share2 *drlwe.RKGShare
	}

	p := new(Party)
	p.RKGProtocol = NewRKGProtocol(params)
	p.sk = sk0Shards[0]
	p.ephSk, p.share1, p.share2 = p.RKGProtocol.AllocateShares()

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
	crp := make([]*ring.Poly, params.Beta())

	for i := 0; i < params.Beta(); i++ {
		crp[i] = crpGenerator.ReadNew()
	}

	b.Run(testString("RelinKeyGen/Round1Gen/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShareRoundOne(p.sk, crp, p.ephSk, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round1Agg/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share1, p.share1, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round2Gen/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShareRoundTwo(p.ephSk, p.sk, p.share1, crp, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Round2Agg/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share2, p.share2, p.share2)
		}
	})

}

func benchKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards
	params := testCtx.params

	ciphertext := ckks.NewCiphertextRandom(testCtx.prng, params, 1, params.MaxLevel(), params.Scale())

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

	b.Run(testString("KeySwitching/Gen/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s0, p.s1, ciphertext.Ciphertext, p.share)
		}
	})

	b.Run(testString("KeySwitching/Agg/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("KeySwitching/KS/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext.Ciphertext, ciphertext.Ciphertext)
		}
	})
}

func benchPublicKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1
	params := testCtx.params

	ciphertext := ckks.NewCiphertextRandom(testCtx.prng, params, 1, params.MaxLevel(), params.Scale())

	type Party struct {
		*PCKSProtocol
		s     *rlwe.SecretKey
		share *drlwe.PCKSShare
	}

	p := new(Party)
	p.PCKSProtocol = NewPCKSProtocol(params, 6.36)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("PublicKeySwitching/Gen/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, pk1, ciphertext.Ciphertext, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/Agg/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/KS/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext.Ciphertext, ciphertext.Ciphertext)
		}
	})
}

func benchRotKeyGen(testCtx *testContext, b *testing.B) {

	ringQP := testCtx.ringQP
	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	type Party struct {
		*RTGProtocol
		s     *rlwe.SecretKey
		share *drlwe.RTGShare
	}

	p := new(Party)
	p.RTGProtocol = NewRotKGProtocol(params)
	p.s = sk0Shards[0]
	p.share = p.AllocateShares()

	crpGenerator := ring.NewUniformSampler(testCtx.prng, ringQP)
	crp := make([]*ring.Poly, params.Beta())

	for i := 0; i < params.Beta(); i++ {
		crp[i] = crpGenerator.ReadNew()
	}
	galEl := params.GaloisElementForRowRotation()
	b.Run(testString("RotKeyGen/Round1/Gen/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, galEl, crp, p.share)
		}
	})

	b.Run(testString("RotKeyGen/Round1/Agg/", testCtx.NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share, p.share, p.share)
		}
	})

	rotKey := ckks.NewSwitchingKey(params)
	b.Run(testString("RotKeyGen/Finalize/", testCtx.NParties, params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenRotationKey(p.share, crp, rotKey)
		}
	})
}

func benchRefresh(testCtx *testContext, b *testing.B) {

	params := testCtx.params

	minLevel, logBound, ok := GetMinimumLevelForBootstrapping(128, params.Scale(), testCtx.NParties, params.Q())

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

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
		crp := crpGenerator.ReadNew()

		ciphertext := ckks.NewCiphertextRandom(testCtx.prng, params, 1, minLevel, params.Scale())

		b.Run(testString("Refresh/Round1/Gen", testCtx.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShares(p.s, logBound, params.LogSlots(), ciphertext, crp, p.share)
			}
		})

		b.Run(testString("Refresh/Round1/Agg", testCtx.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Refresh/Finalize", testCtx.NParties, params), func(b *testing.B) {
			ctOut := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale())
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

	minLevel, logBound, ok := GetMinimumLevelForBootstrapping(128, params.Scale(), testCtx.NParties, params.Q())

	if ok {

		sk0Shards := testCtx.sk0Shards

		type Party struct {
			*MaskedTransformProtocol
			s     *rlwe.SecretKey
			share *MaskedTransformShare
		}

		ciphertext := ckks.NewCiphertextRandom(testCtx.prng, params, 1, minLevel, params.Scale())

		p := new(Party)
		p.MaskedTransformProtocol = NewMaskedTransformProtocol(params, logBound, 3.2)
		p.s = sk0Shards[0]
		p.share = p.AllocateShare(ciphertext.Level(), params.MaxLevel())

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
		crp := crpGenerator.ReadNew()

		permute := func(ptIn, ptOut []*ring.Complex) {
			for i := range ptIn {
				ptOut[i][0].Mul(ptIn[i][0], ring.NewFloat(0.9238795325112867, logBound))
				ptOut[i][1].Mul(ptIn[i][1], ring.NewFloat(0.7071067811865476, logBound))
			}
		}

		b.Run(testString("Refresh&Transform/Round1/Gen", testCtx.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShares(p.s, logBound, params.LogSlots(), ciphertext, crp, permute, p.share)
			}
		})

		b.Run(testString("Refresh&Transform/Round1/Agg", testCtx.NParties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Refresh&Transform/Transform", testCtx.NParties, params), func(b *testing.B) {
			ctOut := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale())
			for i := 0; i < b.N; i++ {
				p.Transform(ciphertext, params.LogSlots(), permute, crp, p.share, ctOut)
			}
		})

	} else {
		b.Log("bench skipped : not enough level to ensure correctness and 128 bit security")
	}
}

func benchThreshold(params ckks.Parameters, NParties, t int, b *testing.B) {

	type Party struct {
		*Thresholdizer
		*Combiner
		*CombinerCache
		gen *drlwe.ShareGenPoly
		s   *rlwe.SecretKey
		tsk *rlwe.SecretKey
	}

	p := new(Party)
	p.s = ckks.NewSecretKey(params)

	b.Run(testString("Thresholdizer/Init/", NParties, params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Thresholdizer = NewThresholdizer(params)
			p.gen = p.Thresholdizer.AllocateShareGenPoly()
			p.Thresholdizer.InitShareGenPoly(p.gen, p.s, t)
			p.tsk = ckks.NewSecretKey(params)
		}
	})

	//Array of all shamir
	shamirPoint := make([]*drlwe.ThreshPublicKey, NParties)
	b.Run(testString("Thresholdizer/KeyGen/", NParties, params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := 0; j < NParties; j++ {
				shamirPoint[j] = p.Thresholdizer.GenKeyFromID(drlwe.PartyID{String: fmt.Sprintf("An arbitrary ID %d", j)})
			}
		}
	})

	b.Run(testString("Thresholdizer/Share/", NParties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			tmpShare := p.Thresholdizer.AllocateSecretShare()

			for j := 0; j < NParties; j++ {
				p.Thresholdizer.GenShareForParty(p.gen, shamirPoint[j], tmpShare)
			}

			for k := 0; k < NParties; k++ {
				p.Thresholdizer.AggregateShares(tmpShare, tmpShare, tmpShare)
			}
			p.Thresholdizer.GenThreshSecretKey(tmpShare, p.tsk)
		}
	})

	activePoints := shamirPoint[:t]
	b.Run(testString("Combiner/Init/", NParties, params)+fmt.Sprintf("threshold=%d", t)+fmt.Sprintf("precomputation=false"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner = NewCombiner(params, t)
		}
	})

	b.Run(testString("Combiner/Init/", NParties, params)+fmt.Sprintf("threshold=%d", t)+fmt.Sprintf("precomputation=true"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.CombinerCache = NewCombinerCache(p.Combiner, activePoints[0], activePoints)
			p.CombinerCache.CacheInverses(activePoints[0], activePoints)
		}
	})
	//Nothing is cached (simulates first decryption)
	b.Run(testString("Combiner/Combine/", NParties, params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.CombinerCache.ClearCache()
			p.Combiner.GenFinalShare(activePoints, activePoints[0], p.tsk, p.tsk)
		}
	})
	// Everything is cached (simulates n-th decryption)
	b.Run(testString("Combiner/CombineCached/", NParties, params)+fmt.Sprintf("threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner.GenFinalShare(activePoints, activePoints[0], p.tsk, p.tsk)
		}
	})
}
