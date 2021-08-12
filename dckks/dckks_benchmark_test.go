package dckks

import (
	"encoding/json"
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

	for _, p := range defaultParams {

		params, err := ckks.NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}

		var testCtx *testContext
		if testCtx, err = genTestParams(params); err != nil {
			panic(err)
		}

		benchPublicKeyGen(testCtx, b)
		benchRelinKeyGen(testCtx, b)
		benchKeySwitching(testCtx, b)
		benchPublicKeySwitching(testCtx, b)
		benchRotKeyGen(testCtx, b)
		benchRefresh(testCtx, b)
		benchMaskedTransform(testCtx, b)
	}
}

func benchPublicKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	type Party struct {
		*CKGProtocol
		s   *rlwe.SecretKey
		s1  *drlwe.CKGShare
		crp drlwe.CKGCRP
	}

	p := new(Party)
	p.CKGProtocol = NewCKGProtocol(params)
	p.s = sk0Shards[0]
	p.s1, p.crp = p.AllocateShares()

	testCtx.crpGenerator.Read(p.crp)

	b.Run(testString("PublicKeyGen/Gen/", parties, params), func(b *testing.B) {

		// Each party creates a new CKGProtocol instance
		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, p.crp, p.s1)
		}
	})

	b.Run(testString("PublicKeyGen/Agg/", parties, params), func(b *testing.B) {

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
		crp    drlwe.RKGCRP
	}

	p := new(Party)
	p.RKGProtocol = NewRKGProtocol(params)
	p.sk = sk0Shards[0]
	p.ephSk, p.share1, p.share2, p.crp = p.RKGProtocol.AllocateShares()

	testCtx.crpGenerator.Read(p.crp)

	b.Run(testString("RelinKeyGen/Round1Gen/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShareRoundOne(p.sk, p.crp, p.ephSk, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round1Agg/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share1, p.share1, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round2Gen/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShareRoundTwo(p.ephSk, p.sk, p.share1, p.crp, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Round2Agg/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share2, p.share2, p.share2)
		}
	})

}

func benchKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards
	params := testCtx.params

	ciphertext := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale())

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

	b.Run(testString("KeySwitching/Gen/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s0, p.s1, ciphertext.Ciphertext, p.share)
		}
	})

	b.Run(testString("KeySwitching/Agg/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("KeySwitching/KS/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext.Ciphertext, ciphertext.Ciphertext)
		}
	})
}

func benchPublicKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1
	params := testCtx.params

	ciphertext := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale())

	type Party struct {
		*PCKSProtocol
		s     *rlwe.SecretKey
		share *drlwe.PCKSShare
	}

	p := new(Party)
	p.PCKSProtocol = NewPCKSProtocol(params, 6.36)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("PublicKeySwitching/Gen/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, pk1, ciphertext.Ciphertext, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/Agg/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/KS/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext.Ciphertext, ciphertext.Ciphertext)
		}
	})
}

func benchRotKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	type Party struct {
		*RTGProtocol
		s     *rlwe.SecretKey
		share *drlwe.RTGShare
		crp   drlwe.RTGCRP
	}

	p := new(Party)
	p.RTGProtocol = NewRotKGProtocol(params)
	p.s = sk0Shards[0]
	p.share, p.crp = p.AllocateShares()

	testCtx.crpGenerator.Read(p.crp)

	galEl := params.GaloisElementForRowRotation()
	b.Run(testString("RotKeyGen/Round1/Gen/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, galEl, p.crp, p.share)
		}
	})

	b.Run(testString("RotKeyGen/Round1/Agg/", parties, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share, p.share, p.share)
		}
	})

	rotKey := ckks.NewSwitchingKey(params)
	b.Run(testString("RotKeyGen/Finalize/", parties, params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenRotationKey(p.share, p.crp, rotKey)
		}
	})
}

func benchRefresh(testCtx *testContext, b *testing.B) {

	params := testCtx.params

	minLevel, logBound, ok := GetMinimumLevelForBootstrapping(128, params.Scale(), parties, params.Q())

	if ok {

		sk0Shards := testCtx.sk0Shards

		type Party struct {
			*RefreshProtocol
			s     *rlwe.SecretKey
			share *RefreshShare
			crp   drlwe.RefreshCRP
		}

		p := new(Party)
		p.RefreshProtocol = NewRefreshProtocol(params, logBound, 3.2)
		p.s = sk0Shards[0]
		p.share, p.crp = p.AllocateShare(minLevel, params.MaxLevel())

		ciphertext := ckks.NewCiphertext(params, 1, minLevel, params.Scale())

		testCtx.crpGenerator.Read(p.crp)

		b.Run(testString("Refresh/Round1/Gen", parties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShares(p.s, logBound, params.LogSlots(), ciphertext, p.crp, p.share)
			}
		})

		b.Run(testString("Refresh/Round1/Agg", parties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Refresh/Finalize", parties, params), func(b *testing.B) {
			ctOut := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale())
			for i := 0; i < b.N; i++ {
				p.Finalize(ciphertext, params.LogSlots(), p.crp, p.share, ctOut)
			}
		})

	} else {
		b.Log("bench skipped : not enough level to ensure correctness and 128 bit security")
	}
}

func benchMaskedTransform(testCtx *testContext, b *testing.B) {

	params := testCtx.params

	minLevel, logBound, ok := GetMinimumLevelForBootstrapping(128, params.Scale(), parties, params.Q())

	if ok {

		sk0Shards := testCtx.sk0Shards

		type Party struct {
			*MaskedTransformProtocol
			s     *rlwe.SecretKey
			share *MaskedTransformShare
			crp   drlwe.RefreshCRP
		}

		ciphertext := ckks.NewCiphertext(params, 1, minLevel, params.Scale())

		p := new(Party)
		p.MaskedTransformProtocol = NewMaskedTransformProtocol(params, logBound, 3.2)
		p.s = sk0Shards[0]
		p.share, p.crp = p.AllocateShare(ciphertext.Level(), params.MaxLevel())

		testCtx.crpGenerator.Read(p.crp)

		permute := func(ptIn, ptOut []*ring.Complex) {
			for i := range ptIn {
				ptOut[i][0].Mul(ptIn[i][0], ring.NewFloat(0.9238795325112867, logBound))
				ptOut[i][1].Mul(ptIn[i][1], ring.NewFloat(0.7071067811865476, logBound))
			}
		}

		b.Run(testString("Refresh&Transform/Round1/Gen", parties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShares(p.s, logBound, params.LogSlots(), ciphertext, p.crp, permute, p.share)
			}
		})

		b.Run(testString("Refresh&Transform/Round1/Agg", parties, params), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Refresh&Transform/Transform", parties, params), func(b *testing.B) {
			ctOut := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale())
			for i := 0; i < b.N; i++ {
				p.Transform(ciphertext, params.LogSlots(), permute, p.crp, p.share, ctOut)
			}
		})

	} else {
		b.Log("bench skipped : not enough level to ensure correctness and 128 bit security")
	}
}
