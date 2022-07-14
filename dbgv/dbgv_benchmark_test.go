package dbgv

import (
	"testing"

	"github.com/tuneinsight/lattigo/v3/bgv"
	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

func Benchmark_Dbgv(b *testing.B) {

	var err error
	for _, p := range TestParams[:] {
		var params bgv.Parameters
		if params, err = bgv.NewParametersFromLiteral(p); err != nil {
			b.Fatal(err)
		}

		var tc *testContext
		if tc, err = gentestContext(params); err != nil {
			b.Fatal(err)
		}

		benchPublicKeyGen(tc, b)
		benchRelinKeyGen(tc, b)
		benchKeyswitching(tc, b)
		benchPublicKeySwitching(tc, b)
		benchRotKeyGen(tc, b)
		benchRefresh(tc, b)
	}
}

func benchPublicKeyGen(tc *testContext, b *testing.B) {

	sk0Shards := tc.sk0Shards

	type Party struct {
		*CKGProtocol
		s  *rlwe.SecretKey
		s1 *drlwe.CKGShare
	}

	p := new(Party)
	p.CKGProtocol = NewCKGProtocol(tc.params)
	p.s = sk0Shards[0]
	p.s1 = p.AllocateShare()

	crp := p.SampleCRP(tc.crs)

	b.Run(testString("PublicKeyGen/Round1/Gen", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, crp, p.s1)
		}
	})

	b.Run(testString("PublicKeyGen/Round1/Agg", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShare(p.s1, p.s1, p.s1)
		}
	})

	pk := bgv.NewPublicKey(tc.params)
	b.Run(testString("PublicKeyGen/Finalize", parties, tc.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenPublicKey(p.s1, crp, pk)
		}
	})
}

func benchRelinKeyGen(tc *testContext, b *testing.B) {

	sk0Shards := tc.sk0Shards

	type Party struct {
		*RKGProtocol
		ephSk  *rlwe.SecretKey
		sk     *rlwe.SecretKey
		share1 *drlwe.RKGShare
		share2 *drlwe.RKGShare

		rlk *rlwe.RelinearizationKey
	}

	p := new(Party)
	p.RKGProtocol = NewRKGProtocol(tc.params)
	p.sk = sk0Shards[0]
	p.ephSk, p.share1, p.share2 = p.RKGProtocol.AllocateShare()
	p.rlk = bgv.NewRelinearizationKey(tc.params, 2)

	crp := p.SampleCRP(tc.crs)

	b.Run(testString("RelinKeyGen/Round1/Gen", parties, tc.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenShareRoundOne(p.sk, crp, p.ephSk, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round1/Agg", parties, tc.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.AggregateShare(p.share1, p.share1, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round2/Gen", parties, tc.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenShareRoundTwo(p.ephSk, p.sk, p.share1, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Round2/Agg", parties, tc.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.AggregateShare(p.share2, p.share2, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Finalize", parties, tc.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenRelinearizationKey(p.share1, p.share2, p.rlk)
		}
	})
}

func benchKeyswitching(tc *testContext, b *testing.B) {

	sk0Shards := tc.sk0Shards
	sk1Shards := tc.sk1Shards

	type Party struct {
		*CKSProtocol
		s0    *rlwe.SecretKey
		s1    *rlwe.SecretKey
		share *drlwe.CKSShare
	}

	ciphertext := bgv.NewCiphertext(tc.params, 1, tc.params.MaxLevel(), 1)

	p := new(Party)
	p.CKSProtocol = NewCKSProtocol(tc.params, 6.36)
	p.s0 = sk0Shards[0]
	p.s1 = sk1Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("Keyswitching/Round1/Gen", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s0, p.s1, ciphertext.Value[1], p.share)
		}
	})

	b.Run(testString("Keyswitching/Round1/Agg", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShare(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Keyswitching/Finalize", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(ciphertext, p.share, ciphertext)
		}
	})
}

func benchPublicKeySwitching(tc *testContext, b *testing.B) {

	sk0Shards := tc.sk0Shards
	pk1 := tc.pk1

	ciphertext := bgv.NewCiphertext(tc.params, 1, tc.params.MaxLevel(), 1)

	type Party struct {
		*PCKSProtocol
		s     *rlwe.SecretKey
		share *drlwe.PCKSShare
	}

	p := new(Party)
	p.PCKSProtocol = NewPCKSProtocol(tc.params, 6.36)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("PublicKeySwitching/Round1/Gen", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, pk1, ciphertext.Value[1], p.share)

		}
	})

	b.Run(testString("PublicKeySwitching/Round1/Agg", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShare(p.share, p.share, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/Finalize", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(ciphertext, p.share, ciphertext)
		}
	})
}

func benchRotKeyGen(tc *testContext, b *testing.B) {

	sk0Shards := tc.sk0Shards

	type Party struct {
		*RTGProtocol
		s     *rlwe.SecretKey
		share *drlwe.RTGShare
	}

	p := new(Party)
	p.RTGProtocol = NewRotKGProtocol(tc.params)
	p.s = sk0Shards[0]
	p.share = p.AllocateShare()

	crp := p.SampleCRP(tc.crs)

	b.Run(testString("RotKeyGen/Round1/Gen", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, tc.params.GaloisElementForRowRotation(), crp, p.share)
		}
	})

	b.Run(testString("RotKeyGen/Round1/Agg", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShare(p.share, p.share, p.share)
		}
	})

	rotKey := bgv.NewSwitchingKey(tc.params)
	b.Run(testString("RotKeyGen/Finalize", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenRotationKey(p.share, crp, rotKey)
		}
	})
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

	ciphertext := bgv.NewCiphertext(tc.params, 1, minLevel, 1)

	crp := p.SampleCRP(maxLevel, tc.crs)

	b.Run(testString("Refresh/Round1/Gen", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, ciphertext.Value[1], ciphertext.Scale, crp, p.share)
		}
	})

	b.Run(testString("Refresh/Round1/Agg", parties, tc.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShare(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Refresh/Finalize", parties, tc.params), func(b *testing.B) {
		ctOut := bgv.NewCiphertext(tc.params, 1, maxLevel, 1)
		for i := 0; i < b.N; i++ {
			p.Finalize(ciphertext, crp, p.share, ctOut)
		}
	})
}
