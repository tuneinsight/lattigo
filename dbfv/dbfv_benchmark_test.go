package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"testing"
)

func Benchmark_DBFV(b *testing.B) {
	b.Run("PublicKeyGen", benchPublicKeyGen)
	b.Run("RelinKeyGen", benchRelinKeyGen)
	b.Run("RelinKeyGenNaive", benchRelinKeyGenNaive)
	b.Run("KeySwitching", benchKeyswitching)
	b.Run("PublicKeySwitching", benchPublicKeySwitching)
	b.Run("RotKeyGen", benchRotKeyGen)
	b.Run("Refresh", benchRefresh)
}

func benchPublicKeyGen(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}
		crpGenerator := ring.NewUniformSampler(prng, testCtx.contextQP)

		crp := crpGenerator.SampleNew()

		type Party struct {
			*CKGProtocol
			s  *ring.Poly
			s1 CKGShare
		}

		p := new(Party)
		p.CKGProtocol = NewCKGProtocol(parameters)
		p.s = sk0Shards[0].Get()
		p.s1 = p.AllocateShares()

		b.Run(testString("Round1/Gen", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, crp, p.s1)
			}
		})

		b.Run(testString("Round1/Agg", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.s1, p.s1, p.s1)
			}
		})

		pk := bfv.NewPublicKey(parameters)
		b.Run(testString("Finalize", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenPublicKey(p.s1, crp, pk)
			}
		})

	}
}

func benchRelinKeyGen(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards

		type Party struct {
			*RKGProtocol
			u      *ring.Poly
			s      *ring.Poly
			share1 RKGShareRoundOne
			share2 RKGShareRoundTwo
			share3 RKGShareRoundThree

			rlk *bfv.EvaluationKey
		}

		p := new(Party)
		p.RKGProtocol = NewEkgProtocol(parameters)
		p.u = p.RKGProtocol.NewEphemeralKey(1.0 / 3.0)
		p.s = sk0Shards[0].Get()
		p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()
		p.rlk = bfv.NewRelinKey(parameters, 2)
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		crpGenerator := ring.NewUniformSampler(prng, testCtx.contextQP)

		crp := make([]*ring.Poly, parameters.Beta)

		for i := uint64(0); i < parameters.Beta; i++ {
			crp[i] = crpGenerator.SampleNew()
		}

		b.Run(testString("Round1/Gen", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenShareRoundOne(p.u, p.s, crp, p.share1)
			}
		})

		b.Run(testString("Round1/Agg", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
			}
		})

		b.Run(testString("Round2/Gen", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenShareRoundTwo(p.share1, p.s, crp, p.share2)
			}
		})

		b.Run(testString("Round2/Agg", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
			}
		})

		b.Run(testString("Round3/Gen", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenShareRoundThree(p.share2, p.u, p.s, p.share3)
			}
		})

		b.Run(testString("Round3/Agg", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundThree(p.share3, p.share3, p.share3)
			}
		})

		b.Run(testString("Finalize", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenRelinearizationKey(p.share2, p.share3, p.rlk)
			}
		})
	}
}

func benchRelinKeyGenNaive(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVTestContext(parameters)

		pk0 := params.pk0
		sk0Shards := params.sk0Shards

		type Party struct {
			*RKGProtocolNaive
			u      *ring.Poly
			s      *ring.Poly
			share1 RKGNaiveShareRoundOne
			share2 RKGNaiveShareRoundTwo

			rlk *bfv.EvaluationKey
		}

		p := new(Party)
		p.RKGProtocolNaive = NewRKGProtocolNaive(parameters)
		p.s = sk0Shards[0].Get()
		p.share1, p.share2 = p.AllocateShares()
		p.rlk = bfv.NewRelinKey(parameters, 2)

		b.Run(testString("Round1/Gen", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundOne(p.s, pk0.Get(), p.share1)
			}
		})

		b.Run(testString("Round1/Agg", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
			}
		})

		b.Run(testString("Round2/Gen", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundTwo(p.share1, p.s, pk0.Get(), p.share2)
			}
		})

		b.Run(testString("Round2/Agg", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
			}
		})

		b.Run(testString("Finalize", parties, parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenRelinearizationKey(p.share2, p.rlk)
			}
		})
	}
}

func benchKeyswitching(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVTestContext(parameters)

		sk0Shards := params.sk0Shards
		sk1Shards := params.sk1Shards

		type Party struct {
			*CKSProtocol
			s0    *ring.Poly
			s1    *ring.Poly
			share CKSShare
		}

		p := new(Party)
		p.CKSProtocol = NewCKSProtocol(parameters, 6.36)
		p.s0 = sk0Shards[0].Get()
		p.s1 = sk1Shards[0].Get()
		p.share = p.AllocateShare()

		ciphertext := bfv.NewCiphertextRandom(parameters, 1)

		b.Run(testString("Round1/Gen", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s0, p.s1, ciphertext, p.share)
			}
		})

		b.Run(testString("Round1/Agg", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Finalize", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.KeySwitch(p.share, ciphertext, ciphertext)
			}
		})
	}
}

func benchPublicKeySwitching(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVTestContext(parameters)

		sk0Shards := params.sk0Shards
		pk1 := params.pk1

		ciphertext := bfv.NewCiphertextRandom(parameters, 1)

		type Party struct {
			*PCKSProtocol
			s     *ring.Poly
			share PCKSShare
		}

		p := new(Party)
		p.PCKSProtocol = NewPCKSProtocol(parameters, 6.36)
		p.s = sk0Shards[0].Get()
		p.share = p.AllocateShares()

		b.Run(testString("Round1/Gen", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, pk1, ciphertext, p.share)

			}
		})

		b.Run(testString("Round1/Agg", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Finalize", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.KeySwitch(p.share, ciphertext, ciphertext)
			}
		})
	}
}

func benchRotKeyGen(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards

		type Party struct {
			*RTGProtocol
			s     *ring.Poly
			share RTGShare
		}

		p := new(Party)
		p.RTGProtocol = NewRotKGProtocol(parameters)
		p.s = sk0Shards[0].Get()
		p.share = p.AllocateShare()
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		crpGenerator := ring.NewUniformSampler(prng, testCtx.contextQP)
		crp := make([]*ring.Poly, parameters.Beta)

		for i := uint64(0); i < parameters.Beta; i++ {
			crp[i] = crpGenerator.SampleNew()
		}

		mask := (testCtx.n >> 1) - 1

		b.Run(testString("Round1/Gen", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(bfv.RotationRight, uint64(i)&mask, sk0Shards[0].Get(), crp, &p.share)
			}
		})

		b.Run(testString("Round1/Agg", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share, p.share, p.share)
			}
		})

		rotKey := bfv.NewRotationKeys()
		b.Run(testString("Finalize", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Finalize(p.share, crp, rotKey)
			}
		})

	}

}

func benchRefresh(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards

		type Party struct {
			*RefreshProtocol
			s     *ring.Poly
			share RefreshShare
		}

		p := new(Party)
		p.RefreshProtocol = NewRefreshProtocol(parameters)
		p.s = sk0Shards[0].Get()
		p.share = p.AllocateShares()
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		crpGenerator := ring.NewUniformSampler(prng, testCtx.contextQP)
		crp := crpGenerator.SampleNew()

		ciphertext := bfv.NewCiphertextRandom(parameters, 1)

		b.Run(testString("Round1/Gen", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShares(p.s, ciphertext, crp, p.share)
			}
		})

		b.Run(testString("Round1/Agg", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share, p.share, p.share)
			}
		})

		b.Run(testString("Finalize", parties, parameters), func(b *testing.B) {
			ctOut := bfv.NewCiphertext(parameters, 1)
			for i := 0; i < b.N; i++ {
				p.Finalize(ciphertext, crp, p.share, ctOut) // TODO: why does this fail ?
			}
		})
	}
}
