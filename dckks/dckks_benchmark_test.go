package dckks

import (
	"testing"

	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

func BenchmarkDCKKS(b *testing.B) {
	b.Run("PublicKeyGen", benchPublicKeyGen)
	b.Run("RelinKeyGen", benchRelinKeyGen)
	b.Run("RelinKeyGenNaive", benchRelinKeyGenNaive)
	b.Run("KeySwitching", benchKeySwitching)
	b.Run("PublicKeySwitching", benchPublicKeySwitching)
	b.Run("RotKeyGen", benchRotKeyGen)
	b.Run("Refresh", benchRefresh)
}

func benchPublicKeyGen(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		sk0Shards := params.sk0Shards
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQP)
		crp := crpGenerator.ReadNew()

		type Party struct {
			*CKGProtocol
			s  *ring.Poly
			s1 CKGShare
		}

		p := new(Party)
		p.CKGProtocol = NewCKGProtocol(parameters)
		p.s = sk0Shards[0].Get()
		p.s1 = p.AllocateShares()

		b.Run(testString("Gen/", parties, parameters), func(b *testing.B) {

			// Each party creates a new CKGProtocol instance
			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, crp, p.s1)
			}
		})

		b.Run(testString("Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.s1, p.s1, p.s1)
			}
		})
	}
}

func benchRelinKeyGen(b *testing.B) {
	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		sk0Shards := params.sk0Shards

		type Party struct {
			*RKGProtocol
			u      *ring.Poly
			s      *ring.Poly
			share1 RKGShareRoundOne
			share2 RKGShareRoundTwo
			share3 RKGShareRoundThree
		}

		p := new(Party)
		p.RKGProtocol = NewEkgProtocol(parameters)
		p.u = p.RKGProtocol.NewEphemeralKey()
		p.s = sk0Shards[0].Get()
		p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}
		crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQP)
		crp := make([]*ring.Poly, parameters.Beta)

		for i := uint64(0); i < parameters.Beta; i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		b.Run(testString("Round1Gen/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundOne(p.u, p.s, crp, p.share1)
			}
		})

		b.Run(testString("Round1Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
			}
		})

		b.Run(testString("Round2Gen/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundTwo(p.share1, p.s, crp, p.share2)
			}
		})

		b.Run(testString("Round2Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
			}
		})

		b.Run(testString("Round3Gen/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundThree(p.share2, p.u, p.s, p.share3)
			}
		})

		b.Run(testString("Round3Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundThree(p.share3, p.share3, p.share3)
			}
		})
	}
}

func benchRelinKeyGenNaive(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		pk0 := params.pk0
		sk0Shards := params.sk0Shards

		type Party struct {
			*RKGProtocolNaive
			u      *ring.Poly
			s      *ring.Poly
			share1 RKGNaiveShareRoundOne
			share2 RKGNaiveShareRoundTwo
		}

		p := new(Party)
		p.RKGProtocolNaive = NewRKGProtocolNaive(parameters)
		p.s = sk0Shards[0].Get()
		p.share1, p.share2 = p.AllocateShares()

		b.Run(testString("Round1Gen/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundOne(p.s, pk0.Get(), p.share1)
			}
		})

		b.Run(testString("Round1Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
			}
		})

		b.Run(testString("Round2Gen/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundTwo(p.share1, p.s, pk0.Get(), p.share2)
			}
		})

		b.Run(testString("Round2Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
			}
		})
	}

}

func benchKeySwitching(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

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

		ciphertext := ckks.NewCiphertextRandom(params.prng, parameters, 1, parameters.MaxLevel, parameters.Scale)

		b.Run(testString("Gen/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s0, p.s1, ciphertext, p.share)
			}
		})

		b.Run(testString("Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(testString("KS/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.KeySwitch(p.share, ciphertext, ciphertext)
			}
		})
	}
}

func benchPublicKeySwitching(b *testing.B) {
	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		sk0Shards := params.sk0Shards
		pk1 := params.pk1

		ciphertext := ckks.NewCiphertextRandom(params.prng, parameters, 1, parameters.MaxLevel, parameters.Scale)

		type Party struct {
			*PCKSProtocol
			s     *ring.Poly
			share PCKSShare
		}

		p := new(Party)
		p.PCKSProtocol = NewPCKSProtocol(parameters, 6.36)
		p.s = sk0Shards[0].Get()
		p.share = p.AllocateShares(ciphertext.Level())

		b.Run(testString("Gen/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, pk1, ciphertext, p.share)
			}
		})

		b.Run(testString("Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(testString("KS/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.KeySwitch(p.share, ciphertext, ciphertext)
			}
		})
	}
}

func benchRotKeyGen(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		contextKeys := params.dckksContext.contextQP
		sk0Shards := params.sk0Shards

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

		crpGenerator := ring.NewUniformSampler(prng, contextKeys)
		crp := make([]*ring.Poly, parameters.Beta)

		for i := uint64(0); i < parameters.Beta; i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		mask := uint64((contextKeys.N >> 1) - 1)

		b.Run(testString("Round1/Gen", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(ckks.RotationRight, uint64(i)&mask, sk0Shards[0].Get(), crp, &p.share)
			}
		})

		b.Run(testString("Round1/Agg", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share, p.share, p.share)
			}
		})

		rotKey := ckks.NewRotationKeys()
		b.Run(testString("Finalize", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Finalize(parameters, p.share, crp, rotKey)
			}
		})
	}
}

func benchRefresh(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		sk0Shards := params.sk0Shards
		contextQ := params.dckksContext.contextQ

		levelStart := uint64(3)

		type Party struct {
			*RefreshProtocol
			s      *ring.Poly
			share1 RefreshShareDecrypt
			share2 RefreshShareRecrypt
		}

		p := new(Party)
		p.RefreshProtocol = NewRefreshProtocol(parameters)
		p.s = sk0Shards[0].Get()
		p.share1, p.share2 = p.AllocateShares(levelStart)
		keyedPRNG, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		crpGenerator := ring.NewUniformSampler(keyedPRNG, contextQ)
		crp := crpGenerator.ReadNew()

		ciphertext := ckks.NewCiphertextRandom(params.prng, parameters, 1, levelStart, parameters.Scale)

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := ring.NewUniformSampler(prng, contextQ)

		uniformSampler.Read(ciphertext.Value()[0])
		uniformSampler.Read(ciphertext.Value()[1])

		b.Run(testString("Gen/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShares(p.s, levelStart, parties, ciphertext, crp, p.share1, p.share2)
			}
		})

		b.Run(testString("Agg/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share1, p.share1, p.share1)
			}
		})

		b.Run(testString("Decrypt/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Decrypt(ciphertext, p.share1)
			}
		})

		b.Run(testString("Recode/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Recode(ciphertext)
			}
		})

		b.Run(testString("Recrypt/", parties, parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Recrypt(ciphertext, crp, p.share2)
			}
		})
	}
}
