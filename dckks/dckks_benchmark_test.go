package dckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"testing"
)

func Benchmark_DCKKS(b *testing.B) {
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

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
		sk0Shards := params.sk0Shards

		crpGenerator := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
		crpGenerator.Seed([]byte{})
		crp := crpGenerator.Clock()

		type Party struct {
			*CKGProtocol
			s  *ring.Poly
			s1 CKGShare
		}

		p := new(Party)
		p.CKGProtocol = NewCKGProtocol(ckksContext)
		p.s = sk0Shards[0].Get()
		p.s1 = p.AllocateShares()

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			// Each party creates a new CKGProtocol instance
			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, crp, p.s1)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.s1, p.s1, p.s1)
			}
		})
	}
}

func benchRelinKeyGen(b *testing.B) {
	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
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
		p.RKGProtocol = NewEkgProtocol(ckksContext)
		p.u = p.RKGProtocol.NewEphemeralKey(1.0 / 3.0)
		p.s = sk0Shards[0].Get()
		p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()
		crpGenerator := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
		crpGenerator.Seed([]byte{})
		crp := make([]*ring.Poly, ckksContext.Beta())

		for i := uint64(0); i < ckksContext.Beta(); i++ {
			crp[i] = crpGenerator.Clock()
		}

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundOne(p.u, p.s, crp, p.share1)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round2Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundTwo(p.share1, p.s, crp, p.share2)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round2Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round3Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundThree(p.share2, p.u, p.s, p.share3)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round3Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundThree(p.share3, p.share3, p.share3)
			}
		})
	}
}

func benchRelinKeyGenNaive(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
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
		p.RKGProtocolNaive = NewRKGProtocolNaive(ckksContext)
		p.s = sk0Shards[0].Get()
		p.share1, p.share2 = p.AllocateShares()

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundOne(p.s, pk0.Get(), p.share1)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round2Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundTwo(p.share1, p.s, pk0.Get(), p.share2)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round2Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
			}
		})
	}

}

func benchKeySwitching(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
		sk0Shards := params.sk0Shards
		sk1Shards := params.sk1Shards

		type Party struct {
			*CKSProtocol
			s0    *ring.Poly
			s1    *ring.Poly
			share CKSShare
		}

		p := new(Party)
		p.CKSProtocol = NewCKSProtocol(ckksContext, 6.36)
		p.s0 = sk0Shards[0].Get()
		p.s1 = sk1Shards[0].Get()
		p.share = p.AllocateShare()

		ciphertext := ckksContext.NewRandomCiphertext(1, ckksContext.Levels()-1, ckksContext.Scale())

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s0, p.s1, ciphertext, p.share)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/KS", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.KeySwitch(p.share, ciphertext, ciphertext)
			}
		})
	}
}

func benchPublicKeySwitching(b *testing.B) {
	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
		sk0Shards := params.sk0Shards
		pk1 := params.pk1

		ciphertext := ckksContext.NewRandomCiphertext(1, ckksContext.Levels()-1, ckksContext.Scale())

		type Party struct {
			*PCKSProtocol
			s     *ring.Poly
			share PCKSShare
		}

		p := new(Party)
		p.PCKSProtocol = NewPCKSProtocol(ckksContext, 6.36)
		p.s = sk0Shards[0].Get()
		p.share = p.AllocateShares(ciphertext.Level())

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, pk1, ciphertext, p.share)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/KS", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.KeySwitch(p.share, ciphertext, ciphertext)
			}
		})
	}
}

func benchRotKeyGen(b *testing.B) {

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
		contextKeys := ckksContext.ContextKeys()
		sk0Shards := params.sk0Shards

		type Party struct {
			*RTGProtocol
			s     *ring.Poly
			share RTGShare
		}

		p := new(Party)
		p.RTGProtocol = NewRotKGProtocol(ckksContext)
		p.s = sk0Shards[0].Get()
		p.share = p.AllocateShare()

		crpGenerator := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
		crpGenerator.Seed([]byte{})
		crp := make([]*ring.Poly, ckksContext.Beta())

		for i := uint64(0); i < ckksContext.Beta(); i++ {
			crp[i] = crpGenerator.Clock()
		}

		mask := uint64((contextKeys.N >> 1) - 1)

		b.Run(testString("Round1/Gen", parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(ckks.RotationRight, uint64(i)&mask, sk0Shards[0].Get(), crp, &p.share)
			}
		})

		b.Run(testString("Round1/Agg", parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share, p.share, p.share)
			}
		})

		rotKey := ckksContext.NewRotationKeys()
		b.Run(testString("Finalize", parameters), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Finalize(p.share, crp, rotKey)
			}
		})
	}
}

func benchRefresh(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
		sk0Shards := params.sk0Shards

		levelStart := uint64(3)

		type Party struct {
			*RefreshProtocol
			s      *ring.Poly
			share1 RefreshShareDecrypt
			share2 RefreshShareRecrypt
		}

		p := new(Party)
		p.RefreshProtocol = NewRefreshProtocol(ckksContext)
		p.s = sk0Shards[0].Get()
		p.share1, p.share2 = p.AllocateShares(levelStart)

		crpGenerator := ring.NewCRPGenerator(nil, ckksContext.ContextQ())
		crpGenerator.Seed([]byte{})
		crp := crpGenerator.Clock()

		ciphertext := ckksContext.NewCiphertext(1, levelStart, ckksContext.Scale())

		ckksContext.ContextQ().UniformPoly(ciphertext.Value()[0])
		ckksContext.ContextQ().UniformPoly(ciphertext.Value()[1])

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Gen", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShares(p.s, levelStart, parties, ciphertext, crp, p.share1, p.share2)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Agg", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Aggregate(p.share1, p.share1, p.share1)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Decrypt", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Decrypt(ciphertext, p.share1)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Recode", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Recode(ciphertext)
			}
		})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Recrypt", parties, ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.Recrypt(ciphertext, crp, p.share2)
			}
		})
	}
}
