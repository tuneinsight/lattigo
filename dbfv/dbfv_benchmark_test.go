package dbfv

import (
	"fmt"
	//"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"testing"
)

func Benchmark_DBFV(b *testing.B) {
	b.Run("PublicKeyGen", benchPublicKeyGen)
	b.Run("RelinKeyGen", benchRelinKeyGen)
	b.Run("RelinKeyGenNaive", benchRelinKeyGenNaive)
	b.Run("KeySwitching", benchKeyswitching)
	//b.Run("PublicKeySwitching", benchPublicKeySwitching)
	//b.Run("RotKeyGen", benchRotKeyGen)
	//b.Run("Refresh", benchRefresh)
}

func benchPublicKeyGen(b *testing.B) {

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		sk0Shards := params.sk0Shards

		crpGenerator, err := ring.NewCRPGenerator(nil, contextKeys)
		if err != nil {
			b.Error(err)
		}

		crpGenerator.Seed([]byte{})
		crp := crpGenerator.Clock()

		type Party struct {
			*CKGProtocol
			s  *ring.Poly
			s1 CKGShare
		}

		p := new(Party)
		p.CKGProtocol = NewCKGProtocol(bfvContext)
		p.s = sk0Shards[0].Get()
		p.s1 = p.AllocateShares()

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Gen", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s, crp, p.s1)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Agg", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.s1, p.s1, p.s1)
			}
		})
	}
}

func benchRelinKeyGen(b *testing.B) {

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
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
		p.RKGProtocol = NewEkgProtocol(bfvContext)
		p.u, err = p.RKGProtocol.NewEphemeralKey(1.0 / 3.0)
		if err != nil {
			b.Error(err)
		}
		p.s = sk0Shards[0].Get()
		p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()

		crpGenerator, err := ring.NewCRPGenerator(nil, bfvContext.ContextKeys())
		if err != nil {
			b.Error(err)
		}
		crpGenerator.Seed([]byte{})
		crp := make([]*ring.Poly, bfvContext.Beta())

		for i := uint64(0); i < bfvContext.Beta(); i++ {
			crp[i] = crpGenerator.Clock()
		}

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round1Gen", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenShareRoundOne(p.u, p.s, crp, p.share1)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round1Agg", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round2Gen", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenShareRoundTwo(p.share1, p.s, crp, p.share2)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round2Agg", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round3Gen", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.GenShareRoundThree(p.share2, p.u, p.s, p.share3)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round3Agg", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundThree(p.share3, p.share3, p.share3)
			}
		})
	}
}

func benchRelinKeyGenNaive(b *testing.B) {

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
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
		p.RKGProtocolNaive = NewRKGProtocolNaive(bfvContext)
		p.s = sk0Shards[0].Get()
		p.share1, p.share2 = p.AllocateShares()

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round1Gen", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundOne(p.s, pk0.Get(), p.share1)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round1Agg", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round2Gen", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShareRoundTwo(p.share1, p.s, pk0.Get(), p.share2)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round2Agg", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
			}
		})
	}
}

func benchKeyswitching(b *testing.B) {

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		sk0Shards := params.sk0Shards
		sk1Shards := params.sk1Shards

		type Party struct {
			*CKSProtocol
			s0    *ring.Poly
			s1    *ring.Poly
			share CKSShare
		}

		p := new(Party)
		p.CKSProtocol = NewCKSProtocol(bfvContext, 6.36)
		p.s0 = sk0Shards[0].Get()
		p.s1 = sk1Shards[0].Get()
		p.share = p.AllocateShare()

		ciphertext := bfvContext.NewRandomCiphertext(1)

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round1Gen", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.GenShare(p.s0, p.s1, ciphertext, p.share)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round1Agg", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.AggregateShares(p.share, p.share, p.share)
			}
		})

		b.Run(fmt.Sprintf("N=%d/logQ=%d/Round2KS", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				p.KeySwitch(p.share, ciphertext, ciphertext)
			}
		})
	}
}
