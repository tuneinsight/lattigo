package dckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"testing"
)

func Benchmark_DCKKS(b *testing.B) {
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

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
		sk0Shards := params.sk0Shards

		crpGenerator, err := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
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
		p.CKGProtocol = NewCKGProtocol(ckksContext)
		p.s = sk0Shards[0].Get()
		p.s1 = p.AllocateShares()

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				// Each party creates a new CKGProtocol instance
				for i := 0; i < b.N; i++ {
					p.GenShare(p.s, crp, p.s1)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round2",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

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
		p.u, err = p.RKGProtocol.NewEphemeralKey(1.0 / 3.0)
		if err != nil {
			b.Error(err)
		}
		p.s = sk0Shards[0].Get()
		p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()

		crpGenerator, err := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
		if err != nil {
			b.Error(err)
		}
		crpGenerator.Seed([]byte{})
		crp := make([]*ring.Poly, ckksContext.Beta())

		for i := uint64(0); i < ckksContext.Beta(); i++ {
			crp[i] = crpGenerator.Clock()
		}

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1Gen",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.GenShareRoundOne(p.u, p.s, crp, p.share1)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1Agg",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round2Gen",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.GenShareRoundTwo(p.share1, p.s, crp, p.share2)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round2Agg",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round3Gen",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.GenShareRoundThree(p.share2, p.u, p.s, p.share3)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round3Agg",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

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

		ekgNaive := NewEkgProtocolNaive(ckksContext)

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					ekgNaive.GenSamples(sk0Shards[0].Get(), pk0.Get())
				}
			})

		samples := ekgNaive.GenSamples(sk0Shards[0].Get(), pk0.Get())

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round2",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					ekgNaive.Aggregate(sk0Shards[0].Get(), pk0.Get(), [][][2]*ring.Poly{samples})
				}
			})

		aggregatedSamples := ekgNaive.Aggregate(sk0Shards[0].Get(), pk0.Get(), [][][2]*ring.Poly{samples})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round3",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					ekgNaive.Finalize([][][2]*ring.Poly{aggregatedSamples})
				}
			})
	}
}

func benchKeyswitching(b *testing.B) {

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

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Gen",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.GenShare(p.s0, p.s1, ciphertext, p.share)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Agg",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.AggregateShares(p.share, p.share, p.share)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/KS",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

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

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Gen",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.GenShare(p.s, pk1, ciphertext, p.share)

				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Agg",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.AggregateShares(p.share, p.share, p.share)
				}
			})

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/KS",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					p.KeySwitch(p.share, ciphertext, ciphertext)
				}
			})
	}
}

func benchRotKeyGen(b *testing.B) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

		ckksContext := params.ckksContext
		sk0Shards := params.sk0Shards

		rkg := NewRKG(ckksContext)

		crpGenerator, err := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
		if err != nil {
			b.Error(err)
		}
		crpGenerator.Seed([]byte{})
		crp := make([]*ring.Poly, ckksContext.Beta())

		for i := uint64(0); i < ckksContext.Beta(); i++ {
			crp[i] = crpGenerator.Clock()
		}

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Gen",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					rkg.GenShareRotRow(sk0Shards[0].Get(), crp)
				}
			})

		shares := rkg.GenShareRotRow(sk0Shards[0].Get(), crp)

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Agg",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					rkg.AggregateRotRow([][]*ring.Poly{shares}, crp)
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

		crp := params.ckksContext.ContextQ().NewPoly()
		ciphertext := params.ckksContext.NewRandomCiphertext(1, 2, params.ckksContext.Scale())

		refreshShares := make([]*RefreshShares, parties)

		startLevel := uint64(2)

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Gen",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {

				for i := 0; i < b.N; i++ {
					refreshShares[0] = GenRefreshShares(sk0Shards[0], startLevel, parties, ckksContext, ciphertext.Value()[0], crp)
				}
			})

		for i := uint64(1); i < parties; i++ {
			refreshShares[i] = GenRefreshShares(sk0Shards[i], startLevel, parties, ckksContext, ciphertext.Value()[0], crp)
		}

		b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Agg",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(b *testing.B) {
				for i := 0; i < b.N; i++ {

					Refresh(ciphertext, refreshShares, ckksContext, crp)

					b.StopTimer()
					ciphertext.Value()[0].Coeffs = ciphertext.Value()[0].Coeffs[:startLevel+1]
					b.StartTimer()
				}
			})
	}
}
