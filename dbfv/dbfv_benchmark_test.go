package dbfv

import (
	"fmt"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"testing"
)

func Benchmark_DBFVScheme(b *testing.B) {

	paramSets := bfv.DefaultParams

	for _, params := range paramSets {

		bfvContext := bfv.NewBfvContext()
		if err := bfvContext.SetParameters(&params); err != nil {
			b.Error(err)
		}

		context := bfvContext.ContextQ()

		kgen := bfvContext.NewKeyGenerator()

		sk0, pk0 := kgen.NewKeyPair()
		sk1, pk1 := kgen.NewKeyPair()

		_, _ = pk0, pk1

		crpGenerator, err := ring.NewCRPGenerator(nil, context)
		if err != nil {
			b.Error(err)
		}

		crpGenerator.Seed([]byte{})

		// EKG
		/*
			parties := 2

			ekgV2Naive := NewEkgProtocolNaive(context)

			// [nParties][CrtDecomp][WDecomp][2]
			samples := make([][][][2]*ring.Poly, parties)
			for i := 0; i < parties; i++ {
				samples[i] = ekgV2Naive.GenSamples(sk0.Get(), pk0.Get())
			}

			aggregatedSamples := make([][][][2]*ring.Poly, parties)
			for i := 0; i < parties; i++ {
				aggregatedSamples[i] = ekgV2Naive.Aggregate(sk0.Get(), pk0.Get(), samples)
			}

			//EKG_V2_Naive_Round_0
			b.Run(fmt.Sprintf("params=%d/decomp=%d/EKG_Naive_Round0", params.N, bitDecomp), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					ekgV2Naive.GenSamples(sk0.Get(), pk1.Get())
				}
			})

			//EKG_V2_Naive_Round_1
			b.Run(fmt.Sprintf("params=%d/decomp=%d/EKG_Naive_Round1", params.N, bitDecomp), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					ekgV2Naive.Aggregate(sk0.Get(), pk1.Get(), samples)
				}
			})

			//EKG_V2_Naive_Round_2
			b.Run(fmt.Sprintf("params=%d/decomp=%d/EKG_Naive_Round2", params.N, bitDecomp), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					ekgV2Naive.Finalize(aggregatedSamples)
				}
			})
		*/
		// EKG
		crp := make([]*ring.Poly, len(context.Modulus))

		// TODO : replace by the .beta
		for i := 0; i < len(context.Modulus); i++ {
			crp[i] = crpGenerator.Clock()
		}

		ekg := NewEkgProtocol(bfvContext)
		share1, share2, share3 := ekg.AllocateShares()

		b.Run(fmt.Sprintf("params=%d//EKG_Round_1", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ekg.GenShareRoundOne(sk0.Get(), sk1.Get(), crp, share1)
			}
		})

		b.Run(fmt.Sprintf("params=%d//EKG_Round_2", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ekg.GenShareRoundTwo(share1, sk0.Get(), crp, share2)
			}
		})

		b.Run(fmt.Sprintf("params=%d//EKG_Round_3", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ekg.GenShareRoundThree(share2, sk1.Get(), sk1.Get(), share3)
			}
		})

		evk := bfvContext.NewRelinKeyEmpty(1)
		b.Run(fmt.Sprintf("params=%d//EKG_GenGey", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ekg.GenRelinearizationKey(share2, share3, evk)
			}
		})

		/*
			// CKG
			ckg := NewCKGProtocol(bfvContext)
			share := ckg.AllocateShares()
			crs := crpGenerator.Clock()

			b.Run(fmt.Sprintf("params=%d/CKG_Round_1", params.N), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					ckg.GenShare(sk0.Get(), crs, share)
				}
			})

			pk := bfvContext.NewPublicKey()
			b.Run(fmt.Sprintf("params=%d/CKG_GenKey", params.N), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					ckg.GenPublicKey(share, crs, pk)
					//ckg.GenPublicKey()
				}
			})

			// CKS
			sigmaSmudging := 3.19
			cks := NewCKSProtocol(bfvContext, sigmaSmudging)
			cksShare := cks.AllocateShare()
			ciphertext := bfvContext.NewRandomCiphertext(1)
			delta := bfvContext.ContextQ().NewPoly()
			bfvContext.ContextQ().Sub(sk0.Get(), sk1.Get(), delta)

			b.Run(fmt.Sprintf("params=%d/sigmaSmudging=%f/CKS_Round_1", params.N, sigmaSmudging), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					cks.GenShareDelta(delta, ciphertext, cksShare)
				}
			})

			b.Run(fmt.Sprintf("params=%d/sigmaSmudging=%f/CKS_KeySwitch", params.N, sigmaSmudging), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					cks.KeySwitch(cksShare, ciphertext, ciphertext)
				}
			})

			// PCKS
			pcks := NewPCKSProtocol(bfvContext, sigmaSmudging)
			pcksShare := pcks.AllocateShares()

			b.Run(fmt.Sprintf("params=%d/sigmaSmudging=%f/PCKS_Round_1", params.N, sigmaSmudging), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pcks.GenShare(sk0.Get(), pk0, ciphertext, pcksShare)
				}
			})

			b.Run(fmt.Sprintf("params=%d/sigmaSmudging=%f/PCKS_KeySwitch", params.N, sigmaSmudging), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pcks.KeySwitch(pcksShare, ciphertext, ciphertext)
				}
			})

		*/
	}
}
