package dbfv

//
//import (
//	"fmt"
//	"github.com/lca1/lattigo/bfv"
//	"github.com/lca1/lattigo/ring"
//	"log"
//	"testing"
//)
//
//func Benchmark_DBFVScheme(b *testing.B) {
//
//	paramSets := bfv.DefaultParams[3:4]
//	bitDecomps := []uint64{60}
//	nParties := []int{2}
//
//	for _, params := range paramSets {
//
//		bfvContext := bfv.NewBfvContext()
//		if err := bfvContext.SetParameters(&params); err != nil {
//			log.Fatal(err)
//		}
//
//		context := bfvContext.ContextQ()
//
//		kgen := bfvContext.NewKeyGenerator()
//
//		sk0, pk0, err := kgen.NewKeyPair(1.0 / 3)
//		if err != nil {
//			log.Fatal(err)
//		}
//
//		sk1, pk1, err := kgen.NewKeyPair(1.0 / 3)
//		if err != nil {
//			log.Fatal(err)
//		}
//
//		crpGenerator, err := NewCRPGenerator(nil, context)
//		if err != nil {
//			log.Fatal(err)
//		}
//
//		crpGenerator.Seed([]byte{})
//
//		// EKG_Naive
//		for _, parties := range nParties {
//
//			for _, bitDecomp := range bitDecomps {
//
//				ekgV2Naive := NewEkgProtocolNaive(context, bitDecomp)
//
//				// [nParties][CrtDecomp][WDecomp][2]
//				samples := make([][][][2]*ring.Poly, parties)
//				for i := 0; i < parties; i++ {
//					samples[i] = ekgV2Naive.GenSamples(sk0.Get(), pk0.Get())
//				}
//
//				aggregatedSamples := make([][][][2]*ring.Poly, parties)
//				for i := 0; i < parties; i++ {
//					aggregatedSamples[i] = ekgV2Naive.AggregateShares(sk0.Get(), pk0.Get(), samples)
//				}
//
//				//EKG_V2_Naive_Round_0
//				b.Run(fmt.Sprintf("params=%d/parties=%d/decomp=%d/EKG_Naive_Round0", params.N, parties, bitDecomp), func(b *testing.B) {
//					for i := 0; i < b.N; i++ {
//						ekgV2Naive.GenSamples(sk0.Get(), pk1.Get())
//					}
//				})
//
//				//EKG_V2_Naive_Round_1
//				b.Run(fmt.Sprintf("params=%d/parties=%d/decomp=%d/EKG_Naive_Round1", params.N, parties, bitDecomp), func(b *testing.B) {
//					for i := 0; i < b.N; i++ {
//						ekgV2Naive.AggregateShares(sk0.Get(), pk1.Get(), samples)
//					}
//				})
//
//				//EKG_V2_Naive_Round_2
//				b.Run(fmt.Sprintf("params=%d/parties=%d/decomp=%d/EKG_Naive_Round2", params.N, parties, bitDecomp), func(b *testing.B) {
//					for i := 0; i < b.N; i++ {
//						ekgV2Naive.Finalize(aggregatedSamples)
//					}
//				})
//			}
//		}
//
//		// EKG
//		for _, parties := range nParties {
//
//			for _, bitDecomp := range bitDecomps {
//
//				bitLog := uint64((60 + (60 % bitDecomp)) / bitDecomp)
//
//				crp := make([][]*ring.Poly, len(context.Modulus))
//
//				for i := 0; i < len(context.Modulus); i++ {
//					crp[i] = make([]*ring.Poly, bitLog)
//					for j := uint64(0); j < bitLog; j++ {
//						crp[i][j] = crpGenerator.Clock()
//					}
//				}
//
//
//				rkgParties := make([]struct{*rkgProtocol
//				share1 rkgShareRoundOne
//				share2 rkgShareRoundTwo
//				share3 rkgShareRoundThree}, parties)
//
//				for i := 0; i < parties; i++ {
//					rkgParties[i].rkgProtocol = NewEkgProtocol(bfvContext, bitDecomp)
//				}
//
//				//EKG_V2_Round_0
//				b.Run(fmt.Sprintf("params=%d/parties=%d/decomp=%d/EKG_Round0", params.N, parties, bitDecomp), func(b *testing.B) {
//					for i := 0; i < b.N; i++ {
//						rkgParties[0].GenShareRoundOne(sk0.Get(), sk1.Get(), crp, rkgParties[0].share1)
//					}
//				})
//				for _, party := range rkgParties[1:] {
//					party.GenShareRoundOne(sk0.Get(), sk1.Get(), crp, party.share1)
//				}
//
//
//				// TODO: benchmark that
//				for _, party := range rkgParties[1:] {
//					rkgParties[0].AggregateShareRoundOne(party.share1, rkgParties[0].share1, rkgParties[0].share1)
//				}
//
//				//EKG_V2_Round_1
//				b.Run(fmt.Sprintf("params=%d/parties=%d/decomp=%d/EKG_Round1", params.N, parties, bitDecomp), func(b *testing.B) {
//					for i := 0; i < b.N; i++ {
//						rkgParties[0].GenShareRoundTwo(sk1.Get(), crp)
//					}
//				})
//				for _, party := range rkgParties[1:] {
//					party.GenShareRoundTwo(sk1.Get(), crp)
//				}
//
//				// TODO: benchmark that
//				for _, party := range rkgParties {
//					rkgParties[0].AggregateShareRoundTwo(party.GetShareRoundTwo())
//				}
//				for _, party := range rkgParties[1:] {
//					party.AggregateShareRoundTwo(rkgParties[0].GetShareRoundTwo())
//				}
//
//
//				//EKG_V2_Round_2
//				b.Run(fmt.Sprintf("params=%d/parties=%d/decomp=%d/EKG_Round2", params.N, parties, bitDecomp), func(b *testing.B) {
//					for i := 0; i < b.N; i++ {
//						rkgParties[0].GenShareRound3(sk1.Get(), sk1.Get())
//					}
//				})
//				for _, party := range rkgParties[1:] {
//					party.GenShareRound3(sk0.Get(), sk1.Get())
//				}
//
//
//				for _, party := range rkgParties[1:] {
//					party.AggregateShareRound3(rkgParties[0].GetShareRoundThree())
//				}
//
//
//
//				//EKG_V2_Round_3
//				b.Run(fmt.Sprintf("params=%d/parties=%d/decomp=%d/EKG_Round3", params.N, parties, bitDecomp), func(b *testing.B) {
//					for i := 0; i < b.N; i++ {
//						// TODO: benchmark that
//						for _, party := range rkgParties {
//							rkgParties[0].AggregateShareRound3(party.GetShareRoundThree())
//						}
//					}
//				})
//			}
//		}
//
//		//ckgProtocol
//		for _, parties := range nParties {
//
//			ckgInstance := NewCKGProtocol(bfvContext, crpGenerator.Clock())
//			ckgInstance.GenShare(sk0.Get())
//
//			shares := make([]*ring.Poly, parties)
//			for i := 0; i < parties; i++ {
//				shares[i] = ckgInstance.GetShare()
//			}
//
//			// CKG_Round_0
//			b.Run(fmt.Sprintf("params=%d/parties=%d/CKG_Round0", params.N, parties), func(b *testing.B) {
//				for i := 0; i < b.N; i++ {
//					ckgInstance.GenShare(sk0.Get())
//				}
//			})
//
//			// CKG_Round_1
//			b.Run(fmt.Sprintf("params=%d/parties=%d/CKG_Round1", params.N, parties), func(b *testing.B) {
//				for i := 0; i < b.N; i++ {
//					ckgInstance.AggregateShares(shares)
//					//ckgInstance.GetPublicKey()
//				}
//			})
//		}
//
//		//cksProtocol
//		sigmaSmudging := 3.19
//		for _, parties := range nParties {
//
//			cksInstance := NewCKSProtocol(sk0.Get(), sk1.Get(), context, sigmaSmudging)
//
//			ciphertext := bfvContext.NewRandomCiphertext(1)
//
//			hi := make([]*ring.Poly, parties)
//			for i := 0; i < parties; i++ {
//				hi[i] = cksInstance.GenShare(ciphertext.Value()[1])
//			}
//
//			// CKS_Round_0
//			b.Run(fmt.Sprintf("params=%d/parties=%d/sigmaSmudging=%f/CKS_Round0", params.N, parties, sigmaSmudging), func(b *testing.B) {
//				for i := 0; i < b.N; i++ {
//					cksInstance.GenShare(ciphertext.Value()[1])
//				}
//			})
//
//			// CKS_Round_1
//			b.Run(fmt.Sprintf("params=%d/parties=%d/sigmaSmudging=%f/CKS_Round1", params.N, parties, sigmaSmudging), func(b *testing.B) {
//				for i := 0; i < b.N; i++ {
//					cksInstance.AggregateShares(ciphertext.Value()[0], hi)
//				}
//			})
//		}
//
//		//CKS_Trustless
//		for _, parties := range nParties {
//
//			pcks := NewPCKSProtocol(sk0.Get(), pk1.Get(), context, sigmaSmudging)
//
//			ciphertext := bfvContext.NewRandomCiphertext(1)
//
//			hi := make([][2]*ring.Poly, parties)
//			for i := 0; i < parties; i++ {
//				hi[i] = pcks.GenShare(ciphertext.Value()[1])
//			}
//
//			// CKS_Trustless_Round_0
//			b.Run(fmt.Sprintf("params=%d/parties=%d/sigmaSmudging=%f/PCKS_Round0", params.N, parties, sigmaSmudging), func(b *testing.B) {
//				for i := 0; i < b.N; i++ {
//					pcks.GenShare(ciphertext.Value()[1])
//				}
//			})
//
//			// CKS_Trustless_Round_1
//			b.Run(fmt.Sprintf("params=%d/parties=%d/sigmaSmudging=%f/PCKS_Round1", params.N, parties, sigmaSmudging), func(b *testing.B) {
//				for i := 0; i < b.N; i++ {
//					pcks.AggregateShares(ciphertext.Value(), hi)
//				}
//			})
//		}
//	}
//}
