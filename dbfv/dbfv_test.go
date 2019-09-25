package dbfv

import (
	"fmt"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"log"
	"testing"
)

func Test_DBFVScheme(t *testing.T) {

	paramSets := bfv.DefaultParams[0:2]
	bitDecomps := []uint64{60}
	nParties := []int{5}

	//sigmaSmudging := 6.36

	for _, params := range paramSets {

		// nParties data indpendant element
		bfvContext := bfv.NewBfvContext()
		if err := bfvContext.SetParameters(&params); err != nil {
			log.Fatal(err)
		}

		kgen := bfvContext.NewKeyGenerator()

		evaluator, err := bfvContext.NewEvaluator()
		if err != nil {
			log.Fatal(err)
		}

		context := bfvContext.ContextQ()

		contextT := bfvContext.ContextT()

		encoder := bfvContext.NewBatchEncoder()

		coeffsWant := contextT.NewUniformPoly()
		plaintextWant := bfvContext.NewPlaintext()
		encoder.EncodeUint(coeffsWant.Coeffs[0], plaintextWant)

		ciphertextTest := bfvContext.NewCiphertext(1)

		for _, parties := range nParties {

			crpGenerators := make([]*CRPGenerator, parties)
			for i := 0; i < parties; i++ {
				crpGenerators[i], err = NewCRPGenerator(nil, context)
				if err != nil {
					log.Fatal(err)
				}
				crpGenerators[i].Seed([]byte{})
			}

			// SecretKeys
			sk0_shards := make([]*bfv.SecretKey, parties)
			sk1_shards := make([]*bfv.SecretKey, parties)
			tmp0 := context.NewPoly()
			tmp1 := context.NewPoly()

			for i := 0; i < parties; i++ {
				sk0_shards[i], _ = kgen.NewSecretKey(1.0 / 3)
				sk1_shards[i], _ = kgen.NewSecretKey(1.0 / 3)
				context.Add(tmp0, sk0_shards[i].Get(), tmp0)
				context.Add(tmp1, sk1_shards[i].Get(), tmp1)
			}

			sk0 := new(bfv.SecretKey)
			sk1 := new(bfv.SecretKey)

			sk0.Set(tmp0)
			sk1.Set(tmp1)

			// Publickeys
			pk0, err := kgen.NewPublicKey(sk0)
			if err != nil {
				log.Fatal(err)
			}

			pk1, err := kgen.NewPublicKey(sk1)
			if err != nil {
				log.Fatal(err)
			}

			// Encryptors
			encryptor_pk0, err := bfvContext.NewEncryptor(pk0, nil)
			if err != nil {
				log.Fatal(err)
			}

			//encryptor_pk1, err := bfvContext.NewEncryptor(pk1)
			//if err != nil {
			//	log.Fatal(err)
			//}

			// Decryptors
			decryptor_sk0, err := bfvContext.NewDecryptor(sk0)
			if err != nil {
				log.Fatal(err)
			}

			decryptor_sk1, err := bfvContext.NewDecryptor(sk1)
			if err != nil {
				log.Fatal(err)
			}

			// Reference ciphertext
			ciphertext, err := encryptor_pk0.EncryptFromPkNew(plaintextWant)
			if err != nil {
				log.Fatal(err)
			}

			coeffsMul := contextT.NewPoly()
			for i := 0; i < 1; i++ {
				res, _ := evaluator.MulNew(ciphertext, ciphertext)
				ciphertext = res.Ciphertext()
				contextT.MulCoeffs(coeffsWant, coeffsWant, coeffsMul)
			}

			t.Run(fmt.Sprintf("N=%d/Qi=%dx%d/CRS_PRNG", context.N, len(context.Modulus), 60), func(t *testing.T) {

				Ha, _ := NewPRNG([]byte{})
				Hb, _ := NewPRNG([]byte{})

				// Random 32 byte seed
				seed1 := []byte{0x48, 0xc3, 0x31, 0x12, 0x74, 0x98, 0xd3, 0xf2,
					0x7b, 0x15, 0x15, 0x9b, 0x50, 0xc4, 0x9c, 0x00,
					0x7d, 0xa5, 0xea, 0x68, 0x1f, 0xed, 0x4f, 0x99,
					0x54, 0xc0, 0x52, 0xc0, 0x75, 0xff, 0xf7, 0x5c}

				// New reseed of the PRNG after one clock cycle with the seed1
				seed2 := []byte{250, 228, 6, 63, 97, 110, 68, 153,
					147, 236, 236, 37, 152, 89, 129, 32,
					185, 5, 221, 180, 160, 217, 247, 201,
					211, 188, 160, 163, 176, 83, 83, 138}

				Ha.Seed(seed1)
				Hb.Seed(append(seed1, seed2...)) //Append works since blake2b hashes blocks of 512 bytes

				Ha.SetClock(256)
				Hb.SetClock(255)

				a := Ha.Clock()
				b := Hb.Clock()

				for i := 0; i < 32; i++ {
					if a[i] != b[i] {
						t.Errorf("error : error prng")
						break
					}
				}

				crs_generator_1, _ := NewCRPGenerator(nil, context)
				crs_generator_2, _ := NewCRPGenerator(nil, context)

				crs_generator_1.Seed(seed1)
				crs_generator_2.Seed(append(seed1, seed2...)) //Append works since blake2b hashes blocks of 512 bytes

				crs_generator_1.SetClock(256)
				crs_generator_2.SetClock(255)

				p0 := crs_generator_1.Clock()
				p1 := crs_generator_2.Clock()

				if bfvContext.ContextQ().Equal(p0, p1) != true {
					t.Errorf("error : crs prng generator")
				}
			})

			// EKG_Naive
			for _, bitDecomp := range bitDecomps {

				t.Run(fmt.Sprintf("N=%d/Qi=%dx%d/bitdecomp=%d/EKG", context.N, len(context.Modulus), 60, bitDecomp), func(t *testing.T) {

					bitLog := uint64((60 + (60 % bitDecomp)) / bitDecomp)

					// Each party instantiate an ekg naive protocole
					ekg := make([]*rkgProtocolState, parties)
					ephemeralKeys := make([]*ring.Poly, parties)

					crp := make([][]*ring.Poly, len(context.Modulus))
					for j := 0; j < len(context.Modulus); j++ {
						crp[j] = make([]*ring.Poly, bitLog)
						for u := uint64(0); u < bitLog; u++ {
							crp[j][u] = crpGenerators[0].Clock()
						}
					}


					for i := 0; i < parties; i++ {
						ekg[i] = NewEkgProtocol(bfvContext, bitDecomp)
						ephemeralKeys[i], _ = ekg[i].NewEphemeralKey(1.0 / 3)
					}

					rlk := test_EKG_Protocol(parties, ekg, sk0_shards, ephemeralKeys, crp)

					if err := evaluator.Relinearize(ciphertext, rlk, ciphertextTest); err != nil {
						log.Fatal(err)
					}

					plaintextTest, err := decryptor_sk0.DecryptNew(ciphertextTest)
					if err != nil {
						log.Fatal(err)
					}

					coeffsTest, err := encoder.DecodeUint(plaintextTest)
					if err != nil {
						log.Fatal(err)
					}

					if equalslice(coeffsMul.Coeffs[0], coeffsTest) != true {
						t.Errorf("error : ekg rlk bad decrypt")
					}

				})
			}

			// EKG_Naive
			for _, bitDecomp := range bitDecomps {

				t.Run(fmt.Sprintf("N=%d/Qi=%dx%d/bitdecomp=%d/EKG_Naive", context.N, len(context.Modulus), 60, bitDecomp), func(t *testing.T) {

					// Each party instantiate an ekg naive protocole
					ekgNaive := make([]*EkgProtocolNaive, parties)
					for i := 0; i < parties; i++ {
						ekgNaive[i] = NewEkgProtocolNaive(context, bitDecomp)
					}

					evk := test_EKG_Protocol_Naive(parties, sk0_shards, pk0, ekgNaive)

					rlk := new(bfv.EvaluationKey)
					rlk.SetRelinKeys([][][][2]*ring.Poly{evk[0]}, bitDecomp)

					if err := evaluator.Relinearize(ciphertext, rlk, ciphertextTest); err != nil {
						log.Fatal(err)
					}

					plaintextTest, err := decryptor_sk0.DecryptNew(ciphertextTest)
					if err != nil {
						log.Fatal(err)
					}

					coeffsTest, err := encoder.DecodeUint(plaintextTest)
					if err != nil {
						log.Fatal(err)
					}

					if equalslice(coeffsMul.Coeffs[0], coeffsTest) != true {
						t.Errorf("error : ekg_naive rlk bad decrypt")
					}
				})
			}

			t.Run(fmt.Sprintf("N=%d/Qi=%dx%d/CKG", context.N, len(context.Modulus), 60), func(t *testing.T) {

				crp := crpGenerators[0].Clock()

				type Party struct {
					*ckgProtocolState
					s *ring.Poly
					s1 ckgShare
				}

				ckgParties := make([]*Party, parties)
				for i := 0; i < parties; i++ {
					p := new(Party)
					p.ckgProtocolState = NewCKGProtocol(bfvContext)
					p.s = sk0_shards[i].Get()
					p.s1 = p.AllocateShares()
					ckgParties[i] = p
				}
				P0 := ckgParties[0]

				// Each party creates a new ckgProtocolState instance
				for i, p := range ckgParties {
					p.GenShare(p.s, crp, p.s1)
					if i > 0 {
						P0.AggregateShare(p.s1, P0.s1, P0.s1)
					}
				}

				pk := &bfv.PublicKey{}
				P0.GetAggregatedKey(P0.s1, crp, pk)



				// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
				encryptorTest, err := bfvContext.NewEncryptor(pk, nil)
				if err != nil {
					log.Fatal(err)
				}

				ciphertextTest, err := encryptorTest.EncryptFromPkNew(plaintextWant)

				if err != nil {
					log.Fatal(err)
				}

				plaintextTest, err := decryptor_sk0.DecryptNew(ciphertextTest)
				if err != nil {
					log.Fatal(err)
				}

				coeffsTest, err := encoder.DecodeUint(plaintextTest)
				if err != nil {
					log.Fatal(err)
				}

				if equalslice(coeffsWant.Coeffs[0], coeffsTest) != true {
					t.Errorf("error : ckg protocol, cpk encrypt/decrypt test")
				}

			})

			t.Run(fmt.Sprintf("N=%d/Qi=%dx%d/CKS", context.N, len(context.Modulus), 60), func(t *testing.T) {

				ciphertext, err := encryptor_pk0.EncryptFromPkNew(plaintextWant)
				if err != nil {
					log.Fatal(err)
				}

				ciphertexts := make([]*bfv.Ciphertext, parties)
				for i := 0; i < parties; i++ {
					ciphertexts[i] = ciphertext.CopyNew().Ciphertext()
				}

				// Each party creates its CKS instance with deltaSk = si-si'
				cks := make([]*CKS, parties)
				for i := 0; i < parties; i++ {
					cks[i] = NewCKS(sk0_shards[i].Get(), sk1_shards[i].Get(), context, 6.36)
				}

				// Each party computes its hi share from the shared ciphertext
				// Each party encodes its share and sends it to the other n-1 parties
				hi := make([]*ring.Poly, parties)
				for i := 0; i < parties; i++ {
					hi[i] = cks[i].KeySwitch(ciphertexts[i].Value()[1])
				}
				// Each party receive the shares n-1 shares from the other parties and decodes them
				for i := 0; i < parties; i++ {
					// Then keyswitch the ciphertext with the decoded shares
					cks[i].Aggregate(ciphertexts[i].Value()[0], hi)
				}

				for i := 0; i < parties; i++ {

					plaintextHave, _ := decryptor_sk1.DecryptNew(ciphertexts[i])

					coeffsTest, err := encoder.DecodeUint(plaintextHave)
					if err != nil {
						log.Fatal(err)
					}

					if equalslice(coeffsWant.Coeffs[0], coeffsTest) != true {
						t.Errorf("error : CKS")
					}

				}
			})

			t.Run(fmt.Sprintf("N=%d/Qi=%dx%d/PCKS", context.N, len(context.Modulus), 60), func(t *testing.T) {

				ciphertext, err := encryptor_pk0.EncryptFromPkNew(plaintextWant)
				if err != nil {
					log.Fatal(err)
				}

				ciphertexts := make([]*bfv.Ciphertext, parties)
				for i := 0; i < parties; i++ {
					ciphertexts[i] = ciphertext.CopyNew().Ciphertext()
				}

				pcks := make([]*PCKS, parties)
				for i := 0; i < parties; i++ {
					pcks[i] = NewPCKS(sk0_shards[i].Get(), pk1.Get(), context, 6.36)
				}

				hi := make([][2]*ring.Poly, parties)
				for i := 0; i < parties; i++ {
					hi[i] = pcks[i].KeySwitch(ciphertexts[i].Value()[1])
				}

				for i := 0; i < parties; i++ {
					pcks[i].Aggregate(ciphertexts[i].Value(), hi)
				}

				for i := 0; i < parties; i++ {
					plaintextHave, _ := decryptor_sk1.DecryptNew(ciphertexts[i])

					coeffsTest, err := encoder.DecodeUint(plaintextHave)
					if err != nil {
						log.Fatal(err)
					}

					if equalslice(coeffsWant.Coeffs[0], coeffsTest) != true {
						t.Errorf("error : PCKS")
					}
				}
			})
		}
	}
}

func test_EKG_Protocol_Naive(parties int, sk []*bfv.SecretKey, collectivePk *bfv.PublicKey, ekgNaive []*EkgProtocolNaive) [][][][2]*ring.Poly {

	// ROUND 0
	// Each party generates its samples
	samples := make([][][][2]*ring.Poly, parties)
	for i := 0; i < parties; i++ {
		samples[i] = ekgNaive[i].GenSamples(sk[i].Get(), collectivePk.Get())
	}

	// ROUND 1
	// Each party aggretates its sample with the other n-1 samples
	aggregatedSamples := make([][][][2]*ring.Poly, parties)
	for i := 0; i < parties; i++ {
		aggregatedSamples[i] = ekgNaive[i].Aggregate(sk[i].Get(), collectivePk.Get(), samples)
	}

	// ROUND 2
	// Each party aggregates sums its aggregatedSample with the other n-1 aggregated samples
	evk := make([][][][2]*ring.Poly, parties)
	for i := 0; i < parties; i++ {
		evk[i] = ekgNaive[i].Finalize(aggregatedSamples)
	}

	return evk
}

func test_EKG_Protocol(parties int, ekgProtocols []*rkgProtocolState, sk []*bfv.SecretKey, ephemeralKeys []*ring.Poly, crp [][]*ring.Poly) *bfv.EvaluationKey {

	type Party struct{*rkgProtocolState
		u *ring.Poly
		s *ring.Poly
		share1 rkgShareRoundOne
		share2 rkgShareRoundTwo
		share3 rkgShareRoundThree}

	rkgParties := make([]*Party, parties)


	for i := range rkgParties {
		p := new(Party)
		p.rkgProtocolState = ekgProtocols[i]
		p.u = ephemeralKeys[i]
		p.s = sk[i].Get()
		p.share1, p.share2, p.share3 = p.rkgProtocolState.AllocateShares()
		rkgParties[i] = p
	}

	P0 := rkgParties[0]

	// ROUND 1
	for i, p := range rkgParties {
		p.GenShareRoundOne(p.u, p.s, crp, p.share1)
		if i > 0 {
			P0.AggregateShareRoundOne(p.share1, P0.share1, P0.share1)
		}
	}


	//ROUND 2
	for i, p := range rkgParties {
		p.GenShareRoundTwo(P0.share1, p.s, crp, p.share2)
		if i > 0 {
			P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
		}
	}

	// ROUND 3
	for i, p := range rkgParties {
		p.GenShareRound3(P0.share2, p.u, p.s, p.share3)
		if i > 0 {
			P0.AggregateShareRound3(p.share3, P0.share3, P0.share3)
		}
	}

	evk := new(bfv.EvaluationKey)
	P0.GenRelinearizationKey(P0.share2, P0.share3, evk)
	return evk
}
