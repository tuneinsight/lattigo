package dbfv

import (
	"fmt"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"log"
	"testing"
)

func check(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
	}
}

func Test_DBFVScheme(t *testing.T) {

	paramSets := bfv.DefaultParams[0:2]
	bitDecomps := []uint64{60}
	nParties := []int{5}

	//sigmaSmudging := 6.36

	for _, params := range paramSets {

		// nParties data indpendant element
		bfvContext := bfv.NewBfvContext()
		if err := bfvContext.SetParameters(&params); err != nil {
			t.Error(err)
		}

		kgen := bfvContext.NewKeyGenerator()

		evaluator := bfvContext.NewEvaluator()

		context := bfvContext.ContextQ()

		contextT := bfvContext.ContextT()

		encoder, err := bfvContext.NewBatchEncoder()
		check(t, err)

		coeffsWant := contextT.NewUniformPoly()
		plaintextWant := bfvContext.NewPlaintext()
		err = encoder.EncodeUint(coeffsWant.Coeffs[0], plaintextWant)
		check(t, err)

		for _, parties := range nParties {

			crpGenerators := make([]*CRPGenerator, parties)
			for i := 0; i < parties; i++ {
				crpGenerators[i], err = NewCRPGenerator(nil, context)
				if err != nil {
					t.Error(err)
				}
				crpGenerators[i].Seed([]byte{})
			}

			// SecretKeys
			sk0_shards := make([]*bfv.SecretKey, parties)
			sk1_shards := make([]*bfv.SecretKey, parties)
			tmp0 := context.NewPoly()
			tmp1 := context.NewPoly()

			for i := 0; i < parties; i++ {
				sk0_shards[i] = kgen.NewSecretKey()
				sk1_shards[i] = kgen.NewSecretKey()
				context.Add(tmp0, sk0_shards[i].Get(), tmp0)
				context.Add(tmp1, sk1_shards[i].Get(), tmp1)
			}

			sk0 := new(bfv.SecretKey)
			sk1 := new(bfv.SecretKey)

			sk0.Set(tmp0)
			sk1.Set(tmp1)

			// Publickeys
			pk0 := kgen.NewPublicKey(sk0)
			pk1 := kgen.NewPublicKey(sk1)

			// Encryptors
			encryptor_pk0, err := bfvContext.NewEncryptorFromPk(pk0)
			check(t, err)

			//encryptor_pk1, err := bfvContext.NewEncryptor(pk1)
			//if err != nil {
			//	t.Error(err)
			//}

			// Decryptors
			decryptor_sk0, err := bfvContext.NewDecryptor(sk0)
			check(t, err)

			decryptor_sk1, err := bfvContext.NewDecryptor(sk1)
			check(t, err)

			// Reference ciphertext
			ciphertext, err := encryptor_pk0.EncryptNew(plaintextWant)
			check(t, err)

			ciphertextMul, _ := evaluator.MulNew(ciphertext, ciphertext)
			coeffsMul := contextT.NewPoly()
			contextT.MulCoeffs(coeffsWant, coeffsWant, coeffsMul)

			t.Run(fmt.Sprintf("N=%d/logQ=%d/CRS_PRNG", context.N, context.ModulusBigint.Value.BitLen()), func(t *testing.T) {

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

			// EKG
			for _, bitDecomp := range bitDecomps {

				t.Run(fmt.Sprintf("N=%d/logQ=%d/bitdecomp=%d/EKG", context.N, context.ModulusBigint.Value.BitLen(), bitDecomp), func(t *testing.T) {

					bitLog := uint64((60 + (60 % bitDecomp)) / bitDecomp)

					// Each party instantiate an ekg naive protocole
					ekg := make([]*RKGProtocol, parties)
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

					rlk := test_EKG_Protocol(bfvContext, parties, bitDecomp, ekg, sk0_shards, ephemeralKeys, crp)
					res := bfvContext.NewCiphertext(1)
					if err := evaluator.Relinearize(ciphertextMul, rlk, res); err != nil {
						t.Error(err)
					}

					if equalslice(coeffsMul.Coeffs[0], encoder.DecodeUint(decryptor_sk0.DecryptNew(res))) != true {
						t.Errorf("error : ekg rlk bad decrypt")
					}

				})
			}

			// EKG_Naive
			for _, bitDecomp := range bitDecomps {

				t.Run(fmt.Sprintf("N=%d/logQ=%d/bitdecomp=%d/EKG_Naive", context.N, context.ModulusBigint.Value.BitLen(), bitDecomp), func(t *testing.T) {

					// Each party instantiate an ekg naive protocole
					ekgNaive := make([]*EkgProtocolNaive, parties)
					for i := 0; i < parties; i++ {
						ekgNaive[i] = NewEkgProtocolNaive(context, bitDecomp)
					}

					evk := test_EKG_Protocol_Naive(parties, sk0_shards, pk0, ekgNaive)

					rlk := new(bfv.EvaluationKey)
					rlk.SetRelinKeys([][][][2]*ring.Poly{evk[0]}, bitDecomp)

					res := bfvContext.NewCiphertext(1)
					if err := evaluator.Relinearize(ciphertextMul, rlk, res); err != nil {
						t.Error(err)
					}

					if equalslice(coeffsMul.Coeffs[0], encoder.DecodeUint(decryptor_sk0.DecryptNew(res))) != true {
						t.Errorf("error : ekg_naive rlk bad decrypt")
					}
				})
			}

			t.Run(fmt.Sprintf("N=%d/logQ=%d/CKG", context.N, context.ModulusBigint.Value.BitLen()), func(t *testing.T) {

				crp := crpGenerators[0].Clock()

				type Party struct {
					*CKGProtocol
					s  *ring.Poly
					s1 CKGShare
				}

				ckgParties := make([]*Party, parties)
				for i := 0; i < parties; i++ {
					p := new(Party)
					p.CKGProtocol = NewCKGProtocol(bfvContext)
					p.s = sk0_shards[i].Get()
					p.s1 = p.AllocateShares()
					ckgParties[i] = p
				}
				P0 := ckgParties[0]

				// Each party creates a new CKGProtocol instance
				for i, p := range ckgParties {
					p.GenShare(p.s, crp, p.s1)
					if i > 0 {
						P0.AggregateShares(p.s1, P0.s1, P0.s1)
					}
				}

				pk := &bfv.PublicKey{}
				P0.GenPublicKey(P0.s1, crp, pk)

				// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
				encryptorTest, err := bfvContext.NewEncryptorFromPk(pk)
				if err != nil {
					t.Error(err)
				}

				ciphertextTest, err := encryptorTest.EncryptNew(plaintextWant)

				if err != nil {
					t.Error(err)
				}

				if equalslice(coeffsWant.Coeffs[0], encoder.DecodeUint(decryptor_sk0.DecryptNew(ciphertextTest))) != true {
					t.Errorf("error : ckg protocol, cpk encrypt/decrypt test")
				}

			})

			t.Run(fmt.Sprintf("N=%d/logQ=%d/CKS", context.N, 60), func(t *testing.T) {

				type Party struct {
					*CKSProtocol
					s0    *ring.Poly
					s1    *ring.Poly
					share CKSShare
				}

				cksParties := make([]*Party, parties)
				for i := 0; i < parties; i++ {
					p := new(Party)
					p.CKSProtocol = NewCKSProtocol(bfvContext, 6.36)
					p.s0 = sk0_shards[i].Get()
					p.s1 = sk1_shards[i].Get()
					p.share = p.AllocateShare()
					cksParties[i] = p
				}
				P0 := cksParties[0]

				ciphertext, err := encryptor_pk0.EncryptNew(plaintextWant)
				if err != nil {
					t.Error(err)
				}

				// Each party creates its CKSProtocol instance with tmp = si-si'
				for i, p := range cksParties {
					p.GenShare(p.s0, p.s1, ciphertext, p.share)
					if i > 0 {
						P0.AggregateShares(p.share, P0.share, P0.share)
					}
				}

				ksCiphertext := bfvContext.NewCiphertext(1)
				P0.KeySwitch(P0.share, ciphertext, ksCiphertext)
				if equalslice(coeffsWant.Coeffs[0], encoder.DecodeUint(decryptor_sk1.DecryptNew(ksCiphertext))) != true {
					t.Errorf("error : CKS")
				}

				P0.KeySwitch(P0.share, ciphertext, ciphertext)

				if equalslice(coeffsWant.Coeffs[0], encoder.DecodeUint(decryptor_sk1.DecryptNew(ciphertext))) != true {
					t.Errorf("error : CKS in place")
				}

			})

			t.Run(fmt.Sprintf("N=%d/logQ=%d/PCKS", context.N, context.ModulusBigint.Value.BitLen()), func(t *testing.T) {

				type Party struct {
					*PCKSProtocol
					s     *ring.Poly
					share PCKSShare
				}

				pcksParties := make([]*Party, parties)
				for i := 0; i < parties; i++ {
					p := new(Party)
					p.PCKSProtocol = NewPCKSProtocol(bfvContext, 6.36)
					p.s = sk0_shards[i].Get()
					p.share = p.AllocateShares()
					pcksParties[i] = p
				}
				P0 := pcksParties[0]

				ciphertext, err := encryptor_pk0.EncryptNew(plaintextWant)
				ciphertextSwitched := bfvContext.NewCiphertext(1)
				if err != nil {
					t.Error(err)
				}

				for i, p := range pcksParties {
					p.GenShare(p.s, pk1, ciphertext, p.share)
					if i > 0 {
						P0.AggregateShares(p.share, P0.share, P0.share)
					}
				}

				P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)
				if equalslice(coeffsWant.Coeffs[0], encoder.DecodeUint(decryptor_sk1.DecryptNew(ciphertextSwitched))) != true {
					t.Errorf("error : PCKS")
				}
			})

			t.Run(fmt.Sprintf("N=%d/Qi=%dx%d/BOOT", context.N, len(context.Modulus), 60), func(t *testing.T) {

				// We store the plaintext coeffs as bigint for a reference to quantify the error
				coeffs_plaintext_bigint_fresh := make([]*ring.Int, bfvContext.N())
				bfvContext.ContextQ().PolyToBigint(plaintextWant.Value()[0], coeffs_plaintext_bigint_fresh)

				// We encrypt the plaintext
				ciphertext, err := encryptor_pk0.EncryptNew(plaintextWant)

				// ===== Boot instance =====
				crp := crpGenerators[0].Clock()

				refreshShares := make([]*RefreshShares, parties)

				for i := 0; i < parties; i++ {
					refreshShares[i] = GenRefreshShares(sk0_shards[i], ciphertext, bfvContext, crp, encoder)
				}
				// =========================

				// ==== We simulated added error of size Q/(T^2) and addit to the fresh ciphertext ====
				coeffs_bigint := make([]*ring.Int, bfvContext.N())
				bfvContext.ContextQ().PolyToBigint(ciphertext.Value()[0], coeffs_bigint)

				error_range := ring.Copy(bfvContext.ContextQ().ModulusBigint)
				error_range.Div(error_range, bfvContext.ContextT().ModulusBigint)
				error_range.Div(error_range, bfvContext.ContextT().ModulusBigint)

				for i := uint64(0); i < bfvContext.N(); i++ {
					coeffs_bigint[i].Add(coeffs_bigint[i], ring.RandInt(error_range))
				}

				bfvContext.ContextQ().SetCoefficientsBigint(coeffs_bigint, ciphertext.Value()[0])

				plaintextHave := decryptor_sk0.DecryptNew(ciphertext)
				coeffs_plaintext_bigint_error := make([]*ring.Int, bfvContext.N())
				bfvContext.ContextQ().PolyToBigint(plaintextHave.Value()[0], coeffs_plaintext_bigint_error)

				average_simulated_error := new(ring.Int)
				for i := uint64(0); i < bfvContext.N(); i++ {
					average_simulated_error.Add(average_simulated_error, coeffs_plaintext_bigint_fresh[i])
					average_simulated_error.Sub(average_simulated_error, coeffs_plaintext_bigint_error[i])
				}

				average_simulated_error.Value.Abs(&average_simulated_error.Value)
				average_simulated_error.Div(average_simulated_error, ring.NewUint(bfvContext.N()))
				// =======================================================================================

				// We boot the ciphertext with the simulated error
				Refresh(ciphertext, sk0.Get(), refreshShares, bfvContext, crp, encoder)

				// We decrypt and compare with the original plaintext
				plaintextHave = decryptor_sk0.DecryptNew(ciphertext)
				coeffs_plaintext_bigint_booted := make([]*ring.Int, bfvContext.N())
				bfvContext.ContextQ().PolyToBigint(plaintextHave.Value()[0], coeffs_plaintext_bigint_booted)

				average_residual_error := new(ring.Int)
				for i := uint64(0); i < bfvContext.N(); i++ {
					average_residual_error.Add(average_residual_error, coeffs_plaintext_bigint_fresh[i])
					average_residual_error.Sub(average_residual_error, coeffs_plaintext_bigint_booted[i])
				}

				average_residual_error.Value.Abs(&average_residual_error.Value)
				average_residual_error.Div(average_residual_error, ring.NewUint(bfvContext.N()))

				coeffsTest := encoder.DecodeUint(plaintextHave)
				if err != nil {
					log.Fatal(err)
				}

				t.Logf("Average simulated error before refresh (log2): %d", average_simulated_error.Value.BitLen())
				t.Logf("Average residual error after refresh (log2): %d", average_residual_error.Value.BitLen())

				if equalslice(coeffsWant.Coeffs[0], coeffsTest) != true {
					t.Errorf("error : BOOT")
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

func test_EKG_Protocol(bfvCtx *bfv.BfvContext, parties int, bitDecomp uint64, ekgProtocols []*RKGProtocol, sk []*bfv.SecretKey, ephemeralKeys []*ring.Poly, crp [][]*ring.Poly) *bfv.EvaluationKey {

	type Party struct {
		*RKGProtocol
		u      *ring.Poly
		s      *ring.Poly
		share1 RKGShareRoundOne
		share2 RKGShareRoundTwo
		share3 RKGShareRoundThree
	}

	rkgParties := make([]*Party, parties)

	for i := range rkgParties {
		p := new(Party)
		p.RKGProtocol = ekgProtocols[i]
		p.u = ephemeralKeys[i]
		p.s = sk[i].Get()
		p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()
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
		p.GenShareRoundThree(P0.share2, p.u, p.s, p.share3)
		if i > 0 {
			P0.AggregateShareRoundThree(p.share3, P0.share3, P0.share3)
		}
	}

	evk := bfvCtx.NewRelinKey(1, bitDecomp)
	P0.GenRelinearizationKey(P0.share2, P0.share3, evk)
	return evk
}

func Test_Marshalling(t *testing.T){
	//verify if the un.marshalling works properly
	log.Print("Verifying marshalling for Key Generation")
	bfvCtx,_ := bfv.NewBfvContextWithParam(&bfv.DefaultParams[0])
	KeyGenerator := bfvCtx.NewKeyGenerator()
	crsGen, _ := NewCRPGenerator([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'}, bfvCtx.ContextQ())
	sk := KeyGenerator.NewSecretKey()
	crs := crsGen.Clock()
	keygenProtocol := NewCKGProtocol(bfvCtx)
	KeyGenShareBefore := keygenProtocol.AllocateShares()
	keygenProtocol.GenShare(sk.Get(),crs, KeyGenShareBefore)
	//now we marshall it
	data , err := KeyGenShareBefore.MarshalBinary()

	if err != nil{
		log.Fatal("Could not marshal the CKGShare : ", err )
		t.Fail()
	}

	KeyGenShareAfter := new(CKGShare)
	err = KeyGenShareAfter.UnmarshalBinary(data)
	if err != nil{
		log.Fatal("Could not unmarshal the CKGShare : ",  err)
		t.Fail()
	}

	//comparing the results
	if KeyGenShareBefore.GetDegree() != KeyGenShareAfter.GetDegree(){
		log.Print("Unmatched degree on key gen shares")
		t.Fail()
	}

	for i := 0 ;i < KeyGenShareBefore.GetLenModuli();i++{
		if !equalslice(KeyGenShareAfter.Coeffs[i],KeyGenShareBefore.Coeffs[i]){
			log.Print("Non equal slices in CKGShare")
			t.Fail()
		}

	}


	log.Print("CKGShare marshalling ok ")

	//Check marshalling for the PCKS
	Ciphertext := bfvCtx.NewRandomCiphertext(1)
	KeySwitchProtocol := NewPCKSProtocol(bfvCtx,bfvCtx.Sigma())
	SwitchShare := KeySwitchProtocol.AllocateShares()
	pk := KeyGenerator.NewPublicKey(sk)
	KeySwitchProtocol.GenShare(sk.Get(),pk,Ciphertext,SwitchShare)


	data, err = SwitchShare.MarshalBinary()
	if err != nil{
		log.Print("Error on PCKSShare marshalling : " , err)
		t.Fail()
	}
	SwitchShareReceiver := new(PCKSShare)
	err = SwitchShareReceiver.UnmarshalBinary(data)
	if err != nil{
		log.Print("Error on PCKSShare unmarshalling : " , err)
	}



	for i := 0 ; i < 2; i ++{
		//compare the shares.
		ringBefore := SwitchShare.share[i]
		ringAfter := SwitchShareReceiver.share[i]
		if ringBefore.GetDegree() != ringAfter.GetDegree(){
			log.Print("Error on degree matching")
			t.Fail()
		}
		for d:= 0 ; d < ringAfter.GetLenModuli(); d++{
			if !equalslice(ringAfter.Coeffs[d],ringBefore.Coeffs[d]){
				log.Print("Non equal slices in PCKSShare")
				t.Fail()
			}
		}

	}

	log.Print("PCKSShare marshalling ok ")

	//Now for CKSShare ~ its similar to PKSShare
	//todo write test for cksshare..











}

func Test_Relin_Marshalling(t *testing.T){
	bfvCtx ,_:= bfv.NewBfvContextWithParam(&bfv.DefaultParams[0])
	modulus := bfvCtx.ContextQ().Modulus
	bitDecomp := 60
	bitLog := uint64((60 + (60 % bitDecomp)) / bitDecomp)
	var err error
	crpGenerator,_ :=  NewCRPGenerator(nil, bfvCtx.ContextQ())


	crp := make([][]*ring.Poly, len(modulus))
	for j := 0; j < len(modulus); j++ {
		crp[j] = make([]*ring.Poly, bitLog)
		for u := uint64(0); u < bitLog; u++ {
			crp[j][u] = crpGenerator.Clock()
		}
	}

	rlk := NewEkgProtocol(bfvCtx,uint64(bitDecomp))
	u,_ := rlk.NewEphemeralKey(1/3.0)
	sk := bfvCtx.NewKeyGenerator().NewSecretKey()
	log.Print("Starting to test marshalling for share one")

	r1,r2,r3 := rlk.AllocateShares()
	rlk.GenShareRoundOne(u,sk.Get(),crp,r1)
	data , err := r1.MarshalBinary()
	if err != nil{
		log.Print("Error in marshalling round 1 key : " , err)
		t.Fail()
	}

	r1After := new(RKGShareRoundOne)
	err = r1After.UnmarshalBinary(data)
	if err != nil{
		log.Print("Error in unmarshalling round 1 key : " , err)
		t.Fail()
	}

	log.Print("Now comparing keys for round 1 ")

	if r1.bitLog != r1After.bitLog || r1.modulus != r1After.modulus{
		log.Print("Error bitlog or modulus are different ")
		t.Fail()
	}


	for i := 0 ; i < int(r1.modulus) ; i ++{
		for j := 0 ; j < int(r1.bitLog); j++{
			a := r1.share[i][j]
			b := r1After.share[i][j]
			for k := 0 ; k < a.GetLenModuli(); k++{
				if !equalslice(a.Coeffs[k],b.Coeffs[k]){
					log.Print("Error : coeffs of rings do not match")
					t.Fail()
				}
			}
		}
	}

	log.Print("Sucess : relin key round 1 ok ")

	log.Print("Starting to test marshalling for share two")
	rlk.GenShareRoundTwo(r1,sk.Get(),crp,r2)

	data , err = r2.MarshalBinary()
	if err != nil{
		log.Print("Error on marshalling relin key round 2 : " , err)
		t.Fail()
	}

	r2After := new(RKGShareRoundTwo)
	err = r2After.UnmarshalBinary(data)
	if err != nil{
		log.Print("Error on unmarshalling relin key round 2 : ", err)
		t.Fail()
	}

	log.Print("Now comparing keys for round 2 ")

	if r2.bitLog != r2After.bitLog || r2.modulus != r2After.modulus{
		log.Print("Error bitlog or modulus are different ")
		t.Fail()
	}


	for i := 0 ; i < int(r2.modulus) ; i ++{
		for j := 0 ; j < int(r2.bitLog); j++{
			for idx := 0 ; idx < 2 ; idx++ {
				a := r2.share[i][j][idx]
				b := r2After.share[i][j][idx]
				for k := 0; k < a.GetLenModuli(); k++ {
					if !equalslice(a.Coeffs[k], b.Coeffs[k]) {
						log.Print("Error : coeffs of rings do not match")
						t.Fail()
					}
				}
			}
		}
	}

	log.Print("Success : reling key round 2 ok ")


	log.Print("Starting to test marshalling for share three")

	rlk.GenShareRoundThree(r2,u,sk.Get(),r3)

	data , err = r3.MarshalBinary()
	if err != nil{
		log.Print("Error in marshalling round 3 key : " , err)
		t.Fail()
	}

	r3After := new(RKGShareRoundThree)
	err = r3After.UnmarshalBinary(data)
	if err != nil{
		log.Print("Error in unmarshalling round 3 key : " , err)
		t.Fail()
	}

	log.Print("Now comparing keys for round 3 ")

	if r3.bitLog != r3After.bitLog || r3.modulus != r3After.modulus{
		log.Print("Error bitlog or modulus are different ")
		t.Fail()
	}


	for i := 0 ; i < int(r3.modulus) ; i ++{
		for j := 0 ; j < int(r3.bitLog); j++{
			a := r3.share[i][j]
			b := r3After.share[i][j]
			for k := 0 ; k < a.GetLenModuli(); k++{
				if !equalslice(a.Coeffs[k],b.Coeffs[k]){
					log.Print("Error : coeffs of rings do not match")
					t.Fail()
				}
			}
		}
	}

	log.Print("Success : relin key for round 3 ok ")


	log.Print("All marshalling is passed for relin keys")









}