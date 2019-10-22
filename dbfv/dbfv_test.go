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

type dbfvparams struct {
	parties int

	bfvcontext        *bfv.BfvContext
	contextKeys           *ring.Context
	contextCiphertext *ring.Context
	contextT          *ring.Context
	encoder           *bfv.BatchEncoder
	kgen              *bfv.KeyGenerator

	sk0_shards []*bfv.SecretKey
	sk0        *bfv.SecretKey

	sk1        *bfv.SecretKey
	sk1_shards []*bfv.SecretKey

	pk0 *bfv.PublicKey
	pk1 *bfv.PublicKey

	encryptor_pk0 *bfv.Encryptor
	encryptor_pk1 *bfv.Encryptor
	decryptor_sk0 *bfv.Decryptor
	decryptor_sk1 *bfv.Decryptor
	evaluator     *bfv.Evaluator
}

func Test_DBFVScheme(t *testing.T) {

	var err error

	paramSets := bfv.DefaultParams[0:1]
	parties := 5

	//sigmaSmudging := 6.36

	for _, params := range paramSets {

		dbfvParams := new(dbfvparams)

		dbfvParams.parties = parties

		// nParties data indpendant element
		dbfvParams.bfvcontext = bfv.NewBfvContext()
		if err := dbfvParams.bfvcontext.SetParameters(&params); err != nil {
			log.Fatal(err)
		}

		dbfvParams.kgen = dbfvParams.bfvcontext.NewKeyGenerator()

		dbfvParams.evaluator = dbfvParams.bfvcontext.NewEvaluator()

		dbfvParams.contextKeys = dbfvParams.bfvcontext.ContextKeys()

		dbfvParams.contextCiphertext = dbfvParams.bfvcontext.ContextQ()

		dbfvParams.contextT = dbfvParams.bfvcontext.ContextT()

		dbfvParams.encoder, err = dbfvParams.bfvcontext.NewBatchEncoder()
		check(t, err)



		// SecretKeys
		dbfvParams.sk0_shards = make([]*bfv.SecretKey, dbfvParams.parties)
		dbfvParams.sk1_shards = make([]*bfv.SecretKey, dbfvParams.parties)
		tmp0 := dbfvParams.contextKeys.NewPoly()
		tmp1 := dbfvParams.contextKeys.NewPoly()

		for i := 0; i < dbfvParams.parties; i++ {
			dbfvParams.sk0_shards[i] = dbfvParams.kgen.NewSecretKey()
			dbfvParams.sk1_shards[i] = dbfvParams.kgen.NewSecretKey()
			dbfvParams.contextKeys.Add(tmp0, dbfvParams.sk0_shards[i].Get(), tmp0)
			dbfvParams.contextKeys.Add(tmp1, dbfvParams.sk1_shards[i].Get(), tmp1)
		}

		dbfvParams.sk0 = new(bfv.SecretKey)
		dbfvParams.sk1 = new(bfv.SecretKey)

		dbfvParams.sk0.Set(tmp0)
		dbfvParams.sk1.Set(tmp1)

		// Publickeys
		dbfvParams.pk0 = dbfvParams.kgen.NewPublicKey(dbfvParams.sk0)
		dbfvParams.pk1 = dbfvParams.kgen.NewPublicKey(dbfvParams.sk1);

		// Encryptors
		dbfvParams.encryptor_pk0, err= dbfvParams.bfvcontext.NewEncryptorFromPk(dbfvParams.pk0)
		check(t, err)

		//encryptor_pk1, err := bfvContext.NewEncryptor(pk1)
		//if err != nil {
		//	log.Fatal(err)
		//}

		// Decryptors
		if dbfvParams.decryptor_sk0, err = dbfvParams.bfvcontext.NewDecryptor(dbfvParams.sk0); err != nil {
			log.Fatal(err)
		}

		if dbfvParams.decryptor_sk1, err = dbfvParams.bfvcontext.NewDecryptor(dbfvParams.sk1); err != nil {
			log.Fatal(err)
		}

		test_EKG(dbfvParams, t)
		test_PRNG(dbfvParams, t)
		test_CKG(dbfvParams, t)
		
		test_EKGNaive(dbfvParams, t)
		test_RKG(dbfvParams, t)
		test_CKS(dbfvParams, t)
		test_PCKS(dbfvParams, t)
		test_BOOT(dbfvParams, t)
	}
}

			

func test_PRNG(params *dbfvparams, t *testing.T) {

	contextKeys := params.contextKeys
	bfvContext := params.bfvcontext


	t.Run(fmt.Sprintf("N=%d/logQ=%d/CRS_PRNG", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

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

		crs_generator_1, _ := NewCRPGenerator(nil, contextKeys)
		crs_generator_2, _ := NewCRPGenerator(nil, contextKeys)

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
}



func test_EKGNaive(params *dbfvparams, t *testing.T) {

	var err error 

	contextKeys := params.contextKeys
	parties := params.parties
	bfvContext := params.bfvcontext
	sk0_shards := params.sk0_shards
	pk0 := params.pk0
	evaluator := params.evaluator

	t.Run(fmt.Sprintf("N=%d/logQ=%d/EKG_Naive", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

		// Each party instantiate an ekg naive protocole
		ekgNaive := make([]*EkgProtocolNaive, parties)
		for i := 0; i < parties; i++ {
			ekgNaive[i] = NewEkgProtocolNaive(bfvContext)
		}

		evk := test_EKG_Protocol_Naive(parties, sk0_shards, pk0, ekgNaive)

		rlk := new(bfv.EvaluationKey)
		rlk.SetRelinKeys([][][2]*ring.Poly{evk[0]})


		coeffs, _, ciphertext, _ := newTestVectors(params)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
			coeffs[i] %= params.contextT.Modulus[0]
		}

		ciphertextMul := bfvContext.NewCiphertext(ciphertext.Degree() * 2)
		err = evaluator.Mul(ciphertext, ciphertext, ciphertextMul)
		check(t, err)

		res := bfvContext.NewCiphertext(1)
		err := evaluator.Relinearize(ciphertextMul, rlk, res)
		check(t, err)

		verifyTestVectors(params.decryptor_sk0, params.encoder, coeffs, ciphertextMul, t)
	})
}

func test_EKG(params *dbfvparams, t *testing.T) {

	var err error
	bfvContext := params.bfvcontext
	contextKeys := params.contextKeys
	parties := params.parties
	sk0_shards := params.sk0_shards
	evaluator := params.evaluator

	t.Run(fmt.Sprintf("N=%d/logQ=%d/EKG", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

		// Each party instantiate an ekg naive protocole
		ekg := make([]*RKGProtocol, parties)
		ephemeralKeys := make([]*ring.Poly, parties)

		crpGenerators := make([]*CRPGenerator, parties)
		crp := make([]*ring.Poly, bfvContext.Beta())
		for j := uint64(0); j < bfvContext.Beta(); j++ {
			crpGenerators[j], err = NewCRPGenerator(nil, contextKeys)

			check(t, err)

			crpGenerators[j].Seed([]byte{})
			crp[j] = crpGenerators[j].Clock()
		}

		for i := 0; i < parties; i++ {
			ekg[i] = NewEkgProtocol(bfvContext)
			ephemeralKeys[i], _ = ekg[i].NewEphemeralKey(1.0 / 3)
		}

		rlk := test_EKG_Protocol(bfvContext, parties, ekg, sk0_shards, ephemeralKeys, crp)

		coeffs, _, ciphertext, _ := newTestVectors(params)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
			coeffs[i] %= params.contextT.Modulus[0]
		}

		ciphertextMul := bfvContext.NewCiphertext(ciphertext.Degree() * 2)
		if err = evaluator.Mul(ciphertext, ciphertext, ciphertextMul); err != nil {
			log.Fatal(err)
		}


		res := bfvContext.NewCiphertext(1)
		if err := evaluator.Relinearize(ciphertextMul, rlk, res); err != nil {
			t.Error(err)
		}

		verifyTestVectors(params.decryptor_sk0, params.encoder, coeffs, ciphertextMul, t)

	})
}

func test_RKG(params *dbfvparams, t *testing.T) {

	var err error

	contextKeys := params.contextKeys
	parties := params.parties
	bfvContext := params.bfvcontext
	keygenerator := params.kgen
	sk0_shards := params.sk0_shards
	evaluator := params.evaluator

	rkg := make([]*RKG, parties)
	for i := 0; i < parties; i++ {
		rkg[i] = NewRKG(bfvContext)
	}

	crpGenerator, err := NewCRPGenerator(nil, contextKeys)
	check(t, err)
	crpGenerator.Seed([]byte{})

	crp := make([]*ring.Poly, bfvContext.Beta())
	for i := uint64(0); i < bfvContext.Beta(); i++ {
		crp[i] = crpGenerator.Clock()
	}

	t.Run(fmt.Sprintf("N=%d/logQ=%d/RKG_rot_rows", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {
		

		coeffs, _, ciphertext, _ := newTestVectors(params)

		shares := make([][]*ring.Poly, parties)
		for i := 0; i < parties; i++ {
			shares[i] = rkg[i].GenShareRotRow(sk0_shards[i].Get(), crp)
		}

		rkg[0].AggregateRotRow(shares, crp)
		rotkey := rkg[0].Finalize(keygenerator)

		receiver := bfvContext.NewCiphertext(1)
		if err = evaluator.RotateRows(ciphertext, rotkey, receiver); err != nil {
			log.Fatal(err)
		}

		coeffs = append(coeffs[contextKeys.N>>1:], coeffs[:contextKeys.N>>1]...)

		verifyTestVectors(params.decryptor_sk0, params.encoder, coeffs, receiver, t)

	})

	t.Run(fmt.Sprintf("N=%d/logQ=%d/RKG_rot_col_pow2", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

		coeffs, _, ciphertext, _ := newTestVectors(params)
		mask := (contextKeys.N >> 1) - 1

		receiver := bfvContext.NewCiphertext(ciphertext.Degree())
		for n := uint64(1); n < contextKeys.N>>1; n <<= 1 {

			shares := make([][]*ring.Poly, parties)
			for i := 0; i < parties; i++ {
				shares[i] = rkg[i].GenShareRotLeft(sk0_shards[i].Get(), n, crp)
			}

			rkg[0].AggregateRotColL(shares, n, crp)
			rotkey := rkg[0].Finalize(keygenerator)

			if err = evaluator.RotateColumns(ciphertext, n, rotkey, receiver); err != nil {
				log.Fatal(err)
			}

			coeffsWant := make([]uint64, contextKeys.N)

			for i := uint64(0); i < contextKeys.N>>1; i++ {
				coeffsWant[i] = coeffs[(i+n)&mask]
				coeffsWant[i+(contextKeys.N>>1)] = coeffs[((i+n)&mask)+(contextKeys.N>>1)]
			}

			verifyTestVectors(params.decryptor_sk0, params.encoder, coeffsWant, receiver, t)
		}

	})
}


func test_CKG(params *dbfvparams, t *testing.T) {

	bfvContext := params.bfvcontext
	contextKeys := params.contextKeys
	parties := params.parties
	sk0_shards := params.sk0_shards

	t.Run(fmt.Sprintf("N=%d/logQ=%d/CKG", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

		crpGenerator, err := NewCRPGenerator(nil, contextKeys)
		check(t, err)
		crpGenerator.Seed([]byte{})
		crp := crpGenerator.Clock()

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

		coeffs, plaintextWant, _, _ := newTestVectors(params)

		ciphertextTest, err := encryptorTest.EncryptNew(plaintextWant)

		if err != nil {
			t.Error(err)
		}

		verifyTestVectors(params.decryptor_sk0, params.encoder, coeffs, ciphertextTest, t)
	})
}



func test_CKS(params *dbfvparams, t *testing.T) {

	bfvContext := params.bfvcontext
	contextKeys := params.contextKeys
	parties := params.parties
	sk0_shards := params.sk0_shards
	sk1_shards := params.sk1_shards
	encoder := params.encoder
	decryptor_sk1 := params.decryptor_sk1


	t.Run(fmt.Sprintf("N=%d/logQ=%d/CKS", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

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

		coeffs, _, ciphertext, _ := newTestVectors(params)


		// Each party creates its CKSProtocol instance with tmp = si-si'
		for i, p := range cksParties {
			p.GenShare(p.s0, p.s1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		

		ksCiphertext := bfvContext.NewCiphertext(1)
		P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

		if equalslice(coeffs, encoder.DecodeUint(decryptor_sk1.DecryptNew(ksCiphertext))) != true {
			t.Errorf("error : CKS")
		}

		P0.KeySwitch(P0.share, ciphertext, ciphertext)

		if equalslice(coeffs, encoder.DecodeUint(decryptor_sk1.DecryptNew(ciphertext))) != true {
			t.Errorf("error : CKS in place")
		}

	})
}



func test_PCKS(params *dbfvparams, t *testing.T) {

	bfvContext := params.bfvcontext
	contextKeys := params.contextKeys
	parties := params.parties
	sk0_shards := params.sk0_shards
	pk1 := params.pk1
	encoder := params.encoder
	decryptor_sk1 := params.decryptor_sk1

	t.Run(fmt.Sprintf("N=%d/logQ=%d/PCKS", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

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

		coeffs, _, ciphertext, _ := newTestVectors(params)

		ciphertextSwitched := bfvContext.NewCiphertext(1)


		for i, p := range pcksParties {
			p.GenShare(p.s, pk1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)
		if equalslice(coeffs, encoder.DecodeUint(decryptor_sk1.DecryptNew(ciphertextSwitched))) != true {
			t.Errorf("error : PCKS")
		}
	})
}


func test_BOOT(params *dbfvparams, t *testing.T) {

	bfvContext := params.bfvcontext
	contextKeys := params.contextKeys
	parties := params.parties
	sk0_shards := params.sk0_shards
	encoder := params.encoder
	decryptor_sk0 := params.decryptor_sk0
	sk0 := params.sk0


	t.Run(fmt.Sprintf("N=%d/logQ=%d/BOOT", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

		crpGenerator, err := NewCRPGenerator(nil, contextKeys)
		check(t, err)
		crpGenerator.Seed([]byte{})
		crp := crpGenerator.Clock()

		coeffs, plaintextWant, ciphertext, _ := newTestVectors(params)

		// We store the plaintext coeffs as bigint for a reference to quantify the error
		coeffs_plaintext_bigint_fresh := make([]*ring.Int, bfvContext.N())
		bfvContext.ContextQ().PolyToBigint(plaintextWant.Value()[0], coeffs_plaintext_bigint_fresh)

		// We encrypt the plaintext

		// ===== Boot instance =====

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

		if equalslice(coeffs, coeffsTest) != true {
			t.Errorf("error : BOOT")
		}
	})
}

func test_EKG_Protocol_Naive(parties int, sk []*bfv.SecretKey, collectivePk *bfv.PublicKey, ekgNaive []*EkgProtocolNaive) [][][2]*ring.Poly {

	// ROUND 0
	// Each party generates its samples
	samples := make([][][2]*ring.Poly, parties)
	for i := 0; i < parties; i++ {
		samples[i] = ekgNaive[i].GenSamples(sk[i].Get(), collectivePk.Get())
	}

	// ROUND 1
	// Each party aggretates its sample with the other n-1 samples
	aggregatedSamples := make([][][2]*ring.Poly, parties)
	for i := 0; i < parties; i++ {
		aggregatedSamples[i] = ekgNaive[i].Aggregate(sk[i].Get(), collectivePk.Get(), samples)
	}

	// ROUND 2
	// Each party aggregates sums its aggregatedSample with the other n-1 aggregated samples
	evk := make([][][2]*ring.Poly, parties)
	for i := 0; i < parties; i++ {
		evk[i] = ekgNaive[i].Finalize(aggregatedSamples)
	}

	return evk
}

func test_EKG_Protocol(bfvCtx *bfv.BfvContext, parties int, ekgProtocols []*RKGProtocol, sk []*bfv.SecretKey, ephemeralKeys []*ring.Poly, crp []*ring.Poly) *bfv.EvaluationKey {

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

	evk := bfvCtx.NewRelinKeyEmpty(1)
	P0.GenRelinearizationKey(P0.share2, P0.share3, evk)
	return evk
}

func newTestVectors(params *dbfvparams) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext, err error) {

	coeffsPol := params.contextT.NewUniformPoly()

	plaintext = params.bfvcontext.NewPlaintext()

	if err = params.encoder.EncodeUint(coeffsPol.Coeffs[0], plaintext); err != nil {
		return nil, nil, nil, err
	}

	if ciphertext, err = params.encryptor_pk0.EncryptNew(plaintext); err != nil {
		return nil, nil, nil, err
	}

	return coeffsPol.Coeffs[0], plaintext, ciphertext, nil
}

func verifyTestVectors(decryptor *bfv.Decryptor, encoder *bfv.BatchEncoder, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {

	if bfv.EqualSlice(coeffs, encoder.DecodeUint(decryptor.DecryptNew(ciphertext))) != true {
		t.Errorf("decryption error")
	}
}
