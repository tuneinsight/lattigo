package dckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"log"
	"math"
	"testing"
)

func Test_DCKKScheme(t *testing.T) {

	//sigmaSmudging := 6.36

	var err error

	var parties uint64
	parties = 5

	params := ckks.DefaultParams[14]

	alpha := uint64(len(params.P))
	beta := uint64(math.Ceil(float64(len(params.Modulichain)) / float64(alpha)))

	_ = beta

	var ckkscontext *ckks.CkksContext

	if ckkscontext, err = ckks.NewCkksContext(params); err != nil {
		log.Fatal(err)
	}

	log.Printf("Generating CkksContext for logN=%d/logQ=%d/levels=%d/scale=%f/sigma=%f", ckkscontext.LogN(), ckkscontext.LogQ(), ckkscontext.Levels(), ckkscontext.Scale(), ckkscontext.Sigma())

	logN := ckkscontext.LogN()
	logQ := ckkscontext.LogQ()
	levels := ckkscontext.Levels()
	scale := ckkscontext.Scale()

	encoder := ckkscontext.NewEncoder()

	kgen := ckkscontext.NewKeyGenerator()

	evaluator := ckkscontext.NewEvaluator()
	_ = evaluator

	context := ckkscontext.ContextKeys()
	contextCiphertexts := ckkscontext.ContextQ()
	_ = contextCiphertexts

	coeffsWant := make([]complex128, ckkscontext.Slots())
	for i := uint64(0); i < ckkscontext.Slots(); i++ {
		coeffsWant[i] = randomComplex(1)
	}

	plaintextWant := ckkscontext.NewPlaintext(levels-1, scale)
	if err = encoder.Encode(plaintextWant, coeffsWant, ckkscontext.Slots()); err != nil {
		log.Fatal(err)
	}

	crpGenerators := make([]*CRPGenerator, parties)
	for i := uint64(0); i < parties; i++ {
		crpGenerators[i], err = NewCRPGenerator(nil, context)
		if err != nil {
			log.Fatal(err)
		}
		crpGenerators[i].Seed([]byte{})
	}

	// SecretKeys
	sk0_shards := make([]*ckks.SecretKey, parties)
	sk1_shards := make([]*ckks.SecretKey, parties)
	tmp0 := context.NewPoly()
	tmp1 := context.NewPoly()

	for i := uint64(0); i < parties; i++ {
		sk0_shards[i] = kgen.NewSecretKey()
		sk1_shards[i] = kgen.NewSecretKey()
		context.Add(tmp0, sk0_shards[i].Get(), tmp0)
		context.Add(tmp1, sk1_shards[i].Get(), tmp1)
	}

	sk0 := new(ckks.SecretKey)
	sk1 := new(ckks.SecretKey)

	sk0.Set(tmp0)
	sk1.Set(tmp1)

	// Publickeys
	pk0 := kgen.NewPublicKey(sk0)
	pk1 := kgen.NewPublicKey(sk1)

	_ = pk1

	// Encryptors
	encryptor_pk0, err := ckkscontext.NewEncryptorFromPk(pk0)
	_ = encryptor_pk0
	if err != nil {
		log.Fatal(err)
	}

	// Decryptors
	decryptor_sk0, err := ckkscontext.NewDecryptor(sk0)
	_ = decryptor_sk0
	if err != nil {
		log.Fatal(err)
	}

	decryptor_sk1, err := ckkscontext.NewDecryptor(sk1)
	if err != nil {
		log.Fatal(err)
	}

	_ = decryptor_sk1

	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/CRS_PRNG", parties, logN, logQ, levels, scale), func(t *testing.T) {

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

		if ckkscontext.ContextKeys().Equal(p0, p1) != true {
			t.Errorf("error : crs prng generator")
		}
	})

	// EKG (OK multiple levels)
	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/EKG", parties, logN, logQ, levels, scale), func(t *testing.T) {

		// Each party instantiate an ekg naive protocole
		ekg := make([]*EkgProtocol, parties)
		ephemeralKeys := make([]*ring.Poly, parties)
		crp := make([][]*ring.Poly, parties)

		for i := uint64(0); i < parties; i++ {

			ekg[i] = NewEkgProtocol(context, ckkscontext.KeySwitchPrimes())
			ephemeralKeys[i], _ = ekg[i].NewEphemeralKey(1.0 / 3.0)
			crp[i] = make([]*ring.Poly, beta)

			for j := uint64(0); j < beta; j++ {
				crp[i][j] = crpGenerators[i].Clock()

			}
		}

		evk := test_EKG_Protocol(parties, ekg, sk0_shards, ephemeralKeys, crp)

		rlk, err := kgen.SetRelinKeys(evk[0])
		if err != nil {
			log.Fatal(err)
		}

		coeffs, _, ciphertext, _ := new_test_vectors(ckkscontext, encoder, encryptor_pk0, 1)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
		}

		if err := evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext); err != nil {
			log.Fatal(err)
		}

		evaluator.Rescale(ciphertext, ckkscontext.Scale(), ciphertext)

		if ciphertext.Degree() != 1 {
			t.Errorf("EKG_NAIVE -> bad relinearize")
		}

		verify_test_vectors(ckkscontext, encoder, decryptor_sk0, coeffs, ciphertext, t)

	})

	// EKG Naive (OK multiple levels)
	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/EKG_NAIVE", parties, logN, logQ, levels, scale), func(t *testing.T) {

		// Each party instantiate an ekg naive protocole
		ekgNaive := make([]*EkgProtocolNaive, parties)
		for i := uint64(0); i < parties; i++ {
			ekgNaive[i] = NewEkgProtocolNaive(context, ckkscontext.KeySwitchPrimes())
		}

		evk := test_EKG_Protocol_Naive(parties, sk0_shards, pk0, ekgNaive)

		rlk, err := kgen.SetRelinKeys(evk[0])
		if err != nil {
			log.Fatal(err)
		}

		coeffs, _, ciphertext, _ := new_test_vectors(ckkscontext, encoder, encryptor_pk0, 1)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
		}

		if err := evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext); err != nil {
			log.Fatal(err)
		}

		evaluator.Rescale(ciphertext, ckkscontext.Scale(), ciphertext)

		if ciphertext.Degree() != 1 {
			t.Errorf("EKG_NAIVE -> bad relinearize")
		}

		verify_test_vectors(ckkscontext, encoder, decryptor_sk0, coeffs, ciphertext, t)
	})

	// RKG conjugate (OK multiple levels)
	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/RKG_conjugate", parties, logN, logQ, levels, scale), func(t *testing.T) {

		rkg := make([]*RKG, parties)
		crp := make([][]*ring.Poly, parties)

		for i := uint64(0); i < parties; i++ {

			rkg[i] = NewRKG(context, ckkscontext.KeySwitchPrimes())
			crp[i] = make([]*ring.Poly, len(context.Modulus))

			for j := 0; j < len(context.Modulus); j++ {
				crp[i][j] = crpGenerators[i].Clock()
			}
		}

		coeffs, _, ciphertext, _ := new_test_vectors(ckkscontext, encoder, encryptor_pk0, 1)

		shares := make([][]*ring.Poly, parties)
		for i := uint64(0); i < parties; i++ {
			shares[i] = rkg[i].GenShareRotRow(sk0_shards[i].Get(), crp[i])
		}

		rkg[0].AggregateRotRow(shares, crp[0])
		rotkey := rkg[0].Finalize(kgen)

		if err = evaluator.Conjugate(ciphertext, rotkey, ciphertext); err != nil {
			log.Fatal(err)
		}

		for i := range coeffs {
			coeffs[i] = complex(real(coeffs[i]), -imag(coeffs[i]))
		}

		verify_test_vectors(ckkscontext, encoder, decryptor_sk0, coeffs, ciphertext, t)

	})

	// RKG rot col pow2 (OK multiple levels)
	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/RKG_rot_col_pow2", parties, logN, logQ, levels, scale), func(t *testing.T) {

		rkg := make([]*RKG, parties)
		crp := make([][]*ring.Poly, parties)

		for i := uint64(0); i < parties; i++ {

			rkg[i] = NewRKG(context, ckkscontext.KeySwitchPrimes())
			crp[i] = make([]*ring.Poly, len(context.Modulus))

			for j := 0; j < len(context.Modulus); j++ {
				crp[i][j] = crpGenerators[i].Clock()
			}
		}

		coeffs, _, ciphertext, _ := new_test_vectors(ckkscontext, encoder, encryptor_pk0, 1)
		mask := ckkscontext.Slots() - 1

		receiver := ckkscontext.NewCiphertext(ciphertext.Degree(), ciphertext.Level(), ciphertext.Scale())
		for n := uint64(1); n < context.N>>1; n <<= 1 {

			shares := make([][]*ring.Poly, parties)
			for i := uint64(0); i < parties; i++ {
				shares[i] = rkg[i].GenShareRotLeft(sk0_shards[i].Get(), n, crp[i])
			}

			rkg[0].AggregateRotColL(shares, n, crp[0])
			rotkey := rkg[0].Finalize(kgen)

			if err = evaluator.RotateColumns(ciphertext, n, rotkey, receiver); err != nil {
				log.Fatal(err)
			}

			coeffsWant := make([]complex128, ckkscontext.Slots())

			for i := uint64(0); i < ckkscontext.Slots(); i++ {
				coeffsWant[i] = coeffs[(i+n)&mask]
			}

			verify_test_vectors(ckkscontext, encoder, decryptor_sk0, coeffsWant, receiver, t)
		}

	})

	// RKG rot col pow2 (OK multiple levels (always assumed to be at max level))
	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/CKG", parties, logN, logQ, levels, scale), func(t *testing.T) {

		crp := make([]*ring.Poly, parties)
		for i := uint64(0); i < parties; i++ {
			crp[i] = crpGenerators[i].Clock()
		}

		ckg := make([]*CKG, parties)
		for i := uint64(0); i < parties; i++ {
			ckg[i] = NewCKG(context, crp[i])
		}

		// Each party creates a new CKG instance
		shares := make([]*ring.Poly, parties)
		for i := uint64(0); i < parties; i++ {
			ckg[i].GenShare(sk0_shards[i].Get())
			shares[i] = ckg[i].GetShare()
		}

		pkTest := make([]*ckks.PublicKey, parties)
		for i := uint64(0); i < parties; i++ {
			ckg[i].AggregateShares(shares)
			pkTest[i], err = ckg[i].Finalize()
			if err != nil {
				log.Fatal(err)
			}
		}

		// Verifies that all parties have the same share collective public key
		for i := uint64(1); i < parties; i++ {
			if context.Equal(pkTest[0].Get()[0], pkTest[i].Get()[0]) != true || ckkscontext.ContextKeys().Equal(pkTest[0].Get()[1], pkTest[i].Get()[1]) != true {
				t.Errorf("error : ckg protocol, cpk establishement")
			}
		}

		// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
		encryptorTest, err := ckkscontext.NewEncryptorFromPk(pkTest[0])
		if err != nil {
			log.Fatal(err)
		}

		coeffs, _, ciphertext, _ := new_test_vectors(ckkscontext, encoder, encryptorTest, 1)

		verify_test_vectors(ckkscontext, encoder, decryptor_sk0, coeffs, ciphertext, t)

	})

	// CKS (OK multiple levels)
	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/CKS", parties, logN, logQ, levels, scale), func(t *testing.T) {

		coeffs, _, ciphertext, _ := new_test_vectors(ckkscontext, encoder, encryptor_pk0, 1)

		ciphertexts := make([]*ckks.Ciphertext, parties)
		for i := uint64(0); i < parties; i++ {
			ciphertexts[i] = ciphertext.CopyNew().Ciphertext()
		}

		// Each party creates its CKS instance with deltaSk = si-si'
		cks := make([]*CKS, parties)
		for i := uint64(0); i < parties; i++ {
			cks[i] = NewCKS(sk0_shards[i].Get(), sk1_shards[i].Get(), ckkscontext, 6.36)
		}

		// Each party computes its hi share from the shared ciphertext
		// Each party encodes its share and sends it to the other n-1 parties
		hi := make([]*ring.Poly, parties)
		for i := uint64(0); i < parties; i++ {
			hi[i] = cks[i].KeySwitch(ciphertexts[i].Value()[1])
		}
		// Each party receive the shares n-1 shares from the other parties and decodes them
		for i := uint64(0); i < parties; i++ {
			// Then keyswitch the ciphertext with the decoded shares
			cks[i].Aggregate(ciphertexts[i].Value()[0], hi)
		}

		for i := uint64(0); i < parties; i++ {

			verify_test_vectors(ckkscontext, encoder, decryptor_sk1, coeffs, ciphertexts[i], t)
		}
	})

	// PCKS (OK multiple levels)
	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/PCKS", parties, logN, logQ, levels, scale), func(t *testing.T) {

		coeffs, _, ciphertext, _ := new_test_vectors(ckkscontext, encoder, encryptor_pk0, 1)

		evaluator.DropLevel(ciphertext.Element(), 1)

		type Party struct {
			*PCKSProtocol
			s     *ring.Poly
			share PCKSShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.PCKSProtocol = NewPCKSProtocol(ckkscontext, 6.36)
			p.s = sk0_shards[i].Get()
			p.share = p.AllocateShares(ciphertext.Level())
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		ciphertextSwitched := ckkscontext.NewCiphertext(1, ciphertext.Level(), ciphertext.Scale())

		for i, p := range pcksParties {
			p.GenShare(p.s, pk1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

		verify_test_vectors(ckkscontext, encoder, decryptor_sk1, coeffs, ciphertextSwitched, t)
	})

	// CBOOT (OK multiple levels)
	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/CBOOT", parties, logN, logQ, levels, scale), func(t *testing.T) {

		coeffs, _, ciphertext, _ := new_test_vectors(ckkscontext, encoder, encryptor_pk0, 1)

		crpGenerator, _ := NewCRPGenerator(nil, contextCiphertexts)

		crp := crpGenerator.Clock()

		levelStart := uint64(3)

		refreshShares := make([]*RefreshShares, parties)
		for i := uint64(0); i < parties; i++ {
			refreshShares[i] = GenRefreshShares(sk0_shards[i], levelStart, parties, ckkscontext, ciphertext.Value()[1], crp)
		}

		for ciphertext.Level() != levelStart {
			evaluator.DropLevel(ciphertext.Element(), 1)
		}

		Refresh(ciphertext, refreshShares, ckkscontext, crp)

		verify_test_vectors(ckkscontext, encoder, decryptor_sk0, coeffs, ciphertext, t)

	})

}

func test_EKG_Protocol_Naive(parties uint64, sk []*ckks.SecretKey, collectivePk *ckks.PublicKey, ekgNaive []*EkgProtocolNaive) [][][2]*ring.Poly {

	// ROUND 0
	// Each party generates its samples
	samples := make([][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		samples[i] = ekgNaive[i].GenSamples(sk[i].Get(), collectivePk.Get())
	}

	// ROUND 1
	// Each party aggretates its sample with the other n-1 samples
	aggregatedSamples := make([][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		aggregatedSamples[i] = ekgNaive[i].Aggregate(sk[i].Get(), collectivePk.Get(), samples)
	}

	// ROUND 2
	// Each party aggregates sums its aggregatedSample with the other n-1 aggregated samples
	evk := make([][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		evk[i] = ekgNaive[i].Finalize(aggregatedSamples)
	}

	return evk
}

func test_EKG_Protocol(parties uint64, ekgProtocols []*EkgProtocol, sk []*ckks.SecretKey, ephemeralKeys []*ring.Poly, crp [][]*ring.Poly) [][][2]*ring.Poly {

	// ROUND 1
	samples := make([][]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		samples[i] = ekgProtocols[i].GenSamples(ephemeralKeys[i], sk[i].Get(), crp[i])
	}

	//ROUND 2
	aggregatedSamples := make([][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		aggregatedSamples[i] = ekgProtocols[i].Aggregate(sk[i].Get(), samples, crp[i])
	}

	// ROUND 3
	keySwitched := make([][]*ring.Poly, parties)
	sum := make([][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		sum[i] = ekgProtocols[i].Sum(aggregatedSamples)
		keySwitched[i] = ekgProtocols[i].KeySwitch(ephemeralKeys[i], sk[i].Get(), sum[i])
	}

	// ROUND 4
	collectiveEvaluationKey := make([][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		collectiveEvaluationKey[i] = ekgProtocols[i].ComputeEVK(keySwitched, sum[i])
	}

	return collectiveEvaluationKey
}

func new_test_vectors(ckkscontext *ckks.CkksContext, encoder *ckks.Encoder, encryptor *ckks.Encryptor, a float64) (values []complex128, plaintext *ckks.Plaintext, ciphertext *ckks.Ciphertext, err error) {

	slots := ckkscontext.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = randomComplex(a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = ckkscontext.NewPlaintext(ckkscontext.Levels()-1, ckkscontext.Scale())

	if err = encoder.Encode(plaintext, values, ckkscontext.Slots()); err != nil {
		return nil, nil, nil, err
	}

	ciphertext, err = encryptor.EncryptNew(plaintext)
	if err != nil {
		return nil, nil, nil, err
	}

	return values, plaintext, ciphertext, nil
}

func verify_test_vectors(ckkscontext *ckks.CkksContext, encoder *ckks.Encoder, decryptor *ckks.Decryptor, valuesWant []complex128, element interface{}, t *testing.T) (err error) {

	var plaintextTest *ckks.Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *ckks.Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*ckks.Ciphertext))
	case *ckks.Plaintext:
		plaintextTest = element.(*ckks.Plaintext)
	}

	valuesTest = encoder.Decode(plaintextTest, ckkscontext.Slots())

	var DeltaReal0, DeltaImag0, DeltaReal1, DeltaImag1 float64

	for i := range valuesWant {

		// Test for big values (> 1)
		DeltaReal0 = real(valuesWant[i]) / real(valuesTest[i])
		DeltaImag0 = imag(valuesWant[i]) / imag(valuesTest[i])

		// Test for small values (< 1)
		DeltaReal1 = real(valuesWant[i]) - real(valuesTest[i])
		DeltaImag1 = imag(valuesWant[i]) - imag(valuesTest[i])

		if DeltaReal1 < 0 {
			DeltaReal1 *= -1
		}
		if DeltaImag1 < 0 {
			DeltaImag1 *= -1
		}

		if (DeltaReal0 < 0.999 || DeltaReal0 > 1.001 || DeltaImag0 < 0.999 || DeltaImag0 > 1.001) && (DeltaReal1 > 0.001 || DeltaImag1 > 0.001) {
			t.Errorf("error : coeff %d, want %f have %f", i, valuesWant[i], valuesTest[i])
			break
		}
	}

	return nil
}
