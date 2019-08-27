package dckks

import (
	"fmt"
	"github.com/lca1/lattigo/ckks"
	"github.com/lca1/lattigo/ring"
	"log"
	"math"
	"math/rand"
	"testing"
	"time"
)

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func randomComplex(min, max float64) complex128 {
	return complex(randomFloat(min, max), randomFloat(min, max))
}

func Test_DBFVScheme(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	//sigmaSmudging := 6.36

	var err error

	var parties, logN, logQ, levels, logScale, bdc uint64
	parties = 5
	logN = 9
	logQ = 49
	levels = 12
	sigma := 3.19
	logScale = 40
	bdc = 20

	var ckkscontext *ckks.CkksContext

	log.Printf("Generating CkksContext for logN=%d/logQ=%d/levels=%d/logScale=%d/sigma=%f", logN, logQ, levels, logScale, sigma)
	if ckkscontext, err = ckks.NewCkksContext(logN, logQ, logScale, levels, sigma); err != nil {
		log.Fatal(err)
	}

	kgen := ckkscontext.NewKeyGenerator()

	evaluator := ckkscontext.NewEvaluator()

	context := ckkscontext.ContextKeys()

	coeffsWant := make([]complex128, ckkscontext.Slots())
	for i := uint64(0); i < ckkscontext.Slots(); i++ {
		coeffsWant[i] = randomComplex(-1, 1)
	}

	plaintextWant := ckkscontext.NewPlaintext(levels-1, logScale)
	if err = plaintextWant.EncodeComplex(ckkscontext, coeffsWant); err != nil {
		log.Fatal(err)
	}

	ciphertextTest := ckkscontext.NewCiphertext(1, levels-1, logScale)

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
	pk0, err := kgen.NewPublicKey(sk0)
	if err != nil {
		log.Fatal(err)
	}

	pk1, err := kgen.NewPublicKey(sk1)
	if err != nil {
		log.Fatal(err)
	}

	_ = pk1

	// Encryptors
	encryptor_pk0, err := ckkscontext.NewEncryptor(pk0)
	if err != nil {
		log.Fatal(err)
	}

	// Decryptors
	decryptor_sk0, err := ckkscontext.NewDecryptor(sk0)
	if err != nil {
		log.Fatal(err)
	}

	decryptor_sk1, err := ckkscontext.NewDecryptor(sk1)
	if err != nil {
		log.Fatal(err)
	}

	_ = decryptor_sk1

	// Reference ciphertext
	ciphertext, err := encryptor_pk0.EncryptNew(plaintextWant)
	if err != nil {
		log.Fatal(err)
	}

	coeffsMul := make([]complex128, ckkscontext.Slots())
	for i := uint64(0); i < ckkscontext.Slots(); i++ {
		coeffsMul[i] = coeffsWant[i] * coeffsWant[i]
	}

	evaluator.MulRelin(ciphertext, ciphertext, nil, ciphertext)

	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/logScale=%d/CRS_PRNG", parties, logN, logQ, levels, logScale), func(t *testing.T) {

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

	// EKG_Naive

	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/logScale=%d/bdc=%d/EKG", parties, logN, logQ, levels, logScale, bdc), func(t *testing.T) {

		bitLog := uint64(math.Ceil(float64(60) / float64(bdc)))

		// Each party instantiate an ekg naive protocole
		ekg := make([]*EkgProtocol, parties)
		ephemeralKeys := make([]*ring.Poly, parties)
		crp := make([][][]*ring.Poly, parties)

		for i := uint64(0); i < parties; i++ {

			ekg[i] = NewEkgProtocol(context, bdc)
			ephemeralKeys[i] = ekg[i].NewEphemeralKey()
			crp[i] = make([][]*ring.Poly, len(context.Modulus))

			for j := 0; j < len(context.Modulus); j++ {
				crp[i][j] = make([]*ring.Poly, bitLog)
				for u := uint64(0); u < bitLog; u++ {
					crp[i][j][u] = crpGenerators[i].Clock()
				}
			}
		}

		evk := test_EKG_Protocol(parties, ekg, sk0_shards, ephemeralKeys, crp)

		rlk, err := kgen.SetRelinKeys(evk[0], bdc)
		if err != nil {
			log.Fatal(err)
		}

		if err := evaluator.Relinearize(ciphertext, rlk, ciphertextTest); err != nil {
			log.Fatal(err)
		}

		verify_test_vectors(ckkscontext, decryptor_sk0, coeffsMul, ciphertextTest, t)

	})

	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/logScale=%d/bdc=%d/EKG_NAIVE", parties, logN, logQ, levels, logScale, bdc), func(t *testing.T) {

		// Each party instantiate an ekg naive protocole
		ekgNaive := make([]*EkgProtocolNaive, parties)
		for i := uint64(0); i < parties; i++ {
			ekgNaive[i] = NewEkgProtocolNaive(context, bdc)
		}

		evk := test_EKG_Protocol_Naive(parties, sk0_shards, pk0, ekgNaive)

		rlk, err := kgen.SetRelinKeys(evk[0], bdc)
		if err != nil {
			log.Fatal(err)
		}

		if err := evaluator.Relinearize(ciphertext, rlk, ciphertextTest); err != nil {
			log.Fatal(err)
		}

		verify_test_vectors(ckkscontext, decryptor_sk0, coeffsMul, ciphertextTest, t)
	})

	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/logScale=%d/CKG", parties, logN, logQ, levels, logScale), func(t *testing.T) {

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
		encryptorTest, err := ckkscontext.NewEncryptor(pkTest[0])
		if err != nil {
			log.Fatal(err)
		}

		ciphertextTest, err := encryptorTest.EncryptNew(plaintextWant)

		if err != nil {
			log.Fatal(err)
		}

		verify_test_vectors(ckkscontext, decryptor_sk0, coeffsWant, ciphertextTest, t)

	})

	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/logScale=%d/CKS", parties, logN, logQ, levels, logScale), func(t *testing.T) {

		ciphertext, err := encryptor_pk0.EncryptNew(plaintextWant)
		if err != nil {
			log.Fatal(err)
		}

		ciphertexts := make([]*ckks.Ciphertext, parties)
		for i := uint64(0); i < parties; i++ {
			ciphertexts[i] = ciphertext.CopyNew().(*ckks.Ciphertext)
		}

		// Each party creates its CKS instance with deltaSk = si-si'
		cks := make([]*CKS, parties)
		for i := uint64(0); i < parties; i++ {
			cks[i] = NewCKS(sk0_shards[i].Get(), sk1_shards[i].Get(), context, 6.36)
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

			verify_test_vectors(ckkscontext, decryptor_sk1, coeffsWant, ciphertexts[i], t)
		}
	})

	t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/logScale=%d/PCKS", parties, logN, logQ, levels, logScale), func(t *testing.T) {

		ciphertext, err := encryptor_pk0.EncryptNew(plaintextWant)
		if err != nil {
			log.Fatal(err)
		}

		ciphertexts := make([]*ckks.Ciphertext, parties)
		for i := uint64(0); i < parties; i++ {
			ciphertexts[i] = ciphertext.CopyNew().(*ckks.Ciphertext)
		}

		pcks := make([]*PCKS, parties)
		for i := uint64(0); i < parties; i++ {
			pcks[i] = NewPCKS(sk0_shards[i].Get(), pk1.Get(), context, 6.36)
		}

		hi := make([][2]*ring.Poly, parties)
		for i := uint64(0); i < parties; i++ {
			hi[i] = pcks[i].KeySwitch(ciphertexts[i].Value()[1])
		}

		for i := uint64(0); i < parties; i++ {
			pcks[i].Aggregate(ciphertexts[i].Value(), hi)
		}

		for i := uint64(0); i < parties; i++ {

			verify_test_vectors(ckkscontext, decryptor_sk1, coeffsWant, ciphertexts[i], t)

		}
	})
}

func test_EKG_Protocol_Naive(parties uint64, sk []*ckks.SecretKey, collectivePk *ckks.PublicKey, ekgNaive []*EkgProtocolNaive) [][][][2]*ring.Poly {

	// ROUND 0
	// Each party generates its samples
	samples := make([][][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		samples[i] = ekgNaive[i].GenSamples(sk[i].Get(), collectivePk.Get())
	}

	// ROUND 1
	// Each party aggretates its sample with the other n-1 samples
	aggregatedSamples := make([][][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		aggregatedSamples[i] = ekgNaive[i].Aggregate(sk[i].Get(), collectivePk.Get(), samples)
	}

	// ROUND 2
	// Each party aggregates sums its aggregatedSample with the other n-1 aggregated samples
	evk := make([][][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		evk[i] = ekgNaive[i].Finalize(aggregatedSamples)
	}

	return evk
}

func test_EKG_Protocol(parties uint64, ekgProtocols []*EkgProtocol, sk []*ckks.SecretKey, ephemeralKeys []*ring.Poly, crp [][][]*ring.Poly) [][][][2]*ring.Poly {

	// ROUND 1
	samples := make([][][]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		samples[i] = ekgProtocols[i].GenSamples(ephemeralKeys[i], sk[i].Get(), crp[i])
	}

	//ROUND 2
	aggregatedSamples := make([][][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		aggregatedSamples[i] = ekgProtocols[i].Aggregate(sk[i].Get(), samples, crp[i])
	}

	// ROUND 3
	keySwitched := make([][][]*ring.Poly, parties)
	sum := make([][][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		sum[i] = ekgProtocols[i].Sum(aggregatedSamples)
		keySwitched[i] = ekgProtocols[i].KeySwitch(ephemeralKeys[i], sk[i].Get(), sum[i])
	}

	// ROUND 4
	collectiveEvaluationKey := make([][][][2]*ring.Poly, parties)
	for i := uint64(0); i < parties; i++ {
		collectiveEvaluationKey[i] = ekgProtocols[i].ComputeEVK(keySwitched, sum[i])
	}

	return collectiveEvaluationKey
}

func verify_test_vectors(ckkscontext *ckks.CkksContext, decryptor *ckks.Decryptor, valuesWant []complex128, element ckks.CkksElement, t *testing.T) (err error) {

	var plaintextTest *ckks.Plaintext
	var valuesTest []complex128

	if element.Degree() == 0 {

		plaintextTest = element.(*ckks.Plaintext)

	} else {

		if plaintextTest, err = decryptor.DecryptNew(element.(*ckks.Ciphertext)); err != nil {
			return err
		}
	}

	valuesTest = plaintextTest.DecodeComplex(ckkscontext)

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
