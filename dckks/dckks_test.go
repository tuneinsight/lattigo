package dckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"log"
	"math"
	"sort"
	"testing"
)

func check(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
	}
}

type dckksContext struct {
	ckksContext *ckks.CkksContext
	encoder     *ckks.Encoder
	evaluator   *ckks.Evaluator

	encryptorPk0 *ckks.Encryptor
	decryptorSk0 *ckks.Decryptor
	decryptorSk1 *ckks.Decryptor

	pk0 *ckks.PublicKey
	pk1 *ckks.PublicKey

	sk0 *ckks.SecretKey
	sk1 *ckks.SecretKey

	sk0Shards []*ckks.SecretKey
	sk1Shards []*ckks.SecretKey
}

type dckksTestParameters struct {
	parties    uint64
	verbose    bool
	medianprec float64

	contexts []*dckksContext
}

var err error
var testParams = new(dckksTestParameters)

func init() {

	testParams.parties = 5

	testParams.medianprec = 15
	testParams.verbose = false

	parameters := []*ckks.Parameters{
		ckks.DefaultParams[13],
		ckks.DefaultParams[14],
		//ckks.DefaultParams[15],
		//ckks.DefaultParams[16],
	}

	testParams.contexts = make([]*dckksContext, len(parameters))

	for i, parameter := range parameters {

		testParams.contexts[i] = new(dckksContext)

		if testParams.contexts[i].ckksContext, err = ckks.NewCkksContext(parameter); err != nil {
			log.Fatal(err)
		}

		testParams.contexts[i].encoder = testParams.contexts[i].ckksContext.NewEncoder()
		testParams.contexts[i].evaluator = testParams.contexts[i].ckksContext.NewEvaluator()

		kgen := testParams.contexts[i].ckksContext.NewKeyGenerator()

		// SecretKeys
		testParams.contexts[i].sk0Shards = make([]*ckks.SecretKey, testParams.parties)
		testParams.contexts[i].sk1Shards = make([]*ckks.SecretKey, testParams.parties)
		tmp0 := testParams.contexts[i].ckksContext.ContextKeys().NewPoly()
		tmp1 := testParams.contexts[i].ckksContext.ContextKeys().NewPoly()

		for j := uint64(0); j < testParams.parties; j++ {
			testParams.contexts[i].sk0Shards[j] = kgen.NewSecretKey()
			testParams.contexts[i].sk1Shards[j] = kgen.NewSecretKey()
			testParams.contexts[i].ckksContext.ContextKeys().Add(tmp0, testParams.contexts[i].sk0Shards[j].Get(), tmp0)
			testParams.contexts[i].ckksContext.ContextKeys().Add(tmp1, testParams.contexts[i].sk1Shards[j].Get(), tmp1)
		}

		testParams.contexts[i].sk0 = new(ckks.SecretKey)
		testParams.contexts[i].sk1 = new(ckks.SecretKey)

		testParams.contexts[i].sk0.Set(tmp0)
		testParams.contexts[i].sk1.Set(tmp1)

		// Publickeys
		testParams.contexts[i].pk0 = kgen.NewPublicKey(testParams.contexts[i].sk0)
		testParams.contexts[i].pk1 = kgen.NewPublicKey(testParams.contexts[i].sk1)

		if testParams.contexts[i].encryptorPk0, err = testParams.contexts[i].ckksContext.NewEncryptorFromPk(testParams.contexts[i].pk0); err != nil {
			log.Fatal(err)
		}

		if testParams.contexts[i].decryptorSk0, err = testParams.contexts[i].ckksContext.NewDecryptor(testParams.contexts[i].sk0); err != nil {
			log.Fatal(err)
		}

		if testParams.contexts[i].decryptorSk1, err = testParams.contexts[i].ckksContext.NewDecryptor(testParams.contexts[i].sk1); err != nil {
			log.Fatal(err)
		}
	}
}

func Test_DCKKS_CRP(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

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

				crs_generator_1, _ := NewCRPGenerator(nil, ckksContext.ContextKeys())
				crs_generator_2, _ := NewCRPGenerator(nil, ckksContext.ContextKeys())

				crs_generator_1.Seed(seed1)
				crs_generator_2.Seed(append(seed1, seed2...)) //Append works since blake2b hashes blocks of 512 bytes

				crs_generator_1.SetClock(256)
				crs_generator_2.SetClock(255)

				p0 := crs_generator_1.Clock()
				p1 := crs_generator_2.Clock()

				if ckksContext.ContextKeys().Equal(p0, p1) != true {
					t.Errorf("error : crs prng generator")
				}
			})
	}
}

func Test_DCKKS_CKG(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

				crpGenerator, err := NewCRPGenerator(nil, ckksContext.ContextKeys())
				check(t, err)
				crpGenerator.Seed([]byte{})
				crp := crpGenerator.Clock()

				type Party struct {
					*CKGProtocol
					s  *ring.Poly
					s1 CKGShare
				}

				ckgParties := make([]*Party, parties)
				for i := uint64(0); i < parties; i++ {
					p := new(Party)
					p.CKGProtocol = NewCKGProtocol(ckksContext)
					p.s = sk0Shards[i].Get()
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

				pk := &ckks.PublicKey{}
				P0.GenPublicKey(P0.s1, crp, pk)

				// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
				encryptorTest, err := ckksContext.NewEncryptorFromPk(pk)
				if err != nil {
					t.Error(err)
				}

				coeffs, _, ciphertext := new_test_vectors(params, encryptorTest, 1, t)

				verify_test_vectors(params, decryptorSk0, coeffs, ciphertext, t)
			})
	}
}

func Test_DCKKS_EKG(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext
		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

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
					p.RKGProtocol = NewEkgProtocol(ckksContext)
					p.u, err = p.RKGProtocol.NewEphemeralKey(1.0 / 3.0)
					check(t, err)
					p.s = sk0Shards[i].Get()
					p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()
					rkgParties[i] = p
				}

				P0 := rkgParties[0]

				crpGenerator, err := NewCRPGenerator(nil, ckksContext.ContextKeys())
				check(t, err)
				crpGenerator.Seed([]byte{})
				crp := make([]*ring.Poly, ckksContext.Beta())

				for i := uint64(0); i < ckksContext.Beta(); i++ {
					crp[i] = crpGenerator.Clock()
				}

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

				evk := ckksContext.NewRelinKeyEmpty()
				P0.GenRelinearizationKey(P0.share2, P0.share3, evk)

				coeffs, _, ciphertext := new_test_vectors(params, encryptorPk0, 1, t)

				for i := range coeffs {
					coeffs[i] *= coeffs[i]
				}

				if err := evaluator.MulRelin(ciphertext, ciphertext, evk, ciphertext); err != nil {
					log.Fatal(err)
				}

				evaluator.Rescale(ciphertext, ckksContext.Scale(), ciphertext)

				if ciphertext.Degree() != 1 {
					t.Errorf("EKG_NAIVE -> bad relinearize")
				}

				verify_test_vectors(params, decryptorSk0, coeffs, ciphertext, t)

			})
	}
}

func Test_DCKKS_EKG_Naive(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext
		evaluator := params.evaluator
		pk0 := params.pk0
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

				// Each party instantiate an ekg naive protocole
				ekgNaive := make([]*EkgProtocolNaive, parties)
				for i := uint64(0); i < parties; i++ {
					ekgNaive[i] = NewEkgProtocolNaive(ckksContext)
				}

				// ROUND 0
				// Each party generates its samples
				samples := make([][][2]*ring.Poly, parties)
				for i := uint64(0); i < parties; i++ {
					samples[i] = ekgNaive[i].GenSamples(sk0Shards[i].Get(), pk0.Get())
				}

				// ROUND 1
				// Each party aggretates its sample with the other n-1 samples
				aggregatedSamples := make([][][2]*ring.Poly, parties)
				for i := uint64(0); i < parties; i++ {
					aggregatedSamples[i] = ekgNaive[i].Aggregate(sk0Shards[i].Get(), pk0.Get(), samples)
				}

				// ROUND 2
				// Each party aggregates sums its aggregatedSample with the other n-1 aggregated samples
				rlk := new(ckks.EvaluationKey)
				rlk.Set(ekgNaive[0].Finalize(aggregatedSamples))

				coeffs, _, ciphertext := new_test_vectors(params, encryptorPk0, 1, t)

				for i := range coeffs {
					coeffs[i] *= coeffs[i]
				}

				if err := evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext); err != nil {
					log.Fatal(err)
				}

				evaluator.Rescale(ciphertext, ckksContext.Scale(), ciphertext)

				if ciphertext.Degree() != 1 {
					t.Errorf("EKG_NAIVE -> bad relinearize")
				}

				verify_test_vectors(params, decryptorSk0, coeffs, ciphertext, t)
			})
	}
}

func Test_DCKKS_RKG_Conjugate(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext
		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

				rkg := make([]*RKG, parties)

				for i := uint64(0); i < parties; i++ {
					rkg[i] = NewRKG(ckksContext)
				}

				crpGenerator, err := NewCRPGenerator(nil, ckksContext.ContextKeys())
				check(t, err)
				crpGenerator.Seed([]byte{})
				crp := make([]*ring.Poly, ckksContext.Beta())

				for i := uint64(0); i < ckksContext.Beta(); i++ {
					crp[i] = crpGenerator.Clock()
				}

				shares := make([][]*ring.Poly, parties)
				for i := uint64(0); i < parties; i++ {
					shares[i] = rkg[i].GenShareRotRow(sk0Shards[i].Get(), crp)
				}

				rkg[0].AggregateRotRow(shares, crp)
				rotkey := rkg[0].Finalize(ckksContext)

				coeffs, _, ciphertext := new_test_vectors(params, encryptorPk0, 1, t)

				err = evaluator.Conjugate(ciphertext, rotkey, ciphertext)
				check(t, err)

				for i := range coeffs {
					coeffs[i] = complex(real(coeffs[i]), -imag(coeffs[i]))
				}

				verify_test_vectors(params, decryptorSk0, coeffs, ciphertext, t)

			})
	}
}

func Test_DCKKS_RKG_RotCols(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext
		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

				rkg := make([]*RKG, parties)

				for i := uint64(0); i < parties; i++ {
					rkg[i] = NewRKG(ckksContext)
				}

				crpGenerator, err := NewCRPGenerator(nil, ckksContext.ContextKeys())
				check(t, err)
				crpGenerator.Seed([]byte{})
				crp := make([]*ring.Poly, ckksContext.Beta())

				for i := uint64(0); i < ckksContext.Beta(); i++ {
					crp[i] = crpGenerator.Clock()
				}

				coeffs, _, ciphertext := new_test_vectors(params, encryptorPk0, 1, t)
				mask := ckksContext.Slots() - 1

				receiver := ckksContext.NewCiphertext(ciphertext.Degree(), ciphertext.Level(), ciphertext.Scale())
				for n := uint64(0); n < ckksContext.LogN(); n++ {

					shares := make([][]*ring.Poly, parties)
					for i := uint64(0); i < parties; i++ {
						shares[i] = rkg[i].GenShareRotLeft(sk0Shards[i].Get(), 1<<n, crp)
					}

					rkg[0].AggregateRotColL(shares, 1<<n, crp)
					rotkey := rkg[0].Finalize(ckksContext)

					err = evaluator.RotateColumns(ciphertext, 1<<n, rotkey, receiver)
					check(t, err)

					coeffsWant := make([]complex128, ckksContext.Slots())

					for i := uint64(0); i < ckksContext.Slots(); i++ {
						coeffsWant[i] = coeffs[(i+(1<<n))&mask]
					}

					verify_test_vectors(params, decryptorSk0, coeffsWant, receiver, t)
				}

			})
	}
}

func Test_DCKKS_CKS(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext
		encryptorPk0 := params.encryptorPk0
		decryptorSk1 := params.decryptorSk1
		sk0Shards := params.sk0Shards
		sk1Shards := params.sk1Shards

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

				type Party struct {
					*CKSProtocol
					s0    *ring.Poly
					s1    *ring.Poly
					share CKSShare
				}

				cksParties := make([]*Party, parties)
				for i := uint64(0); i < parties; i++ {
					p := new(Party)
					p.CKSProtocol = NewCKSProtocol(ckksContext, 6.36)
					p.s0 = sk0Shards[i].Get()
					p.s1 = sk1Shards[i].Get()
					p.share = p.AllocateShare()
					cksParties[i] = p
				}
				P0 := cksParties[0]

				coeffs, _, ciphertext := new_test_vectors(params, encryptorPk0, 1, t)

				// Each party creates its CKSProtocol instance with tmp = si-si'
				for i, p := range cksParties {
					p.GenShare(p.s0, p.s1, ciphertext, p.share)
					if i > 0 {
						P0.AggregateShares(p.share, P0.share, P0.share)
					}
				}

				ksCiphertext := ckksContext.NewCiphertext(1, ciphertext.Level(), ciphertext.Scale())

				P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

				verify_test_vectors(params, decryptorSk1, coeffs, ksCiphertext, t)

				P0.KeySwitch(P0.share, ciphertext, ciphertext)

				verify_test_vectors(params, decryptorSk1, coeffs, ksCiphertext, t)

			})
	}
}

func Test_DCKKS_PCKS(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext
		encryptorPk0 := params.encryptorPk0
		decryptorSk1 := params.decryptorSk1
		sk0Shards := params.sk0Shards
		pk1 := params.pk1

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

				coeffs, _, ciphertext := new_test_vectors(params, encryptorPk0, 1, t)

				type Party struct {
					*PCKSProtocol
					s     *ring.Poly
					share PCKSShare
				}

				pcksParties := make([]*Party, parties)
				for i := uint64(0); i < parties; i++ {
					p := new(Party)
					p.PCKSProtocol = NewPCKSProtocol(ckksContext, 6.36)
					p.s = sk0Shards[i].Get()
					p.share = p.AllocateShares(ciphertext.Level())
					pcksParties[i] = p
				}
				P0 := pcksParties[0]

				ciphertextSwitched := ckksContext.NewCiphertext(1, ciphertext.Level(), ciphertext.Scale())

				for i, p := range pcksParties {
					p.GenShare(p.s, pk1, ciphertext, p.share)
					if i > 0 {
						P0.AggregateShares(p.share, P0.share, P0.share)
					}
				}

				P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

				verify_test_vectors(params, decryptorSk1, coeffs, ciphertextSwitched, t)
			})
	}
}

func Test_DCKKS_Refresh(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		ckksContext := params.ckksContext
		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f",
			parties,
			ckksContext.LogN(),
			ckksContext.LogQ(),
			ckksContext.Levels(),
			ckksContext.Scale()),
			func(t *testing.T) {

				coeffs, _, ciphertext := new_test_vectors(params, encryptorPk0, 1.0, t)

				crpGenerator, _ := NewCRPGenerator(nil, ckksContext.ContextQ())

				crp := crpGenerator.Clock()

				levelStart := uint64(3)

				refreshShares := make([]*RefreshShares, parties)
				for i := uint64(0); i < parties; i++ {
					refreshShares[i] = GenRefreshShares(sk0Shards[i], levelStart, parties, ckksContext, ciphertext.Value()[1], crp)
				}

				for ciphertext.Level() != levelStart {
					evaluator.DropLevel(ciphertext.Element(), 1)
				}

				Refresh(ciphertext, refreshShares, ckksContext, crp)

				verify_test_vectors(params, decryptorSk0, coeffs, ciphertext, t)

			})
	}
}

func new_test_vectors(contextParams *dckksContext, encryptor *ckks.Encryptor, a float64, t *testing.T) (values []complex128, plaintext *ckks.Plaintext, ciphertext *ckks.Ciphertext) {

	slots := contextParams.ckksContext.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = randomComplex(a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = contextParams.ckksContext.NewPlaintext(contextParams.ckksContext.Levels()-1, contextParams.ckksContext.Scale())

	err = contextParams.encoder.Encode(plaintext, values, slots)
	check(t, err)

	ciphertext, err = encryptor.EncryptNew(plaintext)
	check(t, err)

	return values, plaintext, ciphertext
}

func verify_test_vectors(contextParams *dckksContext, decryptor *ckks.Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	var plaintextTest *ckks.Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *ckks.Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*ckks.Ciphertext))
	case *ckks.Plaintext:
		plaintextTest = element.(*ckks.Plaintext)
	}

	valuesTest = contextParams.encoder.Decode(plaintextTest, contextParams.ckksContext.Slots())

	var deltaReal, deltaImag float64

	var minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, contextParams.ckksContext.Slots())

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	distrib_real := make(map[uint64]uint64)
	distrib_imag := make(map[uint64]uint64)

	for i := range valuesWant {

		deltaReal = math.Abs(real(valuesTest[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesTest[i]) - imag(valuesWant[i]))

		diff[i] += complex(deltaReal, 0)
		diff[i] += complex(0, deltaImag)

		meanprec += diff[i]

		if real(diff[i]) > real(minprec) {
			minprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) > imag(minprec) {
			minprec = complex(real(minprec), imag(diff[i]))
		}

		if real(diff[i]) < real(maxprec) {
			maxprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) < imag(maxprec) {
			maxprec = complex(real(maxprec), imag(diff[i]))
		}

		distrib_real[uint64(math.Floor(math.Log2(1/real(diff[i]))))] += 1
		distrib_imag[uint64(math.Floor(math.Log2(1/imag(diff[i]))))] += 1
	}

	meanprec /= complex(float64(contextParams.ckksContext.Slots()), 0)
	medianprec = calcmedian(diff)

	if testParams.verbose {
		t.Logf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
		t.Logf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
		t.Logf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
		t.Logf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
		t.Log()
	}

	if math.Log2(1/real(medianprec)) < testParams.medianprec || math.Log2(1/imag(medianprec)) < testParams.medianprec {
		t.Errorf("Mean precision error : target (%.2f, %.2f) > result (%.2f, %.2f)", testParams.medianprec, testParams.medianprec, math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
	}
}

func calcmedian(values []complex128) (median complex128) {

	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = real(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(tmp[i], imag(values[i]))
	}

	for i := range values {
		tmp[i] = imag(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(real(values[i]), tmp[i])
	}

	index := len(values) / 2

	if len(values)&1 == 1 {
		return values[index]
	}

	if index+1 == len(values) {
		return values[index]
	}

	return (values[index] + values[index+1]) / 2
}
