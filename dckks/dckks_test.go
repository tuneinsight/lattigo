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

	ckksParameters []*ckks.Parameters
}

var err error
var testParams = new(dckksTestParameters)

func init() {

	testParams.parties = 5

	testParams.medianprec = 15
	testParams.verbose = false

	testParams.ckksParameters = []*ckks.Parameters{
		ckks.DefaultParams[13],
		ckks.DefaultParams[14],
		//ckks.DefaultParams[15],
		//ckks.DefaultParams[16],
	}
}

func Test_DCKKS(t *testing.T) {
	t.Run("PublicKeyGen", testPublicKeyGen)
	t.Run("RelinKeyGen", testRelinKeyGen)
	t.Run("RelinKeyGenNaive", testRelinKeyGenNaive)
	t.Run("KeySwitching", testKeyswitching)
	t.Run("PublicKeySwitching", testPublicKeySwitching)
	t.Run("RotKeyGenConjugate", testRotKeyGenConjugate)
	t.Run("RotKeyGenCols", testRotKeyGenCols)
	t.Run("Refresh", testRefresh)
}

func genDCKKSContext(contextParameters *ckks.Parameters) (params *dckksContext) {

	params = new(dckksContext)

	if params.ckksContext, err = ckks.NewCkksContext(contextParameters); err != nil {
		log.Fatal(err)
	}

	params.encoder = params.ckksContext.NewEncoder()
	params.evaluator = params.ckksContext.NewEvaluator()

	kgen := params.ckksContext.NewKeyGenerator()

	// SecretKeys
	params.sk0Shards = make([]*ckks.SecretKey, testParams.parties)
	params.sk1Shards = make([]*ckks.SecretKey, testParams.parties)
	tmp0 := params.ckksContext.ContextKeys().NewPoly()
	tmp1 := params.ckksContext.ContextKeys().NewPoly()

	for j := uint64(0); j < testParams.parties; j++ {
		params.sk0Shards[j] = kgen.NewSecretKey()
		params.sk1Shards[j] = kgen.NewSecretKey()
		params.ckksContext.ContextKeys().Add(tmp0, params.sk0Shards[j].Get(), tmp0)
		params.ckksContext.ContextKeys().Add(tmp1, params.sk1Shards[j].Get(), tmp1)
	}

	params.sk0 = new(ckks.SecretKey)
	params.sk1 = new(ckks.SecretKey)

	params.sk0.Set(tmp0)
	params.sk1.Set(tmp1)

	// Publickeys
	params.pk0 = kgen.NewPublicKey(params.sk0)
	params.pk1 = kgen.NewPublicKey(params.sk1)

	if params.encryptorPk0, err = params.ckksContext.NewEncryptorFromPk(params.pk0); err != nil {
		log.Fatal(err)
	}

	if params.decryptorSk0, err = params.ckksContext.NewDecryptor(params.sk0); err != nil {
		log.Fatal(err)
	}

	if params.decryptorSk1, err = params.ckksContext.NewDecryptor(params.sk1); err != nil {
		log.Fatal(err)
	}

	return
}

func testPublicKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

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

				crpGenerator, err := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
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

func testRelinKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

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
					p.u, err = p.NewEphemeralKey(1.0 / 3.0)
					check(t, err)
					p.s = sk0Shards[i].Get()
					p.share1, p.share2, p.share3 = p.AllocateShares()
					rkgParties[i] = p
				}

				P0 := rkgParties[0]

				crpGenerator, err := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
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

func testRelinKeyGenNaive(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

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

				type Party struct {
					*RKGProtocolNaive
					u      *ring.Poly
					s      *ring.Poly
					share1 RKGNaiveShareRoundOne
					share2 RKGNaiveShareRoundTwo
				}

				rkgParties := make([]*Party, parties)

				for i := range rkgParties {
					p := new(Party)
					p.RKGProtocolNaive = NewRKGProtocolNaive(ckksContext)
					p.s = sk0Shards[i].Get()
					p.share1, p.share2 = p.AllocateShares()
					rkgParties[i] = p
				}

				P0 := rkgParties[0]

				// ROUND 1
				for i, p := range rkgParties {
					rkgParties[i].GenShareRoundOne(p.s, pk0.Get(), p.share1)
					if i > 0 {
						P0.AggregateShareRoundOne(p.share1, P0.share1, P0.share1)
					}
				}

				// ROUND 2
				for i, p := range rkgParties {
					rkgParties[i].GenShareRoundTwo(P0.share1, p.s, pk0.Get(), p.share2)
					if i > 0 {
						P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
					}
				}

				evk := ckksContext.NewRelinKeyEmpty()
				P0.GenRelinearizationKey(P0.share2, evk)

				coeffs, _, ciphertext := new_test_vectors(params, encryptorPk0, 1, t)

				for i := range coeffs {
					coeffs[i] *= coeffs[i]
				}

				if err := evaluator.MulRelin(ciphertext, ciphertext, evk, ciphertext); err != nil {
					log.Fatal(err)
				}

				if ciphertext.Degree() != 1 {
					t.Errorf("EKG_NAIVE -> bad relinearize")
				}

				evaluator.Rescale(ciphertext, ckksContext.Scale(), ciphertext)

				verify_test_vectors(params, decryptorSk0, coeffs, ciphertext, t)
			})
	}
}

func testKeyswitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

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

func testPublicKeySwitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

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

func testRotKeyGenConjugate(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

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

				crpGenerator, err := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
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

func testRotKeyGenCols(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

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

				crpGenerator, err := ring.NewCRPGenerator(nil, ckksContext.ContextKeys())
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

func testRefresh(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := genDCKKSContext(parameters)

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

				crpGenerator, _ := ring.NewCRPGenerator(nil, ckksContext.ContextQ())

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
