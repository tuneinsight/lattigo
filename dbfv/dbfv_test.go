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

type dbfvContext struct {
	bfvContext *bfv.BfvContext
	encoder    *bfv.BatchEncoder
	kgen       *bfv.KeyGenerator

	sk0Shards []*bfv.SecretKey
	sk0       *bfv.SecretKey

	sk1       *bfv.SecretKey
	sk1Shards []*bfv.SecretKey

	pk0 *bfv.PublicKey
	pk1 *bfv.PublicKey

	encryptorPk0 *bfv.Encryptor
	decryptorSk0 *bfv.Decryptor
	decryptorSk1 *bfv.Decryptor
	evaluator    *bfv.Evaluator
}

type dbfvTestParameters struct {
	parties uint64

	contexts []*dbfvContext
}

var err error
var testParams = new(dbfvTestParameters)

func init() {

	var err error

	testParams.parties = 5

	parameters := bfv.DefaultParams[0:2]

	//sigmaSmudging := 6.36

	testParams.contexts = make([]*dbfvContext, len(parameters))

	for i, params := range parameters {

		testParams.contexts[i] = new(dbfvContext)

		// nParties data indpendant element
		testParams.contexts[i].bfvContext = bfv.NewBfvContext()
		if err := testParams.contexts[i].bfvContext.SetParameters(&params); err != nil {
			log.Fatal(err)
		}

		testParams.contexts[i].kgen = testParams.contexts[i].bfvContext.NewKeyGenerator()

		testParams.contexts[i].evaluator = testParams.contexts[i].bfvContext.NewEvaluator()

		if testParams.contexts[i].encoder, err = testParams.contexts[i].bfvContext.NewBatchEncoder(); err != nil {
			log.Fatal(err)
		}

		// SecretKeys
		testParams.contexts[i].sk0Shards = make([]*bfv.SecretKey, testParams.parties)
		testParams.contexts[i].sk1Shards = make([]*bfv.SecretKey, testParams.parties)
		tmp0 := testParams.contexts[i].bfvContext.ContextKeys().NewPoly()
		tmp1 := testParams.contexts[i].bfvContext.ContextKeys().NewPoly()

		for j := uint64(0); j < testParams.parties; j++ {
			testParams.contexts[i].sk0Shards[j] = testParams.contexts[i].kgen.NewSecretKey()
			testParams.contexts[i].sk1Shards[j] = testParams.contexts[i].kgen.NewSecretKey()
			testParams.contexts[i].bfvContext.ContextKeys().Add(tmp0, testParams.contexts[i].sk0Shards[j].Get(), tmp0)
			testParams.contexts[i].bfvContext.ContextKeys().Add(tmp1, testParams.contexts[i].sk1Shards[j].Get(), tmp1)
		}

		testParams.contexts[i].sk0 = new(bfv.SecretKey)
		testParams.contexts[i].sk1 = new(bfv.SecretKey)

		testParams.contexts[i].sk0.Set(tmp0)
		testParams.contexts[i].sk1.Set(tmp1)

		// Publickeys
		testParams.contexts[i].pk0 = testParams.contexts[i].kgen.NewPublicKey(testParams.contexts[i].sk0)
		testParams.contexts[i].pk1 = testParams.contexts[i].kgen.NewPublicKey(testParams.contexts[i].sk1)

		if testParams.contexts[i].encryptorPk0, err = testParams.contexts[i].bfvContext.NewEncryptorFromPk(testParams.contexts[i].pk0); err != nil {
			log.Fatal(err)
		}

		if testParams.contexts[i].decryptorSk0, err = testParams.contexts[i].bfvContext.NewDecryptor(testParams.contexts[i].sk0); err != nil {
			log.Fatal(err)
		}

		if testParams.contexts[i].decryptorSk1, err = testParams.contexts[i].bfvContext.NewDecryptor(testParams.contexts[i].sk1); err != nil {
			log.Fatal(err)
		}

	}
}

func Test_DBFV_CKG(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		sk0Shards := params.sk0Shards
		decryptorSk0 := params.decryptorSk0

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

			crpGenerator, err := ring.NewCRPGenerator(nil, contextKeys)
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
				p.CKGProtocol = NewCKGProtocol(bfvContext)
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

			pk := &bfv.PublicKey{}
			P0.GenPublicKey(P0.s1, crp, pk)

			// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
			encryptorTest, err := bfvContext.NewEncryptorFromPk(pk)
			check(t, err)

			coeffs, _, ciphertext := newTestVectors(params, encryptorTest, t)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertext, t)
		})
	}
}

func Test_DBFV_EKG(t *testing.T) {

	parties := testParams.parties

	for _, params := range testParams.contexts {

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		sk0Shards := params.sk0Shards
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		evaluator := params.evaluator

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

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
				p.RKGProtocol = NewEkgProtocol(bfvContext)
				p.u, err = p.RKGProtocol.NewEphemeralKey(1.0 / 3.0)
				check(t, err)
				p.s = sk0Shards[i].Get()
				p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()
				rkgParties[i] = p
			}

			P0 := rkgParties[0]

			crpGenerator, err := ring.NewCRPGenerator(nil, bfvContext.ContextKeys())
			check(t, err)
			crpGenerator.Seed([]byte{})
			crp := make([]*ring.Poly, bfvContext.Beta())

			for i := uint64(0); i < bfvContext.Beta(); i++ {
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

			evk := bfvContext.NewRelinKeyEmpty(1)
			P0.GenRelinearizationKey(P0.share2, P0.share3, evk)

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, t)

			for i := range coeffs {
				coeffs[i] *= coeffs[i]
				coeffs[i] %= bfvContext.ContextT().Modulus[0]
			}

			ciphertextMul := bfvContext.NewCiphertext(ciphertext.Degree() * 2)
			err = evaluator.Mul(ciphertext, ciphertext, ciphertextMul)
			check(t, err)

			res := bfvContext.NewCiphertext(1)
			err = evaluator.Relinearize(ciphertextMul, evk, res)
			check(t, err)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertextMul, t)
		})
	}
}

/*
func test_EKGNaive(params *testParams.contexts[i], t *testing.T) {

	var err error

	contextKeys := params.contextKeys
	parties := params.parties
	bfvContext := params.bfvContext
	sk0Shards := params.sk0Shards
	pk0 := params.pk0
	evaluator := params.evaluator

	t.Run(fmt.Sprintf("N=%d/logQ=%d/EKG_Naive", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

		// Each party instantiate an ekg naive protocole
		ekgNaive := make([]*EkgProtocolNaive, parties)
		for i := 0; i < parties; i++ {
			ekgNaive[i] = NewEkgProtocolNaive(bfvContext)
		}

		evk := test_EKG_Protocol_Naive(parties, sk0Shards, pk0, ekgNaive)

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

		verifyTestVectors(params.decryptorSk0, params.encoder, coeffs, ciphertextMul, t)
	})
}

func test_RKG(params *testParams.contexts[i], t *testing.T) {

	var err error

	contextKeys := params.contextKeys
	parties := params.parties
	bfvContext := params.bfvContext
	keygenerator := params.kgen
	sk0Shards := params.sk0Shards
	evaluator := params.evaluator

	rkg := make([]*RKG, parties)
	for i := 0; i < parties; i++ {
		rkg[i] = NewRKG(bfvContext)
	}

	crpGenerator, err := ring.NewCRPGenerator(nil, contextKeys)
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
			shares[i] = rkg[i].GenShareRotRow(sk0Shards[i].Get(), crp)
		}

		rkg[0].AggregateRotRow(shares, crp)
		rotkey := rkg[0].Finalize(keygenerator)

		receiver := bfvContext.NewCiphertext(1)
		if err = evaluator.RotateRows(ciphertext, rotkey, receiver); err != nil {
			log.Fatal(err)
		}

		coeffs = append(coeffs[contextKeys.N>>1:], coeffs[:contextKeys.N>>1]...)

		verifyTestVectors(params.decryptorSk0, params.encoder, coeffs, receiver, t)

	})

	t.Run(fmt.Sprintf("N=%d/logQ=%d/RKG_rot_col_pow2", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

		coeffs, _, ciphertext, _ := newTestVectors(params)
		mask := (contextKeys.N >> 1) - 1

		receiver := bfvContext.NewCiphertext(ciphertext.Degree())
		for n := uint64(1); n < contextKeys.N>>1; n <<= 1 {

			shares := make([][]*ring.Poly, parties)
			for i := 0; i < parties; i++ {
				shares[i] = rkg[i].GenShareRotLeft(sk0Shards[i].Get(), n, crp)
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

			verifyTestVectors(params.decryptorSk0, params.encoder, coeffsWant, receiver, t)
		}

	})
}

func test_CKS(params *testParams.contexts[i], t *testing.T) {

	bfvContext := params.bfvContext
	contextKeys := params.contextKeys
	parties := params.parties
	sk0Shards := params.sk0Shards
	sk1Shards := params.sk1Shards
	encoder := params.encoder
	decryptorSk1 := params.decryptorSk1

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
			p.s0 = sk0Shards[i].Get()
			p.s1 = sk1Shards[i].Get()
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

		if equalslice(coeffs, encoder.DecodeUint(decryptorSk1.DecryptNew(ksCiphertext))) != true {
			t.Errorf("error : CKS")
		}

		P0.KeySwitch(P0.share, ciphertext, ciphertext)

		if equalslice(coeffs, encoder.DecodeUint(decryptorSk1.DecryptNew(ciphertext))) != true {
			t.Errorf("error : CKS in place")
		}

	})
}

func test_PCKS(params *testParams.contexts[i], t *testing.T) {

	bfvContext := params.bfvContext
	contextKeys := params.contextKeys
	parties := params.parties
	sk0Shards := params.sk0Shards
	pk1 := params.pk1
	encoder := params.encoder
	decryptorSk1 := params.decryptorSk1

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
			p.s = sk0Shards[i].Get()
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
		if equalslice(coeffs, encoder.DecodeUint(decryptorSk1.DecryptNew(ciphertextSwitched))) != true {
			t.Errorf("error : PCKS")
		}
	})
}

func test_BOOT(params *testParams.contexts[i], t *testing.T) {

	bfvContext := params.bfvContext
	contextKeys := params.contextKeys
	parties := params.parties
	sk0Shards := params.sk0Shards
	encoder := params.encoder
	decryptorSk0 := params.decryptorSk0
	sk0 := params.sk0

	t.Run(fmt.Sprintf("N=%d/logQ=%d/BOOT", contextKeys.N, contextKeys.ModulusBigint.Value.BitLen()), func(t *testing.T) {

		crpGenerator, err := ring.NewCRPGenerator(nil, contextKeys)
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
			refreshShares[i] = GenRefreshShares(sk0Shards[i], ciphertext, bfvContext, crp, encoder)
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

		plaintextHave := decryptorSk0.DecryptNew(ciphertext)
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
		plaintextHave = decryptorSk0.DecryptNew(ciphertext)
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
*/

func newTestVectors(contextParams *dbfvContext, encryptor *bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {

	coeffsPol := contextParams.bfvContext.ContextT().NewUniformPoly()

	plaintext = contextParams.bfvContext.NewPlaintext()

	err = contextParams.encoder.EncodeUint(coeffsPol.Coeffs[0], plaintext)
	check(t, err)

	ciphertext, err = encryptor.EncryptNew(plaintext)
	check(t, err)

	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(contextParams *dbfvContext, decryptor *bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {

	if bfv.EqualSlice(coeffs, contextParams.encoder.DecodeUint(decryptor.DecryptNew(ciphertext))) != true {
		t.Errorf("decryption error")
	}
}
