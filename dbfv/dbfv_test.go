package dbfv

import (
	"fmt"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"log"
	"math/big"
	"testing"
)

func check(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
	}
}

type dbfvTestContext struct {
	*dbfvContext

	params *bfv.Parameters

	encoder    *bfv.Encoder
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

	contexts map[uint64]*bfv.Parameters
}

var err error
var testParams = new(dbfvTestParameters)

func init() {
	testParams.parties = 3

	testParams.contexts = bfv.DefaultParams
}

func Test_DBFV(t *testing.T) {
	t.Run("PublicKeyGen", testPublicKeyGen)
	t.Run("RelinKeyGen", testRelinKeyGen)
	t.Run("RelinKeyGenNaive", testRelinKeyGenNaive)
	t.Run("KeySwitching", testKeyswitching)
	t.Run("PublicKeySwitching", testPublicKeySwitching)
	t.Run("RotKeyGenRotRows", testRotKeyGenRotRows)
	t.Run("RotKeyGenRotCols", testRotKeyGenRotCols)
	t.Run("Refresh", testRefresh)

}

func genDBFVTestContext(params *bfv.Parameters) (testCtx *dbfvTestContext) {

	testCtx = new(dbfvTestContext)

	testCtx.params = params
	testCtx.dbfvContext = newDbfvContext(params)

	testCtx.encoder = bfv.NewEncoder(params)
	testCtx.evaluator = bfv.NewEvaluator(params)

	kgen := bfv.NewKeyGenerator(params)



	// SecretKeys
	testCtx.sk0Shards = make([]*bfv.SecretKey, testParams.parties)
	testCtx.sk1Shards = make([]*bfv.SecretKey, testParams.parties)
	tmp0 := testCtx.contextQ1P.NewPoly()
	tmp1 := testCtx.contextQ1P.NewPoly()

	for j := uint64(0); j < testParams.parties; j++ {
		testCtx.sk0Shards[j] = kgen.NewSecretKey()
		testCtx.sk1Shards[j] = kgen.NewSecretKey()
		testCtx.contextQ1P.Add(tmp0, testCtx.sk0Shards[j].Get(), tmp0)
		testCtx.contextQ1P.Add(tmp1, testCtx.sk1Shards[j].Get(), tmp1)
	}

	testCtx.sk0 = new(bfv.SecretKey)
	testCtx.sk1 = new(bfv.SecretKey)

	testCtx.sk0.Set(tmp0)
	testCtx.sk1.Set(tmp1)

	// Publickeys
	testCtx.pk0 = kgen.NewPublicKey(testCtx.sk0)
	testCtx.pk1 = kgen.NewPublicKey(testCtx.sk1)

	testCtx.encryptorPk0 = bfv.NewEncryptorFromPk(params, testCtx.pk0)
	testCtx.decryptorSk0 = bfv.NewDecryptor(params, testCtx.sk0)
	testCtx.decryptorSk1 = bfv.NewDecryptor(params, testCtx.sk1)

	return
}

func testPublicKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		decryptorSk0 := testCtx.decryptorSk0

		t.Run(fmt.Sprintf("N=%d/logQ=%d", testCtx.n, testCtx.contextQ1P.ModulusBigint.BitLen()), func(t *testing.T) {

			crpGenerator := ring.NewCRPGenerator(nil, testCtx.contextQ1P)
			crpGenerator.Seed([]byte{})
			crp := crpGenerator.ClockNew()

			type Party struct {
				*CKGProtocol
				s  *ring.Poly
				s1 CKGShare
			}

			ckgParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.CKGProtocol = NewCKGProtocol(parameters)
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
			encryptorTest := bfv.NewEncryptorFromPk(parameters, pk)

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorTest, t)

			verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)
		})
	}
}

func testRelinKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk0 := testCtx.decryptorSk0
		evaluator := testCtx.evaluator

		t.Run(fmt.Sprintf("N=%d/logQ=%d", testCtx.n, testCtx.contextQ1P.ModulusBigint.BitLen()), func(t *testing.T) {

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
				p.RKGProtocol = NewEkgProtocol(parameters)
				p.u = p.RKGProtocol.NewEphemeralKey(1.0 / 3.0)
				p.s = sk0Shards[i].Get()
				p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()
				rkgParties[i] = p
			}

			P0 := rkgParties[0]

			crpGenerator := ring.NewCRPGenerator(nil, testCtx.contextQ1P)
			crpGenerator.Seed([]byte{})
			crp := make([]*ring.Poly, testCtx.beta)

			for i := uint64(0); i < testCtx.beta; i++ {
				crp[i] = crpGenerator.ClockNew()
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

			evk := bfv.NewRelinKey(parameters, 1)
			P0.GenRelinearizationKey(P0.share2, P0.share3, evk)

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			for i := range coeffs {
				coeffs[i] *= coeffs[i]
				coeffs[i] %= testCtx.contextT.Modulus[0]
			}

			ciphertextMul := bfv.NewCiphertext(parameters, ciphertext.Degree()*2)
			evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

			res := bfv.NewCiphertext(parameters, 1)
			evaluator.Relinearize(ciphertextMul, evk, res)

			verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertextMul, t)
		})
	}
}

func testRelinKeyGenNaive(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		evaluator := testCtx.evaluator
		pk0 := testCtx.pk0
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk0 := testCtx.decryptorSk0
		sk0Shards := testCtx.sk0Shards

		t.Run(fmt.Sprintf("N=%d/logQ=%d", testCtx.n, testCtx.contextQ1P.ModulusBigint.BitLen()), func(t *testing.T) {

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
				p.RKGProtocolNaive = NewRKGProtocolNaive(parameters)
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

			evk := bfv.NewRelinKey(parameters, 1)
			P0.GenRelinearizationKey(P0.share2, evk)

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			for i := range coeffs {
				coeffs[i] *= coeffs[i]
				coeffs[i] %= testCtx.contextT.Modulus[0]
			}

			ciphertextMul := bfv.NewCiphertext(parameters, ciphertext.Degree()*2)
			evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

			res := bfv.NewCiphertext(parameters, 1)
			evaluator.Relinearize(ciphertextMul, evk, res)

			verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertextMul, t)
		})
	}
}

func testKeyswitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		sk1Shards := testCtx.sk1Shards
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk1 := testCtx.decryptorSk1

		t.Run(fmt.Sprintf("N=%d/logQ=%d", testCtx.n, testCtx.contextQ1P.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*CKSProtocol
				s0    *ring.Poly
				s1    *ring.Poly
				share CKSShare
			}

			cksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.CKSProtocol = NewCKSProtocol(parameters, 6.36)
				p.s0 = sk0Shards[i].Get()
				p.s1 = sk1Shards[i].Get()
				p.share = p.AllocateShare()
				cksParties[i] = p
			}
			P0 := cksParties[0]

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			// Each party creates its CKSProtocol instance with tmp = si-si'
			for i, p := range cksParties {
				p.GenShare(p.s0, p.s1, ciphertext, p.share)
				if i > 0 {
					P0.AggregateShares(p.share, P0.share, P0.share)
				}
			}

			ksCiphertext := bfv.NewCiphertext(parameters, 1)
			P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

			verifyTestVectors(testCtx, decryptorSk1, coeffs, ksCiphertext, t)

			P0.KeySwitch(P0.share, ciphertext, ciphertext)

			verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertext, t)

		})
	}
}

func testPublicKeySwitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		pk1 := testCtx.pk1
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk1 := testCtx.decryptorSk1

		t.Run(fmt.Sprintf("N=%d/logQ=%d", testCtx.n, testCtx.contextQ1P.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*PCKSProtocol
				s     *ring.Poly
				share PCKSShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.PCKSProtocol = NewPCKSProtocol(parameters, 6.36)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShares()
				pcksParties[i] = p
			}
			P0 := pcksParties[0]

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			ciphertextSwitched := bfv.NewCiphertext(parameters, 1)

			for i, p := range pcksParties {
				p.GenShare(p.s, pk1, ciphertext, p.share)
				if i > 0 {
					P0.AggregateShares(p.share, P0.share, P0.share)
				}
			}

			P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

			verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertextSwitched, t)
		})
	}
}

func testRotKeyGenRotRows(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		testCtx := genDBFVTestContext(parameters)

		evaluator := testCtx.evaluator
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk0 := testCtx.decryptorSk0
		sk0Shards := testCtx.sk0Shards

		t.Run(fmt.Sprintf("N=%d/logQ=%d", testCtx.n, testCtx.contextQ1P.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*RTGProtocol
				s     *ring.Poly
				share RTGShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RTGProtocol = NewRotKGProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShare()
				pcksParties[i] = p
			}
			P0 := pcksParties[0]

			crpGenerator := ring.NewCRPGenerator(nil, testCtx.contextQ1P)
			crpGenerator.Seed([]byte{})
			crp := make([]*ring.Poly, testCtx.beta)

			for i := uint64(0); i < testCtx.beta; i++ {
				crp[i] = crpGenerator.ClockNew()
			}

			for i, p := range pcksParties {
				p.GenShare(bfv.RotationRow, 0, p.s, crp, &p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			rotkey := bfv.NewRotationKeys()
			P0.Finalize(P0.share, crp, rotkey)

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			evaluator.RotateRows(ciphertext, rotkey, ciphertext)

			coeffs = append(coeffs[testCtx.n>>1:], coeffs[:testCtx.n>>1]...)

			verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)

		})
	}

}

func testRotKeyGenRotCols(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		evaluator := testCtx.evaluator
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk0 := testCtx.decryptorSk0
		sk0Shards := testCtx.sk0Shards

		t.Run(fmt.Sprintf("N=%d/logQ=%d", testCtx.n, testCtx.contextQ1P.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*RTGProtocol
				s     *ring.Poly
				share RTGShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RTGProtocol = NewRotKGProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShare()
				pcksParties[i] = p
			}

			P0 := pcksParties[0]

			crpGenerator := ring.NewCRPGenerator(nil, testCtx.contextQ1P)
			crpGenerator.Seed([]byte{})
			crp := make([]*ring.Poly, testCtx.beta)

			for i := uint64(0); i < testCtx.beta; i++ {
				crp[i] = crpGenerator.ClockNew()
			}

			mask := (testCtx.n >> 1) - 1

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			receiver := bfv.NewCiphertext(parameters, ciphertext.Degree())

			for k := uint64(1); k < testCtx.n>>1; k <<= 1 {

				for i, p := range pcksParties {
					p.GenShare(bfv.RotationLeft, k, p.s, crp, &p.share)
					if i > 0 {
						P0.Aggregate(p.share, P0.share, P0.share)
					}
				}

				rotkey := bfv.NewRotationKeys()
				P0.Finalize(P0.share, crp, rotkey)

				evaluator.RotateColumns(ciphertext, k, rotkey, receiver)

				coeffsWant := make([]uint64, testCtx.n)

				for i := uint64(0); i < testCtx.n>>1; i++ {
					coeffsWant[i] = coeffs[(i+k)&mask]
					coeffsWant[i+(testCtx.n>>1)] = coeffs[((i+k)&mask)+(testCtx.n>>1)]
				}

				verifyTestVectors(testCtx, decryptorSk0, coeffsWant, receiver, t)
			}
		})
	}
}

func testRefresh(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		encryptorPk0 := testCtx.encryptorPk0
		sk0Shards := testCtx.sk0Shards
		encoder := testCtx.encoder
		decryptorSk0 := testCtx.decryptorSk0

		kgen := bfv.NewKeyGenerator(parameters)

		rlk := kgen.NewRelinKey(testCtx.sk0, 2)

		t.Run(fmt.Sprintf("N=%d/logQ=%d/Refresh", testCtx.n, testCtx.contextQ1P.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*RefreshProtocol
				s       *ring.Poly
				share   RefreshShare
				ptShare *bfv.Plaintext
			}

			RefreshParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RefreshProtocol = NewRefreshProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShares()
				p.ptShare = bfv.NewPlaintext(parameters)
				RefreshParties[i] = p
			}

			P0 := RefreshParties[0]

			crpGenerator := ring.NewCRPGenerator(nil, testCtx.contextQ1P)
			crpGenerator.Seed([]byte{})
			crp := crpGenerator.ClockNew()

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			maxDepth := 0

			ciphertextTmp := ciphertext.CopyNew().Ciphertext()
			coeffsTmp := make([]uint64, len(coeffs))

			for i := range coeffs {
				coeffsTmp[i] = coeffs[i]
			}

			// Finds the maximum multiplicative depth
			for true {

				testCtx.evaluator.Relinearize(testCtx.evaluator.MulNew(ciphertextTmp, ciphertextTmp), rlk, ciphertextTmp)

				for j := range coeffsTmp {
					coeffsTmp[j] = ring.BRed(coeffsTmp[j], coeffsTmp[j], testCtx.contextT.Modulus[0], testCtx.contextT.GetBredParams()[0])
				}

				if utils.EqualSliceUint64(coeffsTmp, encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertextTmp))) {
					maxDepth++
				} else {
					break
				}
			}

			// Simulated added error of size Q/(T^2) and add it to the fresh ciphertext
			coeffsBigint := make([]*big.Int, testCtx.n)
			testCtx.contextQ1.PolyToBigint(ciphertext.Value()[0], coeffsBigint)

			errorRange := new(big.Int).Set(testCtx.contextQ1.ModulusBigint)
			errorRange.Quo(errorRange, testCtx.contextT.ModulusBigint)
			errorRange.Quo(errorRange, testCtx.contextT.ModulusBigint)

			for i := uint64(0); i < testCtx.n; i++ {
				coeffsBigint[i].Add(coeffsBigint[i], ring.RandInt(errorRange))
			}

			testCtx.contextQ1.SetCoefficientsBigint(coeffsBigint, ciphertext.Value()[0])

			for i, p := range RefreshParties {
				p.GenShares(p.s, ciphertext, crp, p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			// We refresh the ciphertext with the simulated error
			P0.Decrypt(ciphertext, P0.share.RefreshShareDecrypt, P0.ptShare.Value()[0])      // Masked decryption
			P0.Recode(P0.ptShare.Value()[0], P0.ptShare.Value()[0])                          // Masked re-encoding
			P0.Recrypt(P0.ptShare.Value()[0], crp, P0.share.RefreshShareRecrypt, ciphertext) // Masked re-encryption$

			// The refresh also be called all at once with P0.Finalize(ciphertext, crp, P0.share, ctOut)

			// Square the refreshed ciphertext up to the maximum depth-1
			for i := 0; i < maxDepth-1; i++ {

				testCtx.evaluator.Relinearize(testCtx.evaluator.MulNew(ciphertext, ciphertext), rlk, ciphertext)

				for j := range coeffs {
					coeffs[j] = ring.BRed(coeffs[j], coeffs[j], testCtx.contextT.Modulus[0], testCtx.contextT.GetBredParams()[0])
				}
			}

			//Decrypts and compare
			if utils.EqualSliceUint64(coeffs, encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertext))) != true {
				t.Errorf("error : BOOT")
			}
		})
	}
}

func newTestVectors(contextParams *dbfvTestContext, encryptor *bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {
	coeffsPol := contextParams.contextT.NewUniformPoly()
	plaintext = bfv.NewPlaintext(contextParams.params)
	contextParams.encoder.EncodeUint(coeffsPol.Coeffs[0], plaintext)
	ciphertext = encryptor.EncryptNew(plaintext)
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(contextParams *dbfvTestContext, decryptor *bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {
	if utils.EqualSliceUint64(coeffs, contextParams.encoder.DecodeUint(decryptor.DecryptNew(ciphertext))) != true {
		t.Errorf("decryption error")
	}
}

func Test_Marshalling(t *testing.T) {
	params := bfv.DefaultParams[14]

	//verify if the un.marshalling works properly
	dbfvCtx := newDbfvContext(params)
	KeyGenerator := bfv.NewKeyGenerator(params)
	crsGen := ring.NewCRPGenerator([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'}, dbfvCtx.contextQ1P)
	sk := KeyGenerator.NewSecretKey()
	crs := crsGen.ClockNew()
	contextQ := dbfvCtx.contextQ1
	contextPKeys := dbfvCtx.contextP

	Ciphertext := bfv.NewCiphertextRandom(params, 1)

	t.Run(fmt.Sprintf("CPK/N=%d/limbQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {
		keygenProtocol := NewCKGProtocol(params)
		KeyGenShareBefore := keygenProtocol.AllocateShares()
		keygenProtocol.GenShare(sk.Get(), crs, KeyGenShareBefore)
		//now we marshall it
		data, err := KeyGenShareBefore.MarshalBinary()

		fmt.Println("CPK Data (kB) :", len(data)/1000)

		if err != nil {
			log.Fatal("Could not marshal the CKGShare : ", err)

		}

		KeyGenShareAfter := new(CKGShare)
		err = KeyGenShareAfter.UnmarshalBinary(data)
		if err != nil {
			log.Fatal("Could not unmarshal the CKGShare : ", err)

		}

		//comparing the results
		if KeyGenShareBefore.GetDegree() != KeyGenShareAfter.GetDegree() || KeyGenShareBefore.GetLenModuli() != KeyGenShareAfter.GetLenModuli() {
			log.Print("Unmatched degree or moduli length on key gen shares")
			t.Fail()
		}

		for i := 0; i < KeyGenShareBefore.GetLenModuli(); i++ {
			if !utils.EqualSliceUint64(KeyGenShareAfter.Coeffs[i], KeyGenShareBefore.Coeffs[i]) {
				log.Print("Non equal slices in CKGShare")
				t.Fail()
			}

		}
	})

	t.Run(fmt.Sprintf("PCKS/N=%d/limbQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {
		//Check marshalling for the PCKS

		KeySwitchProtocol := NewPCKSProtocol(params, dbfvCtx.params.Sigma)
		SwitchShare := KeySwitchProtocol.AllocateShares()
		pk := KeyGenerator.NewPublicKey(sk)
		KeySwitchProtocol.GenShare(sk.Get(), pk, Ciphertext, SwitchShare)

		data, err := SwitchShare.MarshalBinary()
		fmt.Println("PCKS Data (kB) :", len(data)/1000)
		if err != nil {
			log.Print("Error on PCKSShare marshalling : ", err)
			t.Fail()
		}
		SwitchShareReceiver := new(PCKSShare)
		err = SwitchShareReceiver.UnmarshalBinary(data)
		if err != nil {
			log.Print("Error on PCKSShare unmarshalling : ", err)
			t.Fail()
		}

		for i := 0; i < 2; i++ {
			//compare the shares.
			ringBefore := SwitchShare[i]
			ringAfter := SwitchShareReceiver[i]
			if ringBefore.GetDegree() != ringAfter.GetDegree() {
				log.Print("Error on degree matching")
				t.Fail()
			}
			for d := 0; d < ringAfter.GetLenModuli(); d++ {
				if !utils.EqualSliceUint64(ringAfter.Coeffs[d], ringBefore.Coeffs[d]) {
					log.Print("Non equal slices in PCKSShare")
					t.Fail()
				}
			}

		}
	})

	t.Run(fmt.Sprintf("CKS/N=%d/limbQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {

		//Now for CKSShare ~ its similar to PKSShare
		cksp := NewCKSProtocol(params, dbfvCtx.params.Sigma)
		cksshare := cksp.AllocateShare()
		skIn := KeyGenerator.NewSecretKey()
		skOut := KeyGenerator.NewSecretKey()
		cksp.GenShare(skIn.Get(), skOut.Get(), Ciphertext, cksshare)

		data, err := cksshare.MarshalBinary()
		fmt.Println("CKS Data (kB) :", len(data)/1000)
		if err != nil {
			log.Print("Error on marshalling CKSShare : ", err)
			t.Fail()
		}
		cksshareAfter := new(CKSShare)
		err = cksshareAfter.UnmarshalBinary(data)
		if err != nil {
			log.Print("Error on unmarshalling CKSShare : ", err)
			t.Fail()
		}

		//now compare both shares.

		if cksshare.GetDegree() != cksshareAfter.GetDegree() || cksshare.GetLenModuli() != cksshareAfter.GetLenModuli() {
			log.Print("Unmatched degree or moduli lenght on CKSShare")
			t.Fail()
		}

		for i := 0; i < cksshare.GetLenModuli(); i++ {
			if !utils.EqualSliceUint64(cksshare.Coeffs[i], cksshareAfter.Coeffs[i]) {
				log.Print("Non equal slices in CKSShare")
				t.Fail()
			}

		}
	})

	t.Run(fmt.Sprintf("Refresh/N=%d/limbQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {

		//testing refresh shares
		refreshproto := NewRefreshProtocol(params)
		refreshshare := refreshproto.AllocateShares()
		refreshproto.GenShares(sk.Get(), Ciphertext, crs, refreshshare)

		data, err := refreshshare.MarshalBinary()
		fmt.Println("Refresh Data (kB) :", len(data)/1000)
		if err != nil {
			log.Fatal("Could not marshal RefreshShare", err)
		}
		resRefreshShare := new(RefreshShare)
		err = resRefreshShare.UnmarshalBinary(data)

		if err != nil {
			log.Fatal("Could not unmarshal RefreshShare", err)
		}
		for i, r := range refreshshare.RefreshShareDecrypt.Coeffs {
			if !utils.EqualSliceUint64(resRefreshShare.RefreshShareDecrypt.Coeffs[i], r) {
				log.Fatal("Resulting of marshalling not the same as original : RefreshShare")
			}

		}
		for i, r := range refreshshare.RefreshShareRecrypt.Coeffs {
			if !utils.EqualSliceUint64(resRefreshShare.RefreshShareRecrypt.Coeffs[i], r) {
				log.Fatal("Resulting of marshalling not the same as original : RefreshShare")
			}

		}
	})

	t.Run(fmt.Sprintf("RTG/N=%d/limbQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {

		//check RTGShare
		crpGenerator := ring.NewCRPGenerator(nil, dbfvCtx.contextQ1P)
		modulus := (dbfvCtx.contextQ1.Modulus)
		crp := make([]*ring.Poly, len(modulus))
		for j := 0; j < len(modulus); j++ {
			crp[j] = crpGenerator.ClockNew() //make([]*ring.Poly, bitLog)

		}

		rotProto := NewRotKGProtocol(params)
		rtgShare := rotProto.AllocateShare()
		rotProto.GenShare(1, 64, sk.Get(), crp, &rtgShare)

		data, err := rtgShare.MarshalBinary()
		fmt.Println("RTG Data (kB) :", len(data)/1000)
		if err != nil {
			log.Fatal("could not marshal RTGshare :", err)
		}

		resRTGShare := new(RTGShare)
		err = resRTGShare.UnmarshalBinary(data)
		if err != nil {
			log.Fatal("Could not unmarshal RTGShare: ", err)
		}

		if resRTGShare.Type != rtgShare.Type || resRTGShare.K != rtgShare.K || len(resRTGShare.Value) != len(rtgShare.Value) {
			log.Fatal("result after marshalling is not the same as before marshalling for RTGSahre")
		}

		for i, val := range rtgShare.Value {
			ring1 := val
			ring2 := resRTGShare.Value[i]
			if len(ring1.Coeffs) != len(ring2.Coeffs) {
				log.Fatal("result after marshalling is not the same as before marshalling for RTGSahre")
			}

			for j, elem := range ring1.Coeffs {
				if !utils.EqualSliceUint64(ring2.Coeffs[j], elem) {
					log.Fatal("result after marshalling is not the same as before marshalling for RTGSahre")

				}
			}

		}
	})

}

func Test_Relin_Marshalling(t *testing.T) {
	params := bfv.DefaultParams[14]

	dbfvCtx := newDbfvContext(params)
	contextQ := dbfvCtx.contextQ1
	contextPKeys := dbfvCtx.contextP
	modulus := dbfvCtx.contextQ1.Modulus

	crpGenerator := ring.NewCRPGenerator(nil, dbfvCtx.contextQ1P)

	crp := make([]*ring.Poly, len(modulus))
	for j := 0; j < len(modulus); j++ {
		crp[j] = crpGenerator.ClockNew() //make([]*ring.Poly, bitLog)
		//for u := uint64(0); u < bitLog; u++ {
		//	crp[j][u] = crpGenerator.ClockNew()
		//}
	}

	t.Run(fmt.Sprintf("RLKG/N=%d/limbQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {

		rlk := NewEkgProtocol(params)
		u := rlk.NewEphemeralKey(1 / 3.0)
		sk := bfv.NewKeyGenerator(params).NewSecretKey()
		log.Print("Starting to test marshalling for share one")

		r1, r2, r3 := rlk.AllocateShares()
		rlk.GenShareRoundOne(u, sk.Get(), crp, r1)
		data, err := r1.MarshalBinary()
		fmt.Println("RKG R1 Data (kB) :", len(data)/1000)
		if err != nil {
			log.Print("Error in marshalling round 1 key : ", err)
			t.Fail()
		}

		r1After := new(RKGShareRoundOne)
		err = r1After.UnmarshalBinary(data)
		if err != nil {
			log.Print("Error in unmarshalling round 1 key : ", err)
			t.Fail()
		}

		log.Print("Now comparing keys for round 1 ")

		for i := 0; i < (len(r1)); i++ {
			a := r1[i]
			b := (*r1After)[i]
			for k := 0; k < a.GetLenModuli(); k++ {
				if !utils.EqualSliceUint64(a.Coeffs[k], b.Coeffs[k]) {
					log.Print("Error : coeffs of rings do not match")
					t.Fail()
				}
			}
		}

		log.Print("Sucess : relin key round 1 ok ")

		log.Print("Starting to test marshalling for share two")
		rlk.GenShareRoundTwo(r1, sk.Get(), crp, r2)

		data, err = r2.MarshalBinary()
		fmt.Println("RKG R2 Data (kB) :", len(data)/1000)
		if err != nil {
			log.Print("Error on marshalling relin key round 2 : ", err)
			t.Fail()
		}

		r2After := new(RKGShareRoundTwo)
		err = r2After.UnmarshalBinary(data)
		if err != nil {
			log.Print("Error on unmarshalling relin key round 2 : ", err)
			t.Fail()
		}

		log.Print("Now comparing keys for round 2 ")

		for i := 0; i < (len(r2)); i++ {
			for idx := 0; idx < 2; idx++ {
				a := r2[i][idx]
				b := (*r2After)[i][idx]
				for k := 0; k < a.GetLenModuli(); k++ {
					if !utils.EqualSliceUint64(a.Coeffs[k], b.Coeffs[k]) {
						log.Print("Error : coeffs of rings do not match")
						t.Fail()
					}
				}
			}

		}

		log.Print("Success : reling key round 2 ok ")

		log.Print("Starting to test marshalling for share three")

		rlk.GenShareRoundThree(r2, u, sk.Get(), r3)

		data, err = r3.MarshalBinary()
		fmt.Println("RKG R3 Data (kB) :", len(data)/1000)
		if err != nil {
			log.Print("Error in marshalling round 3 key : ", err)
			t.Fail()
		}

		r3After := new(RKGShareRoundThree)
		err = r3After.UnmarshalBinary(data)
		if err != nil {
			log.Print("Error in unmarshalling round 3 key : ", err)
			t.Fail()
		}

		log.Print("Now comparing keys for round 3 ")

		for i := 0; i < (len(r3)); i++ {
			a := r3[i]
			b := (*r3After)[i]
			for k := 0; k < a.GetLenModuli(); k++ {
				if !utils.EqualSliceUint64(a.Coeffs[k], b.Coeffs[k]) {
					log.Print("Error : coeffs of rings do not match")
					t.Fail()
				}
			}

		}

		log.Print("Success : relin key for round 3 ok ")

		log.Print("All marshalling is passed for relin keys")
	})

}
