package dbfv

import (
	"fmt"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math/big"
	"testing"
)

func check(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
	}
}

type dbfvContext struct {
	bfvContext *bfv.BfvContext
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

	contexts []bfv.Parameters
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

func genDBFVContext(contextParameters *bfv.Parameters) (params *dbfvContext) {

	params = new(dbfvContext)

	params.bfvContext = bfv.NewBfvContextWithParam(contextParameters)

	params.encoder = params.bfvContext.NewEncoder()
	params.evaluator = params.bfvContext.NewEvaluator()

	kgen := params.bfvContext.NewKeyGenerator()

	// SecretKeys
	params.sk0Shards = make([]*bfv.SecretKey, testParams.parties)
	params.sk1Shards = make([]*bfv.SecretKey, testParams.parties)
	tmp0 := params.bfvContext.ContextKeys().NewPoly()
	tmp1 := params.bfvContext.ContextKeys().NewPoly()

	for j := uint64(0); j < testParams.parties; j++ {
		params.sk0Shards[j] = kgen.NewSecretKey()
		params.sk1Shards[j] = kgen.NewSecretKey()
		params.bfvContext.ContextKeys().Add(tmp0, params.sk0Shards[j].Get(), tmp0)
		params.bfvContext.ContextKeys().Add(tmp1, params.sk1Shards[j].Get(), tmp1)
	}

	params.sk0 = new(bfv.SecretKey)
	params.sk1 = new(bfv.SecretKey)

	params.sk0.Set(tmp0)
	params.sk1.Set(tmp1)

	// Publickeys
	params.pk0 = kgen.NewPublicKey(params.sk0)
	params.pk1 = kgen.NewPublicKey(params.sk1)

	params.encryptorPk0 = params.bfvContext.NewEncryptorFromPk(params.pk0)
	params.decryptorSk0 = params.bfvContext.NewDecryptor(params.sk0)
	params.decryptorSk1 = params.bfvContext.NewDecryptor(params.sk1)

	return
}

func testPublicKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		sk0Shards := params.sk0Shards
		decryptorSk0 := params.decryptorSk0

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.BitLen()), func(t *testing.T) {

			crpGenerator := ring.NewCRPGenerator(nil, contextKeys)
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
			encryptorTest := bfvContext.NewEncryptorFromPk(pk)

			coeffs, _, ciphertext := newTestVectors(params, encryptorTest, t)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertext, t)
		})
	}
}

func testRelinKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		sk0Shards := params.sk0Shards
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		evaluator := params.evaluator

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.BitLen()), func(t *testing.T) {

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
				p.u = p.RKGProtocol.NewEphemeralKey(1.0 / 3.0)
				p.s = sk0Shards[i].Get()
				p.share1, p.share2, p.share3 = p.RKGProtocol.AllocateShares()
				rkgParties[i] = p
			}

			P0 := rkgParties[0]

			crpGenerator := ring.NewCRPGenerator(nil, bfvContext.ContextKeys())
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
			evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

			res := bfvContext.NewCiphertext(1)
			evaluator.Relinearize(ciphertextMul, evk, res)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertextMul, t)
		})
	}
}

func testRelinKeyGenNaive(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		evaluator := params.evaluator
		pk0 := params.pk0
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.BitLen()), func(t *testing.T) {

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
				p.RKGProtocolNaive = NewRKGProtocolNaive(bfvContext)
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

			evk := bfvContext.NewRelinKeyEmpty(1)
			P0.GenRelinearizationKey(P0.share2, evk)

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, t)

			for i := range coeffs {
				coeffs[i] *= coeffs[i]
				coeffs[i] %= bfvContext.ContextT().Modulus[0]
			}

			ciphertextMul := bfvContext.NewCiphertext(ciphertext.Degree() * 2)
			evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

			res := bfvContext.NewCiphertext(1)
			evaluator.Relinearize(ciphertextMul, evk, res)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertextMul, t)
		})
	}
}

func testKeyswitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		sk0Shards := params.sk0Shards
		sk1Shards := params.sk1Shards
		encryptorPk0 := params.encryptorPk0
		decryptorSk1 := params.decryptorSk1

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*CKSProtocol
				s0    *ring.Poly
				s1    *ring.Poly
				share CKSShare
			}

			cksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.CKSProtocol = NewCKSProtocol(bfvContext, 6.36)
				p.s0 = sk0Shards[i].Get()
				p.s1 = sk1Shards[i].Get()
				p.share = p.AllocateShare()
				cksParties[i] = p
			}
			P0 := cksParties[0]

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, t)

			// Each party creates its CKSProtocol instance with tmp = si-si'
			for i, p := range cksParties {
				p.GenShare(p.s0, p.s1, ciphertext, p.share)
				if i > 0 {
					P0.AggregateShares(p.share, P0.share, P0.share)
				}
			}

			ksCiphertext := bfvContext.NewCiphertext(1)
			P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

			verifyTestVectors(params, decryptorSk1, coeffs, ksCiphertext, t)

			P0.KeySwitch(P0.share, ciphertext, ciphertext)

			verifyTestVectors(params, decryptorSk1, coeffs, ciphertext, t)

		})
	}
}

func testPublicKeySwitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		sk0Shards := params.sk0Shards
		pk1 := params.pk1
		encryptorPk0 := params.encryptorPk0
		decryptorSk1 := params.decryptorSk1

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*PCKSProtocol
				s     *ring.Poly
				share PCKSShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.PCKSProtocol = NewPCKSProtocol(bfvContext, 6.36)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShares()
				pcksParties[i] = p
			}
			P0 := pcksParties[0]

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, t)

			ciphertextSwitched := bfvContext.NewCiphertext(1)

			for i, p := range pcksParties {
				p.GenShare(p.s, pk1, ciphertext, p.share)
				if i > 0 {
					P0.AggregateShares(p.share, P0.share, P0.share)
				}
			}

			P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

			verifyTestVectors(params, decryptorSk1, coeffs, ciphertextSwitched, t)
		})
	}
}

func testRotKeyGenRotRows(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*RTGProtocol
				s     *ring.Poly
				share RTGShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RTGProtocol = NewRotKGProtocol(bfvContext)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShare()
				pcksParties[i] = p
			}
			P0 := pcksParties[0]

			crpGenerator := ring.NewCRPGenerator(nil, bfvContext.ContextKeys())
			crpGenerator.Seed([]byte{})
			crp := make([]*ring.Poly, bfvContext.Beta())

			for i := uint64(0); i < bfvContext.Beta(); i++ {
				crp[i] = crpGenerator.Clock()
			}

			for i, p := range pcksParties {
				p.GenShare(bfv.RotationRow, 0, p.s, crp, &p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			rotkey := bfvContext.NewRotationKeys()
			P0.Finalize(P0.share, crp, rotkey)

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, t)

			evaluator.RotateRows(ciphertext, rotkey, ciphertext)

			coeffs = append(coeffs[contextKeys.N>>1:], coeffs[:contextKeys.N>>1]...)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertext, t)

		})
	}

}

func testRotKeyGenRotCols(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := bfvContext.ContextKeys()
		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(fmt.Sprintf("N=%d/logQ=%d", contextKeys.N, contextKeys.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*RTGProtocol
				s     *ring.Poly
				share RTGShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RTGProtocol = NewRotKGProtocol(bfvContext)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShare()
				pcksParties[i] = p
			}

			P0 := pcksParties[0]

			crpGenerator := ring.NewCRPGenerator(nil, contextKeys)
			crpGenerator.Seed([]byte{})
			crp := make([]*ring.Poly, bfvContext.Beta())

			for i := uint64(0); i < bfvContext.Beta(); i++ {
				crp[i] = crpGenerator.Clock()
			}

			mask := (contextKeys.N >> 1) - 1

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, t)

			receiver := bfvContext.NewCiphertext(ciphertext.Degree())

			for k := uint64(1); k < contextKeys.N>>1; k <<= 1 {

				for i, p := range pcksParties {
					p.GenShare(bfv.RotationLeft, k, p.s, crp, &p.share)
					if i > 0 {
						P0.Aggregate(p.share, P0.share, P0.share)
					}
				}

				rotkey := bfvContext.NewRotationKeys()
				P0.Finalize(P0.share, crp, rotkey)

				evaluator.RotateColumns(ciphertext, k, rotkey, receiver)

				coeffsWant := make([]uint64, contextKeys.N)

				for i := uint64(0); i < contextKeys.N>>1; i++ {
					coeffsWant[i] = coeffs[(i+k)&mask]
					coeffsWant[i+(contextKeys.N>>1)] = coeffs[((i+k)&mask)+(contextKeys.N>>1)]
				}

				verifyTestVectors(params, decryptorSk0, coeffsWant, receiver, t)
			}
		})
	}
}

func testRefresh(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		params := genDBFVContext(&parameters)

		bfvContext := params.bfvContext
		contextKeys := params.bfvContext.ContextKeys()
		encryptorPk0 := params.encryptorPk0
		sk0Shards := params.sk0Shards
		encoder := params.encoder
		decryptorSk0 := params.decryptorSk0

		kgen := bfvContext.NewKeyGenerator()

		rlk := kgen.NewRelinKey(params.sk0, 2)

		t.Run(fmt.Sprintf("N=%d/logQ=%d/Refresh", contextKeys.N, contextKeys.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*RefreshProtocol
				s       *ring.Poly
				share   RefreshShare
				ptShare *bfv.Plaintext
			}

			RefreshParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RefreshProtocol = NewRefreshProtocol(bfvContext)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShares()
				p.ptShare = bfvContext.NewPlaintext()
				RefreshParties[i] = p
			}

			P0 := RefreshParties[0]

			crpGenerator := ring.NewCRPGenerator(nil, bfvContext.ContextKeys())
			crpGenerator.Seed([]byte{})
			crp := crpGenerator.Clock()

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, t)

			maxDepth := 0

			ciphertextTmp := ciphertext.CopyNew().Ciphertext()
			coeffsTmp := make([]uint64, len(coeffs))

			for i := range coeffs {
				coeffsTmp[i] = coeffs[i]
			}

			// Finds the maximum multiplicative depth
			for true {

				params.evaluator.Relinearize(params.evaluator.MulNew(ciphertextTmp, ciphertextTmp), rlk, ciphertextTmp)

				for j := range coeffsTmp {
					coeffsTmp[j] = ring.BRed(coeffsTmp[j], coeffsTmp[j], bfvContext.ContextT().Modulus[0], bfvContext.ContextT().GetBredParams()[0])
				}

				if utils.EqualSliceUint64(coeffsTmp, encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertextTmp))) {
					maxDepth += 1
				} else {
					break
				}
			}

			// Simulated added error of size Q/(T^2) and add it to the fresh ciphertext
			coeffs_bigint := make([]*big.Int, bfvContext.N())
			bfvContext.ContextQ().PolyToBigint(ciphertext.Value()[0], coeffs_bigint)

			error_range := new(big.Int).Set(bfvContext.ContextQ().ModulusBigint)
			error_range.Quo(error_range, bfvContext.ContextT().ModulusBigint)
			error_range.Quo(error_range, bfvContext.ContextT().ModulusBigint)

			for i := uint64(0); i < bfvContext.N(); i++ {
				coeffs_bigint[i].Add(coeffs_bigint[i], ring.RandInt(error_range))
			}

			bfvContext.ContextQ().SetCoefficientsBigint(coeffs_bigint, ciphertext.Value()[0])

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

				params.evaluator.Relinearize(params.evaluator.MulNew(ciphertext, ciphertext), rlk, ciphertext)

				for j := range coeffs {
					coeffs[j] = ring.BRed(coeffs[j], coeffs[j], bfvContext.ContextT().Modulus[0], bfvContext.ContextT().GetBredParams()[0])
				}
			}

			//Decrypts and compare
			if utils.EqualSliceUint64(coeffs, encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertext))) != true {
				t.Errorf("error : BOOT")
			}
		})
	}
}

func newTestVectors(contextParams *dbfvContext, encryptor *bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {
	coeffsPol := contextParams.bfvContext.ContextT().NewUniformPoly()
	plaintext = contextParams.bfvContext.NewPlaintext()
	contextParams.encoder.EncodeUint(coeffsPol.Coeffs[0], plaintext)
	ciphertext = encryptor.EncryptNew(plaintext)
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(contextParams *dbfvContext, decryptor *bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {
	if utils.EqualSliceUint64(coeffs, contextParams.encoder.DecodeUint(decryptor.DecryptNew(ciphertext))) != true {
		t.Errorf("decryption error")
	}
}
