package dbfv

import (
	"flag"
	"fmt"
	"log"
	"math/big"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var parties int = 3

func testString(opname string, parties int, params *bfv.Parameters) string {
	return fmt.Sprintf("%sparties=%d/LogN=%d/logQ=%d", opname, parties, params.LogN(), params.LogQP())
}

type testContext struct {
	params *bfv.Parameters

	dbfvContext *dbfvContext

	prng utils.PRNG

	encoder bfv.Encoder

	sk0Shards []*bfv.SecretKey
	sk0       *bfv.SecretKey

	sk1       *bfv.SecretKey
	sk1Shards []*bfv.SecretKey

	pk0 *bfv.PublicKey
	pk1 *bfv.PublicKey

	encryptorPk0 bfv.Encryptor
	decryptorSk0 bfv.Decryptor
	decryptorSk1 bfv.Decryptor
	evaluator    bfv.Evaluator
}

func Test_DBFV(t *testing.T) {

	var err error

	var defaultParams = bfv.DefaultParams[bfv.PN12QP109 : bfv.PN12QP109+4] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = bfv.DefaultParams[bfv.PN12QP109 : bfv.PN12QP109+2] // the short test runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = bfv.DefaultParams // the long test suite runs for all default parameters
	}
	for _, p := range defaultParams {

		var testCtx *testContext
		if testCtx, err = gentestContext(p); err != nil {
			panic(err)
		}

		testPublicKeyGen(testCtx, t)
		testRelinKeyGen(testCtx, t)
		testKeyswitching(testCtx, t)
		testPublicKeySwitching(testCtx, t)
		testRotKeyGenRotRows(testCtx, t)
		testRotKeyGenRotCols(testCtx, t)
		testRefresh(testCtx, t)
		testRefreshAndPermutation(testCtx, t)
		testMarshalling(testCtx, t)
		testMarshallingRelin(testCtx, t)
	}
}

func gentestContext(defaultParams *bfv.Parameters) (testCtx *testContext, err error) {

	testCtx = new(testContext)

	testCtx.params = defaultParams.Copy()

	testCtx.dbfvContext = newDbfvContext(testCtx.params)

	if testCtx.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testCtx.encoder = bfv.NewEncoder(testCtx.params)
	testCtx.evaluator = bfv.NewEvaluator(testCtx.params, bfv.EvaluationKey{})

	kgen := bfv.NewKeyGenerator(testCtx.params)

	// SecretKeys
	testCtx.sk0Shards = make([]*bfv.SecretKey, parties)
	testCtx.sk1Shards = make([]*bfv.SecretKey, parties)
	tmp0 := testCtx.dbfvContext.ringQP.NewPoly()
	tmp1 := testCtx.dbfvContext.ringQP.NewPoly()

	for j := 0; j < parties; j++ {
		testCtx.sk0Shards[j] = kgen.GenSecretKey()
		testCtx.sk1Shards[j] = kgen.GenSecretKey()
		testCtx.dbfvContext.ringQP.Add(tmp0, testCtx.sk0Shards[j].Value, tmp0)
		testCtx.dbfvContext.ringQP.Add(tmp1, testCtx.sk1Shards[j].Value, tmp1)
	}

	testCtx.sk0 = bfv.NewSecretKey(testCtx.params)
	testCtx.sk1 = bfv.NewSecretKey(testCtx.params)

	testCtx.sk0.Value.Copy(tmp0)
	testCtx.sk1.Value.Copy(tmp1)

	// Publickeys
	testCtx.pk0 = kgen.GenPublicKey(testCtx.sk0)
	testCtx.pk1 = kgen.GenPublicKey(testCtx.sk1)

	testCtx.encryptorPk0 = bfv.NewEncryptorFromPk(testCtx.params, testCtx.pk0)
	testCtx.decryptorSk0 = bfv.NewDecryptor(testCtx.params, testCtx.sk0)
	testCtx.decryptorSk1 = bfv.NewDecryptor(testCtx.params, testCtx.sk1)

	return
}

func testPublicKeyGen(testCtx *testContext, t *testing.T) {

	sk0Shards := testCtx.sk0Shards
	decryptorSk0 := testCtx.decryptorSk0

	t.Run(testString("PublicKeyGen/", parties, testCtx.params), func(t *testing.T) {

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)
		crp := crpGenerator.ReadNew()

		type Party struct {
			*CKGProtocol
			s  *rlwe.SecretKey
			s1 *drlwe.CKGShare
		}

		ckgParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.CKGProtocol = NewCKGProtocol(testCtx.params)
			p.s = &sk0Shards[i].SecretKey
			p.s1 = p.AllocateShares()
			ckgParties[i] = p
		}
		P0 := ckgParties[0]

		// Checks that bfv.CKGProtocol complies to the drlwe.CollectivePublicKeyGenerator interface
		var _ drlwe.CollectivePublicKeyGenerator = P0.CKGProtocol

		// Each party creates a new CKGProtocol instance
		for i, p := range ckgParties {
			p.GenShare(p.s, crp, p.s1)
			if i > 0 {
				P0.AggregateShares(p.s1, P0.s1, P0.s1)
			}
		}

		pk := bfv.NewPublicKey(testCtx.params)
		P0.GenBFVPublicKey(P0.s1, crp, pk)

		// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
		encryptorTest := bfv.NewEncryptorFromPk(testCtx.params, pk)

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorTest, t)

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)
	})
}

func testRelinKeyGen(testCtx *testContext, t *testing.T) {

	sk0Shards := testCtx.sk0Shards
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0

	t.Run(testString("RelinKeyGen/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RKGProtocol
			ephSk  *rlwe.SecretKey
			sk     *rlwe.SecretKey
			share1 *drlwe.RKGShare
			share2 *drlwe.RKGShare
		}

		rkgParties := make([]*Party, parties)

		for i := range rkgParties {
			p := new(Party)
			p.RKGProtocol = NewRKGProtocol(testCtx.params)
			p.sk = &sk0Shards[i].SecretKey
			p.ephSk, p.share1, p.share2 = p.AllocateShares()
			rkgParties[i] = p
		}

		P0 := rkgParties[0]

		// checks that bfv.RKGProtocol complies to the drlwe.RelinearizationKeyGenerator interface
		var _ drlwe.RelinearizationKeyGenerator = P0.RKGProtocol

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := 0; i < testCtx.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		// ROUND 1
		for i, p := range rkgParties {
			p.GenShareRoundOne(p.sk, crp, p.ephSk, p.share1)
			if i > 0 {
				P0.AggregateShares(p.share1, P0.share1, P0.share1)
			}
		}

		//ROUND 2
		for i, p := range rkgParties {
			p.GenShareRoundTwo(p.ephSk, p.sk, P0.share1, crp, p.share2)
			if i > 0 {
				P0.AggregateShares(p.share2, P0.share2, P0.share2)
			}
		}

		evk := bfv.NewRelinearizationKey(testCtx.params, 1)
		P0.GenBFVRelinearizationKey(P0.share1, P0.share2, evk)

		evaluator := testCtx.evaluator.WithKey(bfv.EvaluationKey{Rlk: evk, Rtks: nil})

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)
		for i := range coeffs {
			coeffs[i] *= coeffs[i]
			coeffs[i] %= testCtx.dbfvContext.ringT.Modulus[0]
		}

		ciphertextMul := bfv.NewCiphertext(testCtx.params, ciphertext.Degree()*2)
		evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

		res := bfv.NewCiphertext(testCtx.params, 1)
		evaluator.Relinearize(ciphertextMul, res)

		verifyTestVectors(testCtx, decryptorSk0, coeffs, res, t)
	})

}

func testKeyswitching(testCtx *testContext, t *testing.T) {

	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk1 := testCtx.decryptorSk1

	t.Run(testString("Keyswitching/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*CKSProtocol
			s0    *ring.Poly
			s1    *ring.Poly
			share CKSShare
		}

		cksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.CKSProtocol = NewCKSProtocol(testCtx.params, 6.36)
			p.s0 = sk0Shards[i].Value
			p.s1 = sk1Shards[i].Value
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

		ksCiphertext := bfv.NewCiphertext(testCtx.params, 1)
		P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

		verifyTestVectors(testCtx, decryptorSk1, coeffs, ksCiphertext, t)

		P0.KeySwitch(P0.share, ciphertext, ciphertext)

		verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertext, t)

	})
}

func testPublicKeySwitching(testCtx *testContext, t *testing.T) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk1 := testCtx.decryptorSk1

	t.Run(testString("PublicKeySwitching/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*PCKSProtocol
			s     *ring.Poly
			share PCKSShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.PCKSProtocol = NewPCKSProtocol(testCtx.params, 6.36)
			p.s = sk0Shards[i].Value
			p.share = p.AllocateShares()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		ciphertextSwitched := bfv.NewCiphertext(testCtx.params, 1)

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

func testRotKeyGenRotRows(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	t.Run(testString("RotKeyGenRotRows/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *rlwe.SecretKey
			share *drlwe.RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(testCtx.params)
			p.s = &sk0Shards[i].SecretKey
			p.share = p.AllocateShares()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		// Checks that bfv.RTGProtocol complies to the drlwe.RotationKeyGenerator interface
		var _ drlwe.RotationKeyGenerator = P0.RTGProtocol

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := 0; i < testCtx.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		galEl := testCtx.params.GaloisElementForRowRotation()
		rotKeySet := bfv.NewRotationKeySet(testCtx.params, []uint64{galEl})

		for i, p := range pcksParties {
			p.GenShare(p.s, galEl, crp, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		P0.GenRotationKey(P0.share, crp, rotKeySet.Keys[galEl])

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		evaluator := testCtx.evaluator.WithKey(bfv.EvaluationKey{Rlk: nil, Rtks: rotKeySet})
		result := evaluator.RotateRowsNew(ciphertext)
		coeffsWant := append(coeffs[testCtx.params.N()>>1:], coeffs[:testCtx.params.N()>>1]...)

		verifyTestVectors(testCtx, decryptorSk0, coeffsWant, result, t)

	})
}

func testRotKeyGenRotCols(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	t.Run(testString("RotKeyGenRotCols/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *rlwe.SecretKey
			share *drlwe.RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(testCtx.params)
			p.s = &sk0Shards[i].SecretKey
			p.share = p.AllocateShares()
			pcksParties[i] = p
		}

		P0 := pcksParties[0]

		// Checks that bfv.RTGProtocol complies to the drlwe.RotationKeyGenerator interface
		var _ drlwe.RotationKeyGenerator = P0.RTGProtocol

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := 0; i < testCtx.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		galEls := testCtx.params.GaloisElementsForRowInnerSum()
		rotKeySet := bfv.NewRotationKeySet(testCtx.params, galEls)

		for _, galEl := range galEls {

			for i, p := range pcksParties {
				p.GenShare(p.s, galEl, crp, p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			P0.GenRotationKey(P0.share, crp, rotKeySet.Keys[galEl])
		}

		evaluator := testCtx.evaluator.WithKey(bfv.EvaluationKey{Rlk: nil, Rtks: rotKeySet})
		for k := 1; k < testCtx.params.N()>>1; k <<= 1 {

			result := evaluator.RotateColumnsNew(ciphertext, int(k))
			coeffsWant := utils.RotateUint64Slots(coeffs, int(k))
			verifyTestVectors(testCtx, decryptorSk0, coeffsWant, result, t)
		}
	})
}

func testRefresh(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	sk0Shards := testCtx.sk0Shards
	encoder := testCtx.encoder
	decryptorSk0 := testCtx.decryptorSk0

	kgen := bfv.NewKeyGenerator(testCtx.params)

	rlk := kgen.GenRelinearizationKey(testCtx.sk0, 2)

	t.Run(testString("Refresh/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RefreshProtocol
			s       *ring.Poly
			share   RefreshShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RefreshProtocol = NewRefreshProtocol(testCtx.params)
			p.s = sk0Shards[i].Value
			p.share = p.AllocateShares()
			p.ptShare = bfv.NewPlaintext(testCtx.params)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)
		crp := crpGenerator.ReadNew()

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		maxDepth := 0

		ciphertextTmp := ciphertext.CopyNew().Ciphertext()
		coeffsTmp := make([]uint64, len(coeffs))

		copy(coeffsTmp, coeffs)

		evaluator := testCtx.evaluator.WithKey(bfv.EvaluationKey{Rlk: rlk, Rtks: nil})
		// Finds the maximum multiplicative depth
		for {

			evaluator.Relinearize(testCtx.evaluator.MulNew(ciphertextTmp, ciphertextTmp), ciphertextTmp)

			for j := range coeffsTmp {
				coeffsTmp[j] = ring.BRed(coeffsTmp[j], coeffsTmp[j], testCtx.dbfvContext.ringT.Modulus[0], testCtx.dbfvContext.ringT.BredParams[0])
			}

			if utils.EqualSliceUint64(coeffsTmp, encoder.DecodeUintNew(decryptorSk0.DecryptNew(ciphertextTmp))) {
				maxDepth++
			} else {
				break
			}
		}

		// Simulated added error of size Q/(T^2) and add it to the fresh ciphertext
		coeffsBigint := make([]*big.Int, testCtx.params.N())
		testCtx.dbfvContext.ringQ.PolyToBigint(ciphertext.Value()[0], coeffsBigint)

		errorRange := new(big.Int).Set(testCtx.dbfvContext.ringQ.ModulusBigint)
		errorRange.Quo(errorRange, testCtx.dbfvContext.ringT.ModulusBigint)
		errorRange.Quo(errorRange, testCtx.dbfvContext.ringT.ModulusBigint)

		for i := 0; i < testCtx.params.N(); i++ {
			coeffsBigint[i].Add(coeffsBigint[i], ring.RandInt(errorRange))
		}

		testCtx.dbfvContext.ringQ.SetCoefficientsBigint(coeffsBigint, ciphertext.Value()[0])

		for i, p := range RefreshParties {
			p.GenShares(p.s, ciphertext, crp, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		P0.Finalize(ciphertext, crp, P0.share, ciphertext)

		// Square the refreshed ciphertext up to the maximum depth-1
		for i := 0; i < maxDepth-1; i++ {

			evaluator.Relinearize(testCtx.evaluator.MulNew(ciphertext, ciphertext), ciphertext)

			for j := range coeffs {
				coeffs[j] = ring.BRed(coeffs[j], coeffs[j], testCtx.dbfvContext.ringT.Modulus[0], testCtx.dbfvContext.ringT.BredParams[0])
			}
		}

		//Decrypts and compare
		require.True(t, utils.EqualSliceUint64(coeffs, encoder.DecodeUintNew(decryptorSk0.DecryptNew(ciphertext))))
	})
}

func testRefreshAndPermutation(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	sk0Shards := testCtx.sk0Shards
	encoder := testCtx.encoder
	decryptorSk0 := testCtx.decryptorSk0

	t.Run(testString("RefreshAndPermutation/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*PermuteProtocol
			s       *ring.Poly
			share   RefreshShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.PermuteProtocol = NewPermuteProtocol(testCtx.params)
			p.s = sk0Shards[i].Value
			p.share = p.AllocateShares()
			p.ptShare = bfv.NewPlaintext(testCtx.params)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)
		crp := crpGenerator.ReadNew()

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		permutation := make([]uint64, len(coeffs))

		for i := range permutation {
			permutation[i] = ring.RandUniform(testCtx.prng, uint64(testCtx.params.N()), uint64(testCtx.params.N()-1))
		}

		for i, p := range RefreshParties {
			p.GenShares(p.s, ciphertext, crp, permutation, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		// We refresh the ciphertext with the simulated error
		P0.Decrypt(ciphertext, P0.share.RefreshShareDecrypt, P0.ptShare.Value()[0])      // Masked decryption
		P0.Permute(P0.ptShare.Value()[0], permutation, P0.ptShare.Value()[0])            // Masked re-encoding
		P0.Recrypt(P0.ptShare.Value()[0], crp, P0.share.RefreshShareRecrypt, ciphertext) // Masked re-encryption$

		// The refresh also be called all at once with P0.Finalize(ciphertext, crp, P0.share, ctOut)

		coeffsPermute := make([]uint64, len(coeffs))

		for i := range coeffs {
			coeffsPermute[i] = coeffs[permutation[i]]
		}

		coeffsHave := encoder.DecodeUintNew(decryptorSk0.DecryptNew(ciphertext))

		//Decrypts and compare
		require.True(t, utils.EqualSliceUint64(coeffsPermute, coeffsHave))
	})
}

func newTestVectors(testCtx *testContext, encryptor bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {

	uniformSampler := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringT)

	coeffsPol := uniformSampler.ReadNew()
	plaintext = bfv.NewPlaintext(testCtx.params)
	testCtx.encoder.EncodeUint(coeffsPol.Coeffs[0], plaintext)
	ciphertext = encryptor.EncryptNew(plaintext)
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(testCtx *testContext, decryptor bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {
	require.True(t, utils.EqualSliceUint64(coeffs, testCtx.encoder.DecodeUintNew(decryptor.DecryptNew(ciphertext))))
}

func testMarshalling(testCtx *testContext, t *testing.T) {

	//verify if the un.marshalling works properly

	crsGen := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)
	crs := crsGen.ReadNew()
	ringQ := testCtx.dbfvContext.ringQ
	ringP := testCtx.dbfvContext.ringP

	Ciphertext := bfv.NewCiphertextRandom(testCtx.prng, testCtx.params, 1)

	t.Run(fmt.Sprintf("Marshalling/CPK/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {
		keygenProtocol := NewCKGProtocol(testCtx.params)
		KeyGenShareBefore := keygenProtocol.AllocateShares()
		keygenProtocol.GenShare(&testCtx.sk0.SecretKey, crs, KeyGenShareBefore)
		//now we marshall it
		data, err := KeyGenShareBefore.MarshalBinary()

		if err != nil {
			log.Fatal("Could not marshal the CKGShare : ", err)

		}

		KeyGenShareAfter := new(drlwe.CKGShare)
		err = KeyGenShareAfter.UnmarshalBinary(data)
		if err != nil {
			log.Fatal("Could not unmarshal the CKGShare : ", err)

		}

		//comparing the results
		require.Equal(t, KeyGenShareBefore.Degree(), KeyGenShareAfter.Degree())
		require.Equal(t, KeyGenShareBefore.LenModuli(), KeyGenShareAfter.LenModuli())

		moduli := KeyGenShareBefore.LenModuli()
		require.Equal(t, KeyGenShareAfter.Coeffs[:moduli], KeyGenShareBefore.Coeffs[:moduli])
	})

	t.Run(fmt.Sprintf("Marshalling/PCKS/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {
		//Check marshalling for the PCKS

		KeySwitchProtocol := NewPCKSProtocol(testCtx.params, testCtx.params.Sigma())
		SwitchShare := KeySwitchProtocol.AllocateShares()
		KeySwitchProtocol.GenShare(testCtx.sk0.Value, testCtx.pk0, Ciphertext, SwitchShare)

		data, err := SwitchShare.MarshalBinary()
		require.NoError(t, err)

		SwitchShareReceiver := new(PCKSShare)
		err = SwitchShareReceiver.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := 0; i < 2; i++ {
			//compare the shares.
			ringBefore := SwitchShare[i]
			ringAfter := SwitchShareReceiver[i]
			require.Equal(t, ringBefore.Degree(), ringAfter.Degree())
			moduli := ringAfter.LenModuli()
			require.Equal(t, ringAfter.Coeffs[:moduli], ringBefore.Coeffs[:moduli])
		}
	})

	t.Run(fmt.Sprintf("Marshalling/CKS/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {

		//Now for CKSShare ~ its similar to PKSShare
		cksp := NewCKSProtocol(testCtx.params, testCtx.params.Sigma())
		cksshare := cksp.AllocateShare()
		cksp.GenShare(testCtx.sk0.Value, testCtx.sk1.Value, Ciphertext, cksshare)

		data, err := cksshare.MarshalBinary()
		require.NoError(t, err)
		cksshareAfter := new(CKSShare)
		err = cksshareAfter.UnmarshalBinary(data)
		require.NoError(t, err)

		//now compare both shares.

		require.Equal(t, cksshare.Degree(), cksshareAfter.Degree())
		require.Equal(t, cksshare.LenModuli(), cksshareAfter.LenModuli())

		moduli := cksshare.LenModuli()
		require.Equal(t, cksshare.Coeffs[:moduli], cksshareAfter.Coeffs[:moduli])
	})

	t.Run(fmt.Sprintf("Marshalling/Refresh/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {

		//testing refresh shares
		refreshproto := NewRefreshProtocol(testCtx.params)
		refreshshare := refreshproto.AllocateShares()
		refreshproto.GenShares(testCtx.sk0.Value, Ciphertext, crs, refreshshare)

		data, err := refreshshare.MarshalBinary()
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

	t.Run(fmt.Sprintf("Marshalling/RTG/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {

		//check RTGShare
		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)
		modulus := (testCtx.dbfvContext.ringQ.Modulus)
		crp := make([]*ring.Poly, len(modulus))
		for j := 0; j < len(modulus); j++ {
			crp[j] = crpGenerator.ReadNew()

		}

		galEl := testCtx.params.GaloisElementForColumnRotationBy(64)

		rotProto := NewRotKGProtocol(testCtx.params)
		rtgShare := rotProto.AllocateShares()
		rotProto.GenShare(&testCtx.sk1.SecretKey, galEl, crp, rtgShare)

		data, err := rtgShare.MarshalBinary()
		if err != nil {
			log.Fatal("could not marshal RTGshare :", err)
		}

		resRTGShare := new(drlwe.RTGShare)
		err = resRTGShare.UnmarshalBinary(data)
		if err != nil {
			log.Fatal("Could not unmarshal RTGShare: ", err)
		}

		if len(resRTGShare.Value) != len(rtgShare.Value) {
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

func testMarshallingRelin(testCtx *testContext, t *testing.T) {

	ringQ := testCtx.dbfvContext.ringQ
	ringP := testCtx.dbfvContext.ringP
	modulus := testCtx.dbfvContext.ringQ.Modulus

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dbfvContext.ringQP)

	crp := make([]*ring.Poly, len(modulus))
	for j := 0; j < len(modulus); j++ {
		crp[j] = crpGenerator.ReadNew()
	}

	t.Run(fmt.Sprintf("Marshalling/RLKG/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {

		rlk := NewRKGProtocol(testCtx.params)

		u, r1, r2 := rlk.AllocateShares()
		rlk.GenShareRoundOne(&testCtx.sk0.SecretKey, crp, u, r1)
		data, err := r1.MarshalBinary()
		require.NoError(t, err)

		r1After := new(drlwe.RKGShare)
		err = r1After.UnmarshalBinary(data)
		require.NoError(t, err)

		// for i := 0; i < (len(r1)); i++ { // TODO test in drlwe
		// 	a := r1[i][0]
		// 	b := (*r1After)[i][0]
		// 	moduli := a.LenModuli()
		// 	require.Equal(t, a.Coeffs[:moduli], b.Coeffs[:moduli])
		// }

		rlk.GenShareRoundTwo(u, &testCtx.sk0.SecretKey, r1, crp, r2)

		data, err = r2.MarshalBinary()
		require.NoError(t, err)

		r2After := new(drlwe.RKGShare)
		err = r2After.UnmarshalBinary(data)
		require.NoError(t, err)

		// for i := 0; i < (len(r2)); i++ { // TODO test in drlwe
		// 	for idx := 0; idx < 2; idx++ {
		// 		a := r2[i][idx]
		// 		b := (*r2After)[i][idx]
		// 		moduli := a.LenModuli()
		// 		require.Equal(t, a.Coeffs[:moduli], b.Coeffs[:moduli])
		// 	}

		// }
	})

}
