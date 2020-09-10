package dbfv

import (
	"fmt"
	"log"
	"math/big"
	"testing"

	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"github.com/stretchr/testify/require"
)

var err error
var params = new(testParams)
var parties uint64 = 3

func testString(opname string, parties uint64, params *bfv.Parameters) string {
	return fmt.Sprintf("%sparties=%d/LogN=%d/logQ=%d", opname, parties, params.LogN(), params.LogQP())
}

type testParams struct {
	params *bfv.Parameters

	dbfvContext *dbfvContext

	prng utils.PRNG

	encoder bfv.Encoder
	kgen    *bfv.KeyGenerator

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

	var defaultParams []*bfv.Parameters

	if testing.Short() {
		defaultParams = bfv.DefaultParams[bfv.PN12QP109 : bfv.PN12QP109+3]
	} else {
		defaultParams = bfv.DefaultParams
	}

	for _, p := range defaultParams {

		if err = genTestParams(p); err != nil {
			panic(err)
		}

		t.Run("PublicKeyGen", testPublicKeyGen)
		t.Run("RelinKeyGen", testRelinKeyGen)
		t.Run("RelinKeyGenNaive", testRelinKeyGenNaive)
		t.Run("KeySwitching", testKeyswitching)
		t.Run("PublicKeySwitching", testPublicKeySwitching)
		t.Run("RotKeyGenRotRows", testRotKeyGenRotRows)
		t.Run("RotKeyGenRotCols", testRotKeyGenRotCols)
		t.Run("Refresh", testRefresh)
		t.Run("RefreshAndPermute", testRefreshAndPermutation)
	}
}

func genTestParams(defaultParams *bfv.Parameters) (err error) {

	params.params = defaultParams.Copy()

	params.dbfvContext = newDbfvContext(params.params)

	params.prng, err = utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	params.encoder = bfv.NewEncoder(params.params)
	params.evaluator = bfv.NewEvaluator(params.params)

	kgen := bfv.NewKeyGenerator(params.params)

	// SecretKeys
	params.sk0Shards = make([]*bfv.SecretKey, parties)
	params.sk1Shards = make([]*bfv.SecretKey, parties)
	tmp0 := params.dbfvContext.ringQP.NewPoly()
	tmp1 := params.dbfvContext.ringQP.NewPoly()

	for j := uint64(0); j < parties; j++ {
		params.sk0Shards[j] = kgen.GenSecretKey()
		params.sk1Shards[j] = kgen.GenSecretKey()
		params.dbfvContext.ringQP.Add(tmp0, params.sk0Shards[j].Get(), tmp0)
		params.dbfvContext.ringQP.Add(tmp1, params.sk1Shards[j].Get(), tmp1)
	}

	params.sk0 = new(bfv.SecretKey)
	params.sk1 = new(bfv.SecretKey)

	params.sk0.Set(tmp0)
	params.sk1.Set(tmp1)

	// Publickeys
	params.pk0 = kgen.GenPublicKey(params.sk0)
	params.pk1 = kgen.GenPublicKey(params.sk1)

	params.encryptorPk0 = bfv.NewEncryptorFromPk(params.params, params.pk0)
	params.decryptorSk0 = bfv.NewDecryptor(params.params, params.sk0)
	params.decryptorSk1 = bfv.NewDecryptor(params.params, params.sk1)

	return
}

func testPublicKeyGen(t *testing.T) {

	sk0Shards := params.sk0Shards
	decryptorSk0 := params.decryptorSk0

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		crpGenerator := ring.NewUniformSampler(params.prng, params.dbfvContext.ringQP)
		crp := crpGenerator.ReadNew()

		type Party struct {
			*CKGProtocol
			s  *ring.Poly
			s1 CKGShare
		}

		ckgParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.CKGProtocol = NewCKGProtocol(params.params)
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
		encryptorTest := bfv.NewEncryptorFromPk(params.params, pk)

		coeffs, _, ciphertext := newTestVectors(encryptorTest, t)

		verifyTestVectors(decryptorSk0, coeffs, ciphertext, t)
	})
}

func testRelinKeyGen(t *testing.T) {

	sk0Shards := params.sk0Shards
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	evaluator := params.evaluator

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*RKGProtocol
			u      *ring.Poly
			s      *ring.Poly
			share1 RKGShare
			share2 RKGShare
		}

		rkgParties := make([]*Party, parties)

		for i := range rkgParties {
			p := new(Party)
			p.RKGProtocol = NewEkgProtocol(params.params)
			p.u = p.RKGProtocol.NewEphemeralKey()
			p.s = sk0Shards[i].Get()
			p.share1, p.share2 = p.RKGProtocol.AllocateShares()
			rkgParties[i] = p
		}

		P0 := rkgParties[0]

		crpGenerator := ring.NewUniformSampler(params.prng, params.dbfvContext.ringQP)
		crp := make([]*ring.Poly, params.params.Beta())

		for i := uint64(0); i < params.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
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
			p.GenShareRoundTwo(P0.share1, p.u, p.s, crp, p.share2)
			if i > 0 {
				P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
			}
		}

		evk := bfv.NewRelinKey(params.params, 1)
		P0.GenRelinearizationKey(P0.share1, P0.share2, evk)

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, t)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
			coeffs[i] %= params.dbfvContext.ringT.Modulus[0]
		}

		ciphertextMul := bfv.NewCiphertext(params.params, ciphertext.Degree()*2)
		evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

		res := bfv.NewCiphertext(params.params, 1)
		evaluator.Relinearize(ciphertextMul, evk, res)

		verifyTestVectors(decryptorSk0, coeffs, ciphertextMul, t)
	})

}

func testRelinKeyGenNaive(t *testing.T) {

	evaluator := params.evaluator
	pk0 := params.pk0
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards

	t.Run(testString("", parties, params.params), func(t *testing.T) {

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
			p.RKGProtocolNaive = NewRKGProtocolNaive(params.params)
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

		evk := bfv.NewRelinKey(params.params, 1)
		P0.GenRelinearizationKey(P0.share2, evk)

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, t)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
			coeffs[i] %= params.dbfvContext.ringT.Modulus[0]
		}

		ciphertextMul := bfv.NewCiphertext(params.params, ciphertext.Degree()*2)
		evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

		res := bfv.NewCiphertext(params.params, 1)
		evaluator.Relinearize(ciphertextMul, evk, res)

		verifyTestVectors(decryptorSk0, coeffs, ciphertextMul, t)
	})
}

func testKeyswitching(t *testing.T) {

	sk0Shards := params.sk0Shards
	sk1Shards := params.sk1Shards
	encryptorPk0 := params.encryptorPk0
	decryptorSk1 := params.decryptorSk1

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*CKSProtocol
			s0    *ring.Poly
			s1    *ring.Poly
			share CKSShare
		}

		cksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.CKSProtocol = NewCKSProtocol(params.params, 6.36)
			p.s0 = sk0Shards[i].Get()
			p.s1 = sk1Shards[i].Get()
			p.share = p.AllocateShare()
			cksParties[i] = p
		}
		P0 := cksParties[0]

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, t)

		// Each party creates its CKSProtocol instance with tmp = si-si'
		for i, p := range cksParties {
			p.GenShare(p.s0, p.s1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		ksCiphertext := bfv.NewCiphertext(params.params, 1)
		P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

		verifyTestVectors(decryptorSk1, coeffs, ksCiphertext, t)

		P0.KeySwitch(P0.share, ciphertext, ciphertext)

		verifyTestVectors(decryptorSk1, coeffs, ciphertext, t)

	})
}

func testPublicKeySwitching(t *testing.T) {

	sk0Shards := params.sk0Shards
	pk1 := params.pk1
	encryptorPk0 := params.encryptorPk0
	decryptorSk1 := params.decryptorSk1

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*PCKSProtocol
			s     *ring.Poly
			share PCKSShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.PCKSProtocol = NewPCKSProtocol(params.params, 6.36)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShares()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, t)

		ciphertextSwitched := bfv.NewCiphertext(params.params, 1)

		for i, p := range pcksParties {
			p.GenShare(p.s, pk1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

		verifyTestVectors(decryptorSk1, coeffs, ciphertextSwitched, t)
	})
}

func testRotKeyGenRotRows(t *testing.T) {

	evaluator := params.evaluator
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *ring.Poly
			share RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(params.params)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		crpGenerator := ring.NewUniformSampler(params.prng, params.dbfvContext.ringQP)
		crp := make([]*ring.Poly, params.params.Beta())

		for i := uint64(0); i < params.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		for i, p := range pcksParties {
			p.GenShare(bfv.RotationRow, 0, p.s, crp, &p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		rotkey := bfv.NewRotationKeys()
		P0.Finalize(P0.share, crp, rotkey)

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, t)

		evaluator.RotateRows(ciphertext, rotkey, ciphertext)

		coeffs = append(coeffs[params.params.N()>>1:], coeffs[:params.params.N()>>1]...)

		verifyTestVectors(decryptorSk0, coeffs, ciphertext, t)

	})
}

func testRotKeyGenRotCols(t *testing.T) {

	evaluator := params.evaluator
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *ring.Poly
			share RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(params.params)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}

		P0 := pcksParties[0]

		crpGenerator := ring.NewUniformSampler(params.prng, params.dbfvContext.ringQP)
		crp := make([]*ring.Poly, params.params.Beta())

		for i := uint64(0); i < params.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		mask := (params.params.N() >> 1) - 1

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, t)

		receiver := bfv.NewCiphertext(params.params, ciphertext.Degree())

		for k := uint64(1); k < params.params.N()>>1; k <<= 1 {

			for i, p := range pcksParties {
				p.GenShare(bfv.RotationLeft, k, p.s, crp, &p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			rotkey := bfv.NewRotationKeys()
			P0.Finalize(P0.share, crp, rotkey)

			evaluator.RotateColumns(ciphertext, k, rotkey, receiver)

			coeffsWant := make([]uint64, params.params.N())

			for i := uint64(0); i < params.params.N()>>1; i++ {
				coeffsWant[i] = coeffs[(i+k)&mask]
				coeffsWant[i+(params.params.N()>>1)] = coeffs[((i+k)&mask)+(params.params.N()>>1)]
			}

			verifyTestVectors(decryptorSk0, coeffsWant, receiver, t)
		}
	})
}

func testRefresh(t *testing.T) {

	encryptorPk0 := params.encryptorPk0
	sk0Shards := params.sk0Shards
	encoder := params.encoder
	decryptorSk0 := params.decryptorSk0

	kgen := bfv.NewKeyGenerator(params.params)

	rlk := kgen.GenRelinKey(params.sk0, 2)

	t.Run(fmt.Sprintf("N=%d/logQ=%d/Refresh", params.params.N(), params.dbfvContext.ringQP.ModulusBigint.BitLen()), func(t *testing.T) {

		type Party struct {
			*RefreshProtocol
			s       *ring.Poly
			share   RefreshShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RefreshProtocol = NewRefreshProtocol(params.params)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShares()
			p.ptShare = bfv.NewPlaintext(params.params)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(params.prng, params.dbfvContext.ringQP)
		crp := crpGenerator.ReadNew()

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, t)

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
				coeffsTmp[j] = ring.BRed(coeffsTmp[j], coeffsTmp[j], params.dbfvContext.ringT.Modulus[0], params.dbfvContext.ringT.GetBredParams()[0])
			}

			if utils.EqualSliceUint64(coeffsTmp, encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertextTmp))) {
				maxDepth++
			} else {
				break
			}
		}

		// Simulated added error of size Q/(T^2) and add it to the fresh ciphertext
		coeffsBigint := make([]*big.Int, params.params.N())
		params.dbfvContext.ringQ.PolyToBigint(ciphertext.Value()[0], coeffsBigint)

		errorRange := new(big.Int).Set(params.dbfvContext.ringQ.ModulusBigint)
		errorRange.Quo(errorRange, params.dbfvContext.ringT.ModulusBigint)
		errorRange.Quo(errorRange, params.dbfvContext.ringT.ModulusBigint)

		for i := uint64(0); i < params.params.N(); i++ {
			coeffsBigint[i].Add(coeffsBigint[i], ring.RandInt(errorRange))
		}

		params.dbfvContext.ringQ.SetCoefficientsBigint(coeffsBigint, ciphertext.Value()[0])

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
				coeffs[j] = ring.BRed(coeffs[j], coeffs[j], params.dbfvContext.ringT.Modulus[0], params.dbfvContext.ringT.GetBredParams()[0])
			}
		}

		//Decrypts and compare
		require.True(t, utils.EqualSliceUint64(coeffs, encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertext))))
	})
}

func testRefreshAndPermutation(t *testing.T) {

	encryptorPk0 := params.encryptorPk0
	sk0Shards := params.sk0Shards
	encoder := params.encoder
	decryptorSk0 := params.decryptorSk0

	t.Run(fmt.Sprintf("N=%d/logQ=%d/RefreshAndPermute", params.params.N(), params.dbfvContext.ringQP.ModulusBigint.BitLen()), func(t *testing.T) {

		type Party struct {
			*PermuteProtocol
			s       *ring.Poly
			share   RefreshShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.PermuteProtocol = NewPermuteProtocol(params.params)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShares()
			p.ptShare = bfv.NewPlaintext(params.params)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		crpGenerator := ring.NewUniformSampler(params.prng, params.dbfvContext.ringQP)
		crp := crpGenerator.ReadNew()

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, t)

		permutation := make([]uint64, len(coeffs))

		for i := range permutation {
			permutation[i] = ring.RandUniform(prng, params.params.N(), params.params.N()-1)
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

		coeffsHave := encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertext))

		//Decrypts and compare
		require.True(t, utils.EqualSliceUint64(coeffsPermute, coeffsHave))
	})
}

func newTestVectors(encryptor bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {

	uniformSampler := ring.NewUniformSampler(params.prng, params.dbfvContext.ringT)

	coeffsPol := uniformSampler.ReadNew()
	plaintext = bfv.NewPlaintext(params.params)
	params.encoder.EncodeUint(coeffsPol.Coeffs[0], plaintext)
	ciphertext = encryptor.EncryptNew(plaintext)
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(decryptor bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {
	require.True(t, utils.EqualSliceUint64(coeffs, params.encoder.DecodeUint(decryptor.DecryptNew(ciphertext))))
}

func Test_Marshalling(t *testing.T) {

	//verify if the un.marshalling works properly

	crsGen := ring.NewUniformSampler(params.prng, params.dbfvContext.ringQP)
	crs := crsGen.ReadNew()
	ringQ := params.dbfvContext.ringQ
	ringP := params.dbfvContext.ringP

	Ciphertext := bfv.NewCiphertextRandom(params.prng, params.params, 1)

	t.Run(fmt.Sprintf("CPK/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {
		keygenProtocol := NewCKGProtocol(params.params)
		KeyGenShareBefore := keygenProtocol.AllocateShares()
		keygenProtocol.GenShare(params.sk0.Get(), crs, KeyGenShareBefore)
		//now we marshall it
		data, err := KeyGenShareBefore.MarshalBinary()

		if err != nil {
			log.Fatal("Could not marshal the CKGShare : ", err)

		}

		KeyGenShareAfter := new(CKGShare)
		err = KeyGenShareAfter.UnmarshalBinary(data)
		if err != nil {
			log.Fatal("Could not unmarshal the CKGShare : ", err)

		}

		//comparing the results
		require.Equal(t, KeyGenShareBefore.GetDegree(), KeyGenShareAfter.GetDegree())
		require.Equal(t, KeyGenShareBefore.GetLenModuli(), KeyGenShareAfter.GetLenModuli())

		moduli := KeyGenShareBefore.GetLenModuli()
		require.Equal(t, KeyGenShareAfter.Coeffs[:moduli], KeyGenShareBefore.Coeffs[:moduli])
	})

	t.Run(fmt.Sprintf("PCKS/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {
		//Check marshalling for the PCKS

		KeySwitchProtocol := NewPCKSProtocol(params.params, params.dbfvContext.params.Sigma())
		SwitchShare := KeySwitchProtocol.AllocateShares()
		KeySwitchProtocol.GenShare(params.sk0.Get(), params.pk0, Ciphertext, SwitchShare)

		data, err := SwitchShare.MarshalBinary()
		require.NoError(t, err)

		SwitchShareReceiver := new(PCKSShare)
		err = SwitchShareReceiver.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := 0; i < 2; i++ {
			//compare the shares.
			ringBefore := SwitchShare[i]
			ringAfter := SwitchShareReceiver[i]
			require.Equal(t, ringBefore.GetDegree(), ringAfter.GetDegree())
			moduli := ringAfter.GetLenModuli()
			require.Equal(t, ringAfter.Coeffs[:moduli], ringBefore.Coeffs[:moduli])
		}
	})

	t.Run(fmt.Sprintf("CKS/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {

		//Now for CKSShare ~ its similar to PKSShare
		cksp := NewCKSProtocol(params.params, params.dbfvContext.params.Sigma())
		cksshare := cksp.AllocateShare()
		cksp.GenShare(params.sk0.Get(), params.sk1.Get(), Ciphertext, cksshare)

		data, err := cksshare.MarshalBinary()
		require.NoError(t, err)
		cksshareAfter := new(CKSShare)
		err = cksshareAfter.UnmarshalBinary(data)
		require.NoError(t, err)

		//now compare both shares.

		require.Equal(t, cksshare.GetDegree(), cksshareAfter.GetDegree())
		require.Equal(t, cksshare.GetLenModuli(), cksshareAfter.GetLenModuli())

		moduli := cksshare.GetLenModuli()
		require.Equal(t, cksshare.Coeffs[:moduli], cksshareAfter.Coeffs[:moduli])
	})

	t.Run(fmt.Sprintf("Refresh/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {

		//testing refresh shares
		refreshproto := NewRefreshProtocol(params.params)
		refreshshare := refreshproto.AllocateShares()
		refreshproto.GenShares(params.sk0.Get(), Ciphertext, crs, refreshshare)

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

	t.Run(fmt.Sprintf("RTG/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {

		//check RTGShare
		crpGenerator := ring.NewUniformSampler(params.prng, params.dbfvContext.ringQP)
		modulus := (params.dbfvContext.ringQ.Modulus)
		crp := make([]*ring.Poly, len(modulus))
		for j := 0; j < len(modulus); j++ {
			crp[j] = crpGenerator.ReadNew() //make([]*ring.Poly, bitLog)

		}

		rotProto := NewRotKGProtocol(params.params)
		rtgShare := rotProto.AllocateShare()
		rotProto.GenShare(1, 64, params.sk1.Get(), crp, &rtgShare)

		data, err := rtgShare.MarshalBinary()
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

	ringQ := params.dbfvContext.ringQ
	ringP := params.dbfvContext.ringP
	modulus := params.dbfvContext.ringQ.Modulus
	prng, err := utils.NewKeyedPRNG(nil)
	if err != nil {
		panic(err)
	}

	crpGenerator := ring.NewUniformSampler(prng, params.dbfvContext.ringQP)

	crp := make([]*ring.Poly, len(modulus))
	for j := 0; j < len(modulus); j++ {
		crp[j] = crpGenerator.ReadNew() //make([]*ring.Poly, bitLog)
		//for u := uint64(0); u < bitLog; u++ {
		//	crp[j][u] = crpGenerator.ClockUniformNew()
		//}
	}

	t.Run(fmt.Sprintf("RLKG/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(t *testing.T) {

		rlk := NewEkgProtocol(params.params)
		u := rlk.NewEphemeralKey()

		r1, r2 := rlk.AllocateShares()
		rlk.GenShareRoundOne(u, params.sk0.Get(), crp, r1)
		data, err := r1.MarshalBinary()
		require.NoError(t, err)

		r1After := new(RKGShare)
		err = r1After.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := 0; i < (len(r1)); i++ {
			a := r1[i][0]
			b := (*r1After)[i][0]
			moduli := a.GetLenModuli()
			require.Equal(t, a.Coeffs[:moduli], b.Coeffs[:moduli])
		}

		rlk.GenShareRoundTwo(r1, u, params.sk0.Get(), crp, r2)

		data, err = r2.MarshalBinary()
		require.NoError(t, err)

		r2After := new(RKGShare)
		err = r2After.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := 0; i < (len(r2)); i++ {
			for idx := 0; idx < 2; idx++ {
				a := r2[i][idx]
				b := (*r2After)[i][idx]
				moduli := a.GetLenModuli()
				require.Equal(t, a.Coeffs[:moduli], b.Coeffs[:moduli])
			}

		}
	})

}
