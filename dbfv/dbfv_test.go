package dbfv

import (
	"encoding/json"
	"flag"
	"fmt"
	"math/big"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v3/bfv"
	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var parties int = 3

func testString(opname string, parties int, params bfv.Parameters) string {
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/parties=%d", opname, params.LogN(), params.LogQP(), parties)
}

type testContext struct {
	params bfv.Parameters

	// Polynomial degree
	n int

	// Polynomial contexts
	ringT *ring.Ring
	ringQ *ring.Ring
	ringP *ring.Ring

	prng utils.PRNG

	encoder bfv.Encoder

	sk0Shards []*rlwe.SecretKey
	sk0       *rlwe.SecretKey

	sk1       *rlwe.SecretKey
	sk1Shards []*rlwe.SecretKey

	pk0 *rlwe.PublicKey
	pk1 *rlwe.PublicKey

	encryptorPk0 bfv.Encryptor
	decryptorSk0 bfv.Decryptor
	decryptorSk1 bfv.Decryptor
	evaluator    bfv.Evaluator

	crs            drlwe.CRS
	uniformSampler *ring.UniformSampler
}

func Test_DBFV(t *testing.T) {

	defaultParams := bfv.DefaultParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = bfv.DefaultParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = append(defaultParams, bfv.DefaultPostQuantumParams...) // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams bfv.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []bfv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {

		params, err := bfv.NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}
		var tc *testContext
		if tc, err = gentestContext(params); err != nil {
			panic(err)
		}
		for _, testSet := range []func(tc *testContext, t *testing.T){

			testPublicKeyGen,
			testRelinKeyGen,
			testKeyswitching,
			testPublicKeySwitching,
			testRotKeyGenRotRows,
			testRotKeyGenRotCols,
			testEncToShares,
			testRefresh,
			testRefreshAndPermutation,
			testMarshalling,
		} {
			testSet(tc, t)
			runtime.GC()
		}
	}
}

func gentestContext(params bfv.Parameters) (tc *testContext, err error) {

	tc = new(testContext)

	tc.params = params

	tc.n = params.N()

	tc.ringT = params.RingT()
	tc.ringQ = params.RingQ()
	tc.ringP = params.RingP()

	prng, _ := utils.NewKeyedPRNG([]byte{'t', 'e', 's', 't'})
	tc.crs = prng
	tc.uniformSampler = ring.NewUniformSampler(prng, params.RingQ())

	tc.encoder = bfv.NewEncoder(tc.params)
	tc.evaluator = bfv.NewEvaluator(tc.params, rlwe.EvaluationKey{})

	kgen := bfv.NewKeyGenerator(tc.params)

	// SecretKeys
	tc.sk0Shards = make([]*rlwe.SecretKey, parties)
	tc.sk1Shards = make([]*rlwe.SecretKey, parties)

	tc.sk0 = bfv.NewSecretKey(tc.params)
	tc.sk1 = bfv.NewSecretKey(tc.params)

	ringQP, levelQ, levelP := params.RingQP(), params.QCount()-1, params.PCount()-1
	for j := 0; j < parties; j++ {
		tc.sk0Shards[j] = kgen.GenSecretKey()
		tc.sk1Shards[j] = kgen.GenSecretKey()
		ringQP.AddLvl(levelQ, levelP, tc.sk0.Value, tc.sk0Shards[j].Value, tc.sk0.Value)
		ringQP.AddLvl(levelQ, levelP, tc.sk1.Value, tc.sk1Shards[j].Value, tc.sk1.Value)
	}

	// Publickeys
	tc.pk0 = kgen.GenPublicKey(tc.sk0)
	tc.pk1 = kgen.GenPublicKey(tc.sk1)

	tc.encryptorPk0 = bfv.NewEncryptor(tc.params, tc.pk0)
	tc.decryptorSk0 = bfv.NewDecryptor(tc.params, tc.sk0)
	tc.decryptorSk1 = bfv.NewDecryptor(tc.params, tc.sk1)

	return
}

func testPublicKeyGen(tc *testContext, t *testing.T) {

	sk0Shards := tc.sk0Shards
	decryptorSk0 := tc.decryptorSk0

	t.Run(testString("PublicKeyGen", parties, tc.params), func(t *testing.T) {

		type Party struct {
			*CKGProtocol
			s  *rlwe.SecretKey
			s1 *drlwe.CKGShare
		}

		ckgParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.CKGProtocol = NewCKGProtocol(tc.params)
			p.s = sk0Shards[i]
			p.s1 = p.AllocateShare()
			ckgParties[i] = p
		}
		P0 := ckgParties[0]

		crp := P0.SampleCRP(tc.crs)

		// Checks that dbfv.CKGProtocol complies to the drlwe.CollectivePublicKeyGenerator interface
		var _ drlwe.CollectivePublicKeyGenerator = P0.CKGProtocol

		// Each party creates a new CKGProtocol instance
		for i, p := range ckgParties {
			p.GenShare(p.s, crp, p.s1)
			if i > 0 {
				P0.AggregateShare(p.s1, P0.s1, P0.s1)
			}
		}

		pk := bfv.NewPublicKey(tc.params)
		P0.GenPublicKey(P0.s1, crp, pk)

		// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
		encryptorTest := bfv.NewEncryptor(tc.params, pk)

		coeffs, _, ciphertext := newTestVectors(tc, encryptorTest, t)

		verifyTestVectors(tc, decryptorSk0, coeffs, ciphertext, t)
	})
}

func testRelinKeyGen(tc *testContext, t *testing.T) {

	sk0Shards := tc.sk0Shards
	encryptorPk0 := tc.encryptorPk0
	decryptorSk0 := tc.decryptorSk0

	t.Run(testString("RelinKeyGen", parties, tc.params), func(t *testing.T) {

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
			p.RKGProtocol = NewRKGProtocol(tc.params)
			p.sk = sk0Shards[i]
			p.ephSk, p.share1, p.share2 = p.AllocateShare()
			rkgParties[i] = p
		}

		P0 := rkgParties[0]

		// Checks that bfv.RKGProtocol complies to the drlwe.RelinearizationKeyGenerator interface
		var _ drlwe.RelinearizationKeyGenerator = P0.RKGProtocol

		crp := P0.SampleCRP(tc.crs)

		// ROUND 1
		for i, p := range rkgParties {
			p.GenShareRoundOne(p.sk, crp, p.ephSk, p.share1)
			if i > 0 {
				P0.AggregateShare(p.share1, P0.share1, P0.share1)
			}
		}

		//ROUND 2
		for i, p := range rkgParties {
			p.GenShareRoundTwo(p.ephSk, p.sk, P0.share1, p.share2)
			if i > 0 {
				P0.AggregateShare(p.share2, P0.share2, P0.share2)
			}
		}

		evk := bfv.NewRelinearizationKey(tc.params, 1)
		P0.GenRelinearizationKey(P0.share1, P0.share2, evk)

		evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: evk, Rtks: nil})

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)
		for i := range coeffs {
			coeffs[i] *= coeffs[i]
			coeffs[i] %= tc.ringT.Modulus[0]
		}

		ciphertextMul := bfv.NewCiphertext(tc.params, ciphertext.Degree()*2)
		evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

		res := bfv.NewCiphertext(tc.params, 1)
		evaluator.Relinearize(ciphertextMul, res)

		verifyTestVectors(tc, decryptorSk0, coeffs, res, t)
	})

}

func testKeyswitching(tc *testContext, t *testing.T) {

	sk0Shards := tc.sk0Shards
	sk1Shards := tc.sk1Shards
	encryptorPk0 := tc.encryptorPk0
	decryptorSk1 := tc.decryptorSk1

	t.Run(testString("Keyswitching", parties, tc.params), func(t *testing.T) {

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		type Party struct {
			cks   *CKSProtocol
			s0    *rlwe.SecretKey
			s1    *rlwe.SecretKey
			share *drlwe.CKSShare
		}

		cksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.cks = NewCKSProtocol(tc.params, 6.36)
			p.s0 = sk0Shards[i]
			p.s1 = sk1Shards[i]
			p.share = p.cks.AllocateShare()
			cksParties[i] = p
		}
		P0 := cksParties[0]

		// Checks that the protocol complies to the drlwe.PublicKeySwitchingProtocol interface
		var _ drlwe.KeySwitchingProtocol = &P0.cks.CKSProtocol

		// Each party creates its CKSProtocol instance with tmp = si-si'
		for i, p := range cksParties {
			p.cks.GenShare(p.s0, p.s1, ciphertext.Value[1], p.share)
			if i > 0 {
				P0.cks.AggregateShare(p.share, P0.share, P0.share)
			}
		}

		ksCiphertext := bfv.NewCiphertext(tc.params, 1)
		P0.cks.KeySwitch(ciphertext, P0.share, ksCiphertext)

		verifyTestVectors(tc, decryptorSk1, coeffs, ksCiphertext, t)

		P0.cks.KeySwitch(ciphertext, P0.share, ciphertext)

		verifyTestVectors(tc, decryptorSk1, coeffs, ciphertext, t)

	})
}

func testPublicKeySwitching(tc *testContext, t *testing.T) {

	sk0Shards := tc.sk0Shards
	pk1 := tc.pk1
	encryptorPk0 := tc.encryptorPk0
	decryptorSk1 := tc.decryptorSk1

	t.Run(testString("PublicKeySwitching", parties, tc.params), func(t *testing.T) {

		type Party struct {
			*PCKSProtocol
			s     *rlwe.SecretKey
			share *drlwe.PCKSShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.PCKSProtocol = NewPCKSProtocol(tc.params, 6.36)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		// Checks that the protocol complies to the drlwe.PublicKeySwitchingProtocol interface
		var _ drlwe.PublicKeySwitchingProtocol = &P0.PCKSProtocol.PCKSProtocol

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		ciphertextSwitched := bfv.NewCiphertext(tc.params, 1)

		for i, p := range pcksParties {
			p.GenShare(p.s, pk1, ciphertext.Value[1], p.share)
			if i > 0 {
				P0.AggregateShare(p.share, P0.share, P0.share)
			}
		}

		P0.KeySwitch(ciphertext, P0.share, ciphertextSwitched)

		verifyTestVectors(tc, decryptorSk1, coeffs, ciphertextSwitched, t)
	})
}

func testRotKeyGenRotRows(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	decryptorSk0 := tc.decryptorSk0
	sk0Shards := tc.sk0Shards

	t.Run(testString("RotKeyGenRotRows", parties, tc.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *rlwe.SecretKey
			share *drlwe.RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(tc.params)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		// Checks that bfv.RTGProtocol complies to the drlwe.RotationKeyGenerator interface
		var _ drlwe.RotationKeyGenerator = P0.RTGProtocol

		crp := P0.SampleCRP(tc.crs)

		galEl := tc.params.GaloisElementForRowRotation()
		rotKeySet := bfv.NewRotationKeySet(tc.params, []uint64{galEl})

		for i, p := range pcksParties {
			p.GenShare(p.s, galEl, crp, p.share)
			if i > 0 {
				P0.AggregateShare(p.share, P0.share, P0.share)
			}
		}

		P0.GenRotationKey(P0.share, crp, rotKeySet.Keys[galEl])

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: nil, Rtks: rotKeySet})
		result := evaluator.RotateRowsNew(ciphertext)
		coeffsWant := append(coeffs[tc.params.N()>>1:], coeffs[:tc.params.N()>>1]...)

		verifyTestVectors(tc, decryptorSk0, coeffsWant, result, t)

	})
}

func testRotKeyGenRotCols(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	decryptorSk0 := tc.decryptorSk0
	sk0Shards := tc.sk0Shards

	t.Run(testString("RotKeyGenRotCols", parties, tc.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *rlwe.SecretKey
			share *drlwe.RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(tc.params)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}

		P0 := pcksParties[0]

		// Checks that bfv.RTGProtocol complies to the drlwe.RotationKeyGenerator interface
		var _ drlwe.RotationKeyGenerator = P0.RTGProtocol

		crp := P0.SampleCRP(tc.crs)

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		galEls := tc.params.GaloisElementsForRowInnerSum()
		rotKeySet := bfv.NewRotationKeySet(tc.params, galEls)

		for _, galEl := range galEls {

			for i, p := range pcksParties {
				p.GenShare(p.s, galEl, crp, p.share)
				if i > 0 {
					P0.AggregateShare(p.share, P0.share, P0.share)
				}
			}

			P0.GenRotationKey(P0.share, crp, rotKeySet.Keys[galEl])
		}

		evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: nil, Rtks: rotKeySet})
		for k := 1; k < tc.params.N()>>1; k <<= 1 {
			result := evaluator.RotateColumnsNew(ciphertext, int(k))
			coeffsWant := utils.RotateUint64Slots(coeffs, int(k))
			verifyTestVectors(tc, decryptorSk0, coeffsWant, result, t)
		}
	})
}

func testEncToShares(tc *testContext, t *testing.T) {

	coeffs, _, ciphertext := newTestVectors(tc, tc.encryptorPk0, t)

	type Party struct {
		e2s         *E2SProtocol
		s2e         *S2EProtocol
		sk          *rlwe.SecretKey
		publicShare *drlwe.CKSShare
		secretShare *rlwe.AdditiveShare
	}

	params := tc.params
	P := make([]Party, parties)

	for i := range P {
		if i == 0 {
			P[i].e2s = NewE2SProtocol(params, 3.2)
			P[i].s2e = NewS2EProtocol(params, 3.2)
		} else {
			P[i].e2s = P[0].e2s.ShallowCopy()
			P[i].s2e = P[0].s2e.ShallowCopy()
		}

		P[i].sk = tc.sk0Shards[i]
		P[i].publicShare = P[i].e2s.AllocateShare()
		P[i].secretShare = rlwe.NewAdditiveShare(params.Parameters)
	}

	// The E2S protocol is run in all tests, as a setup to the S2E test.
	for i, p := range P {

		p.e2s.GenShare(p.sk, ciphertext.Value[1], p.secretShare, p.publicShare)
		if i > 0 {
			p.e2s.AggregateShare(P[0].publicShare, p.publicShare, P[0].publicShare)
		}
	}

	P[0].e2s.GetShare(P[0].secretShare, P[0].publicShare, ciphertext, P[0].secretShare)

	t.Run(testString("E2SProtocol", parties, tc.params), func(t *testing.T) {

		rec := rlwe.NewAdditiveShare(params.Parameters)
		for _, p := range P {
			//fmt.Println("P[", i, "] share:", p.secretShare.Value.Coeffs[0][:see])
			tc.ringT.Add(&rec.Value, &p.secretShare.Value, &rec.Value)
		}

		ptRt := bfv.NewPlaintextRingT(tc.params)
		ptRt.Value.Copy(&rec.Value)

		assert.True(t, utils.EqualSliceUint64(coeffs, tc.encoder.DecodeUintNew(ptRt)))

	})

	crp := P[0].e2s.SampleCRP(params.MaxLevel(), tc.crs)

	t.Run(testString("S2EProtocol", parties, tc.params), func(t *testing.T) {

		for i, p := range P {
			p.s2e.GenShare(p.sk, crp, p.secretShare, p.publicShare)
			if i > 0 {
				p.s2e.AggregateShare(P[0].publicShare, p.publicShare, P[0].publicShare)
			}
		}

		ctRec := bfv.NewCiphertext(tc.params, 1)
		P[0].s2e.GetEncryption(P[0].publicShare, crp, ctRec)

		verifyTestVectors(tc, tc.decryptorSk0, coeffs, ctRec, t)
	})
}

func testRefresh(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	encoder := tc.encoder
	decryptorSk0 := tc.decryptorSk0

	kgen := bfv.NewKeyGenerator(tc.params)

	rlk := kgen.GenRelinearizationKey(tc.sk0, 1)

	t.Run(testString("Refresh", parties, tc.params), func(t *testing.T) {

		type Party struct {
			*RefreshProtocol
			s       *rlwe.SecretKey
			share   *RefreshShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			if i == 0 {
				p.RefreshProtocol = NewRefreshProtocol(tc.params, 3.2)
			} else {
				p.RefreshProtocol = RefreshParties[0].RefreshProtocol.ShallowCopy()
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			p.ptShare = bfv.NewPlaintext(tc.params)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crp := P0.SampleCRP(tc.params.MaxLevel(), tc.crs)

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		maxDepth := 0

		ciphertextTmp := ciphertext.CopyNew()
		coeffsTmp := make([]uint64, len(coeffs))

		copy(coeffsTmp, coeffs)

		evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: rlk, Rtks: nil})
		// Finds the maximum multiplicative depth
		for {

			evaluator.Relinearize(tc.evaluator.MulNew(ciphertextTmp, ciphertextTmp), ciphertextTmp)

			for j := range coeffsTmp {
				coeffsTmp[j] = ring.BRed(coeffsTmp[j], coeffsTmp[j], tc.ringT.Modulus[0], tc.ringT.BredParams[0])
			}

			if utils.EqualSliceUint64(coeffsTmp, encoder.DecodeUintNew(decryptorSk0.DecryptNew(ciphertextTmp))) {
				maxDepth++
			} else {
				break
			}
		}

		// Simulated added error of size Q/(T^2) and add it to the fresh ciphertext
		coeffsBigint := make([]*big.Int, tc.params.N())
		tc.ringQ.PolyToBigint(ciphertext.Value[0], 1, coeffsBigint)

		errorRange := new(big.Int).Set(tc.ringQ.ModulusAtLevel[tc.params.MaxLevel()])
		errorRange.Quo(errorRange, tc.ringT.ModulusAtLevel[0])
		errorRange.Quo(errorRange, tc.ringT.ModulusAtLevel[0])

		for i := 0; i < tc.params.N(); i++ {
			coeffsBigint[i].Add(coeffsBigint[i], ring.RandInt(errorRange))
		}

		tc.ringQ.SetCoefficientsBigint(coeffsBigint, ciphertext.Value[0])

		for i, p := range RefreshParties {
			p.GenShare(p.s, ciphertext.Value[1], crp, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}

		}

		ctRes := bfv.NewCiphertext(tc.params, 1)
		P0.Finalize(ciphertext, crp, P0.share, ctRes)

		// Squares the refreshed ciphertext up to the maximum depth-1
		for i := 0; i < maxDepth-1; i++ {

			evaluator.Relinearize(tc.evaluator.MulNew(ctRes, ctRes), ctRes)

			for j := range coeffs {
				coeffs[j] = ring.BRed(coeffs[j], coeffs[j], tc.ringT.Modulus[0], tc.ringT.BredParams[0])
			}
		}

		//Decrypts and compare
		require.True(t, utils.EqualSliceUint64(coeffs, encoder.DecodeUintNew(decryptorSk0.DecryptNew(ctRes))))
	})
}

func testRefreshAndPermutation(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	encoder := tc.encoder
	decryptorSk0 := tc.decryptorSk0

	t.Run(testString("RefreshAndPermutation", parties, tc.params), func(t *testing.T) {

		type Party struct {
			*MaskedTransformProtocol
			s       *rlwe.SecretKey
			share   *MaskedTransformShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			if i == 0 {
				p.MaskedTransformProtocol = NewMaskedTransformProtocol(tc.params, 3.2)
			} else {
				p.MaskedTransformProtocol = NewMaskedTransformProtocol(tc.params, 3.2)
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			p.ptShare = bfv.NewPlaintext(tc.params)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crp := P0.SampleCRP(tc.params.MaxLevel(), tc.crs)

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		permutation := make([]uint64, len(coeffs))
		N := uint64(tc.params.N())
		prng, _ := utils.NewPRNG()
		for i := range permutation {
			permutation[i] = ring.RandUniform(prng, N, N-1)
		}

		permute := func(coeffs []uint64) {
			coeffsPerm := make([]uint64, len(coeffs))
			for i := range coeffs {
				coeffsPerm[i] = coeffs[permutation[i]]
			}
			copy(coeffs, coeffsPerm)
		}

		for i, p := range RefreshParties {
			p.GenShare(p.s, ciphertext.Value[1], crp, permute, p.share)
			if i > 0 {
				P0.Aggregate(P0.share, p.share, P0.share)
			}
		}

		P0.Transform(ciphertext, permute, crp, P0.share, ciphertext)

		coeffsPermute := make([]uint64, len(coeffs))
		for i := range coeffsPermute {
			coeffsPermute[i] = coeffs[permutation[i]]
		}

		coeffsHave := encoder.DecodeUintNew(decryptorSk0.DecryptNew(ciphertext))

		//Decrypts and compares
		require.True(t, utils.EqualSliceUint64(coeffsPermute, coeffsHave))
	})
}

func newTestVectors(tc *testContext, encryptor bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {

	prng, _ := utils.NewPRNG()
	uniformSampler := ring.NewUniformSampler(prng, tc.ringT)
	coeffsPol := uniformSampler.ReadNew()
	plaintext = bfv.NewPlaintext(tc.params)
	tc.encoder.Encode(coeffsPol.Coeffs[0], plaintext)
	ciphertext = encryptor.EncryptNew(plaintext)
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(tc *testContext, decryptor bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {
	require.True(t, utils.EqualSliceUint64(coeffs, tc.encoder.DecodeUintNew(decryptor.DecryptNew(ciphertext))))
}

func testMarshalling(tc *testContext, t *testing.T) {
	ciphertext := bfv.NewCiphertext(tc.params, 1)
	tc.uniformSampler.Read(ciphertext.Value[0])
	tc.uniformSampler.Read(ciphertext.Value[1])

	t.Run(testString("MarshallingRefresh", parties, tc.params), func(t *testing.T) {

		// Testing refresh shares
		refreshproto := NewRefreshProtocol(tc.params, 3.2)
		refreshshare := refreshproto.AllocateShare()

		crp := refreshproto.SampleCRP(tc.params.MaxLevel(), tc.crs)

		refreshproto.GenShare(tc.sk0, ciphertext.Value[1], crp, refreshshare)

		data, err := refreshshare.MarshalBinary()
		if err != nil {
			t.Fatal("Could not marshal RefreshShare", err)
		}
		resRefreshShare := new(MaskedTransformShare)
		err = resRefreshShare.UnmarshalBinary(data)

		if err != nil {
			t.Fatal("Could not unmarshal RefreshShare", err)
		}
		for i, r := range refreshshare.e2sShare.Value.Coeffs {
			if !utils.EqualSliceUint64(resRefreshShare.e2sShare.Value.Coeffs[i], r) {
				t.Fatal("Result of marshalling not the same as original : RefreshShare")
			}

		}
		for i, r := range refreshshare.s2eShare.Value.Coeffs {
			if !utils.EqualSliceUint64(resRefreshShare.s2eShare.Value.Coeffs[i], r) {
				t.Fatal("Result of marshalling not the same as original : RefreshShare")
			}
		}
	})
}
