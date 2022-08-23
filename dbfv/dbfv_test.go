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

func testString(opname string, parties int, params bfv.Parameters) string {
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/parties=%d", opname, params.LogN(), params.LogQP(), parties)
}

type testContext struct {
	params bfv.Parameters

	NParties int

	// Polynomial degree
	n int

	// Polynomial contexts
	ringT *ring.Ring
	ringQ *ring.Ring
	ringP *ring.Ring

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

func TestDBFV(t *testing.T) {

	var err error

	defaultParams := bfv.DefaultParams[:] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = bfv.DefaultParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = append(defaultParams, bfv.DefaultPostQuantumParams...) // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams bfv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		defaultParams = []bfv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {

		var params bfv.Parameters
		if params, err = bfv.NewParametersFromLiteral(p); err != nil {
			t.Fatal(err)
		}

		var tc *testContext
		N := 3
		if tc, err = gentestContext(params, N); err != nil {
			t.Fatal(err)
		}
		for _, testSet := range []func(tc *testContext, t *testing.T){

			testKeyswitching,
			testPublicKeySwitching,
			testEncToShares,
			testRefresh,
			testRefreshAndPermutation,
			testMarshalling,
		} {
			testSet(tc, t)
			runtime.GC()
		}

		for _, N := range []int{2, 3, 5, 7, 20} {
			if tc, err = gentestContext(params, N); err != nil {
				panic(err)
			}
			testThreshold(tc, t)
		}

	}
}

func gentestContext(params bfv.Parameters, parties int) (tc *testContext, err error) {

	tc = new(testContext)

	tc.params = params

	tc.NParties = parties

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

func testKeyswitching(tc *testContext, t *testing.T) {

	sk0Shards := tc.sk0Shards
	sk1Shards := tc.sk1Shards
	encryptorPk0 := tc.encryptorPk0
	decryptorSk1 := tc.decryptorSk1

	t.Run(testString("Keyswitching", tc.NParties, tc.params), func(t *testing.T) {

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		type Party struct {
			cks   *CKSProtocol
			s0    *rlwe.SecretKey
			s1    *rlwe.SecretKey
			share *drlwe.CKSShare
		}

		cksParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			p.cks = NewCKSProtocol(tc.params, 6.36)
			p.s0 = sk0Shards[i]
			p.s1 = sk1Shards[i]
			p.share = p.cks.AllocateShare()
			cksParties[i] = p
		}
		P0 := cksParties[0]

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

	t.Run(testString("PublicKeySwitching", tc.NParties, tc.params), func(t *testing.T) {

		type Party struct {
			*PCKSProtocol
			s     *rlwe.SecretKey
			share *drlwe.PCKSShare
		}

		pcksParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			p.PCKSProtocol = NewPCKSProtocol(tc.params, 6.36)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

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
	P := make([]Party, tc.NParties)

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

	t.Run(testString("E2SProtocol", tc.NParties, tc.params), func(t *testing.T) {

		rec := rlwe.NewAdditiveShare(params.Parameters)
		for _, p := range P {
			tc.ringT.Add(&rec.Value, &p.secretShare.Value, &rec.Value)
		}

		ptRt := bfv.NewPlaintextRingT(tc.params)
		ptRt.Value.Copy(&rec.Value)

		assert.True(t, utils.EqualSliceUint64(coeffs, tc.encoder.DecodeUintNew(ptRt)))
	})

	crp := P[0].e2s.SampleCRP(params.MaxLevel(), tc.crs)

	t.Run(testString("S2EProtocol", tc.NParties, tc.params), func(t *testing.T) {
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

	t.Run(testString("Refresh", tc.NParties, tc.params), func(t *testing.T) {

		type Party struct {
			*RefreshProtocol
			s       *rlwe.SecretKey
			share   *RefreshShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
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
				P0.AggregateShare(p.share, P0.share, P0.share)
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

	t.Run(testString("RefreshAndPermutation", tc.NParties, tc.params), func(t *testing.T) {

		type Party struct {
			*MaskedTransformProtocol
			s       *rlwe.SecretKey
			share   *MaskedTransformShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
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
				P0.AggregateShare(P0.share, p.share, P0.share)
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

	t.Run(testString("MarshallingRefresh", tc.NParties, tc.params), func(t *testing.T) {

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

//Tests if the shares generated by the thresholdizer yield a correct key switch.
func testThreshold(tc *testContext, t *testing.T) {
	sk0Shards := tc.sk0Shards

	threshold := tc.NParties / 2

	t.Run(testString("Threshold", tc.NParties, tc.params)+fmt.Sprintf("/threshold=%d", threshold), func(t *testing.T) {

		type Party struct {
			*drlwe.Thresholdizer
			drlwe.Combiner
			*drlwe.CachedCombiner
			gen       *drlwe.ShamirPolynomial
			sk        *rlwe.SecretKey
			tsks      *drlwe.ShamirSecretShare
			tsk       *rlwe.SecretKey
			tpk       drlwe.ShamirPublicKey
			pcksShare *drlwe.PCKSShare
		}

		pcksPhase := func(params bfv.Parameters, tpk *rlwe.PublicKey, ct *bfv.Ciphertext, P []*Party) (encOut *bfv.Ciphertext) {

			// Collective key switching from the collective secret key to
			// the target public key

			pcks := NewPCKSProtocol(params, 3.19)

			for _, pi := range P {
				pi.pcksShare = pcks.AllocateShare()
			}

			for _, pi := range P {
				pcks.GenShare(pi.tsk, tpk, ct.Ciphertext.Value[1], pi.pcksShare)
			}

			pcksCombined := pcks.AllocateShare()
			encOut = bfv.NewCiphertext(params, 1)
			for _, pi := range P {
				pcks.AggregateShare(pi.pcksShare, pcksCombined, pcksCombined)
			}
			pcks.KeySwitch(ct, pcksCombined, encOut)

			return

		}

		P := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			p.Thresholdizer = drlwe.NewThresholdizer(tc.params.Parameters)
			p.Combiner = drlwe.NewCombiner(tc.params.Parameters, threshold)
			p.CachedCombiner = drlwe.NewCachedCombiner(tc.params.Parameters, threshold)
			p.sk = sk0Shards[i]
			p.tsk = bfv.NewSecretKey(tc.params)
			p.tpk = drlwe.ShamirPublicKey(i + 1)
			p.tsks = p.Thresholdizer.AllocateThresholdSecretShare()
			P[i] = p
		}

		shares := make(map[*Party]map[*Party]*drlwe.ShamirSecretShare, tc.NParties)
		var err error
		// Every party generates a share for every other party
		for _, pi := range P {

			pi.gen, err = pi.Thresholdizer.GenShamirPolynomial(threshold, pi.sk)
			if err != nil {
				t.Error(err)
			}

			shares[pi] = make(map[*Party]*drlwe.ShamirSecretShare)
			for _, pj := range P {
				shares[pi][pj] = pi.Thresholdizer.AllocateThresholdSecretShare()
				pi.Thresholdizer.GenShamirSecretShare(pj.tpk, pi.gen, shares[pi][pj])
			}
		}

		//Each party aggregates what it has received into a secret key
		for _, pi := range P {
			for _, pj := range P {
				pi.Thresholdizer.AggregateShares(pi.tsks, shares[pj][pi], pi.tsks)
			}
		}

		// Determining which parties are active. In a distributed context, a party
		// would receive the ids of active players and retrieve (or compute) the corresponding keys.
		activeParties := P[:threshold]
		activeShamirPks := make([]drlwe.ShamirPublicKey, threshold)
		for i, p := range activeParties {
			activeShamirPks[i] = p.tpk
		}

		// Combining
		// Slow because each party has to generate its public key on-the-fly. In
		// practice the public key could be precomputed from an id by parties during setup
		for _, pi := range activeParties {
			pi.Combiner.GenAdditiveShare(activeShamirPks, pi.tpk, pi.tsks, pi.tsk)
			tsk := pi.tsk.Value.CopyNew()
			pi.CachedCombiner.Precompute(activeShamirPks, pi.tpk)
			pi.CachedCombiner.GenAdditiveShare(activeShamirPks, pi.tpk, pi.tsks, pi.tsk)
			//the cached and non-cached combiners should yield the same results
			require.True(t, tc.ringQ.Equal(tsk.Q, pi.tsk.Value.Q))
			if tc.ringP != nil {
				require.True(t, tc.ringP.Equal(tsk.P, pi.tsk.Value.P))
			}
		}

		//Clearing caches
		for _, pi := range activeParties {
			pi.CachedCombiner.ClearCache()
		}

		coeffs, _, ciphertext := newTestVectors(tc, tc.encryptorPk0, t)

		ciphertextSwitched := pcksPhase(tc.params, tc.pk1, ciphertext, activeParties)

		verifyTestVectors(tc, tc.decryptorSk1, coeffs, ciphertextSwitched, t)

	})
}
