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
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
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

	encryptorPk0 rlwe.Encryptor
	decryptorSk0 rlwe.Decryptor
	decryptorSk1 rlwe.Decryptor
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
			testEncToShares,
			testRefresh,
			testRefreshAndTransform,
			testRefreshAndTransformSwitchParams,
			testMarshalling,
		} {
			testSet(tc, t)
			runtime.GC()
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

	tc.sk0 = rlwe.NewSecretKey(tc.params.Parameters)
	tc.sk1 = rlwe.NewSecretKey(tc.params.Parameters)

	ringQP, levelQ, levelP := params.RingQP(), params.QCount()-1, params.PCount()-1
	for j := 0; j < parties; j++ {
		tc.sk0Shards[j] = kgen.GenSecretKey()
		tc.sk1Shards[j] = kgen.GenSecretKey()
		ringQP.AddLvl(levelQ, levelP, tc.sk0.Poly, tc.sk0Shards[j].Poly, tc.sk0.Poly)
		ringQP.AddLvl(levelQ, levelP, tc.sk1.Poly, tc.sk1Shards[j].Poly, tc.sk1.Poly)
	}

	// Publickeys
	tc.pk0 = kgen.GenPublicKey(tc.sk0)
	tc.pk1 = kgen.GenPublicKey(tc.sk1)

	tc.encryptorPk0 = bfv.NewEncryptor(tc.params, tc.pk0)
	tc.decryptorSk0 = bfv.NewDecryptor(tc.params, tc.sk0)
	tc.decryptorSk1 = bfv.NewDecryptor(tc.params, tc.sk1)

	return
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
		P[i].publicShare = P[i].e2s.AllocateShare(ciphertext.Level())
		P[i].secretShare = rlwe.NewAdditiveShare(params.Parameters)
	}

	// The E2S protocol is run in all tests, as a setup to the S2E test.
	for i, p := range P {

		p.e2s.GenShare(p.sk, ciphertext.Value[1], p.secretShare, p.publicShare)
		if i > 0 {
			p.e2s.AggregateShares(P[0].publicShare, p.publicShare, P[0].publicShare)
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
				p.s2e.AggregateShares(P[0].publicShare, p.publicShare, P[0].publicShare)
			}
		}

		ctRec := bfv.NewCiphertext(tc.params, 1, tc.params.MaxLevel())
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
			s     *rlwe.SecretKey
			share *RefreshShare
		}

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			if i == 0 {
				p.RefreshProtocol = NewRefreshProtocol(tc.params, 3.2)
			} else {
				p.RefreshProtocol = RefreshParties[0].RefreshProtocol.ShallowCopy()
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare(ciphertext.Level(), tc.params.MaxLevel())
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crp := P0.SampleCRP(tc.params.MaxLevel(), tc.crs)

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
				P0.AggregateShares(p.share, P0.share, P0.share)
			}

		}

		ctRes := bfv.NewCiphertext(tc.params, 1, tc.params.MaxLevel())
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

func testRefreshAndTransform(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	encoder := tc.encoder
	decryptorSk0 := tc.decryptorSk0

	t.Run(testString("RefreshAndPermutation", tc.NParties, tc.params), func(t *testing.T) {

		var err error

		type Party struct {
			*MaskedTransformProtocol
			s     *rlwe.SecretKey
			share *MaskedTransformShare
		}

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			if i == 0 {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(tc.params, tc.params, 3.2); err != nil {
					t.Fatal(err)
				}
			} else {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(tc.params, tc.params, 3.2); err != nil {
					t.Fatal(err)
				}
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare(ciphertext.Level(), tc.params.MaxLevel())
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crp := P0.SampleCRP(tc.params.MaxLevel(), tc.crs)

		permutation := make([]uint64, len(coeffs))
		N := uint64(tc.params.N())
		prng, _ := utils.NewPRNG()
		for i := range permutation {
			permutation[i] = ring.RandUniform(prng, N, N-1)
		}

		transform := &MaskedTransformFunc{
			Decode: true,
			Func: func(coeffs []uint64) {
				coeffsPerm := make([]uint64, len(coeffs))
				for i := range coeffs {
					coeffsPerm[i] = coeffs[permutation[i]]
				}
				copy(coeffs, coeffsPerm)
			},
			Encode: true,
		}

		for i, p := range RefreshParties {
			p.GenShare(p.s, p.s, ciphertext.Value[1], crp, transform, p.share)
			if i > 0 {
				P0.AggregateShares(P0.share, p.share, P0.share)
			}
		}

		P0.Transform(ciphertext, transform, crp, P0.share, ciphertext)

		coeffsPermute := make([]uint64, len(coeffs))
		for i := range coeffsPermute {
			coeffsPermute[i] = coeffs[permutation[i]]
		}

		coeffsHave := encoder.DecodeUintNew(decryptorSk0.DecryptNew(ciphertext))

		//Decrypts and compares
		require.True(t, utils.EqualSliceUint64(coeffsPermute, coeffsHave))
	})
}

func testRefreshAndTransformSwitchParams(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	paramsIn := tc.params

	t.Run(testString("RefreshAndTransformSwitchparams", tc.NParties, tc.params), func(t *testing.T) {

		// Checks that T is also a valid modulus for the next ring degree
		if paramsIn.T()&uint64(4*paramsIn.N()-1) != 1 {
			t.Skip("modulus T is not congruent to 1 mod 4N")
		}

		var paramsOut bfv.Parameters
		var err error
		paramsOut, err = bfv.NewParametersFromLiteral(bfv.ParametersLiteral{
			LogN: paramsIn.LogN(),
			LogQ: []int{54, 49, 49, 49},
			LogP: []int{52, 52},
			T:    paramsIn.T(),
		})

		require.Nil(t, err)

		type Party struct {
			*MaskedTransformProtocol
			sIn   *rlwe.SecretKey
			sOut  *rlwe.SecretKey
			share *MaskedTransformShare
		}

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		RefreshParties := make([]*Party, tc.NParties)
		kgenParamsOut := rlwe.NewKeyGenerator(paramsOut.Parameters)
		skIdealOut := rlwe.NewSecretKey(paramsOut.Parameters)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			if i == 0 {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(paramsIn, paramsOut, 3.2); err != nil {
					t.Fatal(err)
				}
			} else {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(paramsIn, paramsOut, 3.2); err != nil {
					t.Fatal(err)
				}
			}

			p.sIn = sk0Shards[i]

			p.sOut = kgenParamsOut.GenSecretKey() // New shared secret key in target parameters
			paramsOut.RingQ().Add(skIdealOut.Q, p.sOut.Q, skIdealOut.Q)

			p.share = p.AllocateShare(ciphertext.Level(), paramsOut.MaxLevel())

			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crp := P0.SampleCRP(paramsOut.MaxLevel(), tc.crs)

		permutation := make([]uint64, len(coeffs))
		N := uint64(tc.params.N())
		prng, _ := utils.NewPRNG()
		for i := range permutation {
			permutation[i] = ring.RandUniform(prng, N, N-1)
		}

		transform := &MaskedTransformFunc{
			Decode: true,
			Func: func(coeffs []uint64) {
				coeffsPerm := make([]uint64, len(coeffs))
				for i := range coeffs {
					coeffsPerm[i] = coeffs[permutation[i]]
				}
				copy(coeffs, coeffsPerm)
			},
			Encode: true,
		}

		for i, p := range RefreshParties {
			p.GenShare(p.sIn, p.sOut, ciphertext.Value[1], crp, transform, p.share)
			if i > 0 {
				P0.AggregateShares(P0.share, p.share, P0.share)
			}
		}

		P0.Transform(ciphertext, transform, crp, P0.share, ciphertext)

		transform.Func(coeffs)

		coeffsHave := bfv.NewEncoder(paramsOut).DecodeUintNew(bfv.NewDecryptor(paramsOut, skIdealOut).DecryptNew(ciphertext))

		//Decrypts and compares
		require.True(t, utils.EqualSliceUint64(coeffs, coeffsHave))
	})
}

func newTestVectors(tc *testContext, encryptor rlwe.Encryptor, t *testing.T) (coeffs []uint64, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	prng, _ := utils.NewPRNG()
	uniformSampler := ring.NewUniformSampler(prng, tc.ringT)
	coeffsPol := uniformSampler.ReadNew()
	pt = bfv.NewPlaintext(tc.params, tc.params.MaxLevel())
	tc.encoder.Encode(coeffsPol.Coeffs[0], pt)
	return coeffsPol.Coeffs[0], pt, encryptor.EncryptNew(pt)
}

func verifyTestVectors(tc *testContext, decryptor rlwe.Decryptor, coeffs []uint64, ct *rlwe.Ciphertext, t *testing.T) {
	require.True(t, utils.EqualSliceUint64(coeffs, tc.encoder.DecodeUintNew(decryptor.DecryptNew(ct))))
}

func testMarshalling(tc *testContext, t *testing.T) {
	ciphertext := bfv.NewCiphertext(tc.params, 1, tc.params.MaxLevel())
	tc.uniformSampler.Read(ciphertext.Value[0])
	tc.uniformSampler.Read(ciphertext.Value[1])

	t.Run(testString("MarshallingRefresh", tc.NParties, tc.params), func(t *testing.T) {

		// Testing refresh shares
		refreshproto := NewRefreshProtocol(tc.params, 3.2)
		refreshshare := refreshproto.AllocateShare(ciphertext.Level(), tc.params.MaxLevel())

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
