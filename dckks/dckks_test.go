package dckks

import (
	"encoding/json"
	"flag"
	"fmt"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure refresh). Overrides -short and requires -timeout=0.")
var flagPostQuantum = flag.Bool("pq", false, "run post quantum test suite (does not run non-PQ parameters).")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")
var minPrec float64 = 15.0

func testString(opname string, parties int, params ckks.Parameters) string {
	return fmt.Sprintf("%s/RingType=%s/logN=%d/logSlots=%d/logQ=%d/levels=%d/#Pi=%d/Decomp=%d/parties=%d",
		opname,
		params.RingType(),
		params.LogN(),
		params.LogSlots(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.PCount(),
		params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP()),
		parties)
}

type testContext struct {
	params   ckks.Parameters
	NParties int

	ringQ *ring.Ring
	ringP *ring.Ring

	encoder   ckks.Encoder
	evaluator ckks.Evaluator

	encryptorPk0 rlwe.Encryptor
	decryptorSk0 rlwe.Decryptor
	decryptorSk1 rlwe.Decryptor

	pk0 *rlwe.PublicKey
	pk1 *rlwe.PublicKey

	sk0 *rlwe.SecretKey
	sk1 *rlwe.SecretKey

	sk0Shards []*rlwe.SecretKey
	sk1Shards []*rlwe.SecretKey

	crs            drlwe.CRS
	uniformSampler *ring.UniformSampler
}

func TestDCKKS(t *testing.T) {

	var err error

	var testParams []ckks.ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, ckks.ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			t.Fatal(err)
		}
	case *flagLongTest:
		for _, pls := range [][]ckks.ParametersLiteral{
			ckks.DefaultParams,
			ckks.DefaultConjugateInvariantParams,
			ckks.DefaultPostQuantumParams,
			ckks.DefaultPostQuantumConjugateInvariantParams} {
			testParams = append(testParams, pls...)
		}
	case *flagPostQuantum && testing.Short():
		testParams = append(ckks.DefaultPostQuantumParams[:2], ckks.DefaultPostQuantumConjugateInvariantParams[:2]...)
	case *flagPostQuantum:
		testParams = append(ckks.DefaultPostQuantumParams[:4], ckks.DefaultPostQuantumConjugateInvariantParams[:4]...)
	case testing.Short():
		testParams = append(ckks.DefaultParams[:2], ckks.DefaultConjugateInvariantParams[:2]...)
	default:
		testParams = append(ckks.DefaultParams[:4], ckks.DefaultConjugateInvariantParams[:4]...)
	}

	for _, paramsLiteral := range testParams[:] {

		var params ckks.Parameters
		if params, err = ckks.NewParametersFromLiteral(paramsLiteral); err != nil {
			t.Fatal(err)
		}
		N := 3
		var tc *testContext
		if tc, err = genTestParams(params, N); err != nil {
			t.Fatal(err)
		}

		for _, testSet := range []func(tc *testContext, t *testing.T){
			testE2SProtocol,
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

func genTestParams(params ckks.Parameters, NParties int) (tc *testContext, err error) {

	tc = new(testContext)

	tc.params = params

	tc.NParties = NParties

	tc.ringQ = params.RingQ()
	tc.ringP = params.RingP()

	prng, _ := sampling.NewKeyedPRNG([]byte{'t', 'e', 's', 't'})
	tc.crs = prng
	tc.uniformSampler = ring.NewUniformSampler(prng, params.RingQ())

	tc.encoder = ckks.NewEncoder(tc.params)
	tc.evaluator = ckks.NewEvaluator(tc.params, nil)

	kgen := ckks.NewKeyGenerator(tc.params)

	// SecretKeys
	tc.sk0Shards = make([]*rlwe.SecretKey, NParties)
	tc.sk1Shards = make([]*rlwe.SecretKey, NParties)
	tc.sk0 = rlwe.NewSecretKey(tc.params.Parameters)
	tc.sk1 = rlwe.NewSecretKey(tc.params.Parameters)

	ringQP := params.RingQP()
	for j := 0; j < NParties; j++ {
		tc.sk0Shards[j] = kgen.GenSecretKeyNew()
		tc.sk1Shards[j] = kgen.GenSecretKeyNew()
		ringQP.Add(tc.sk0.Value, tc.sk0Shards[j].Value, tc.sk0.Value)
		ringQP.Add(tc.sk1.Value, tc.sk1Shards[j].Value, tc.sk1.Value)
	}

	// Publickeys
	tc.pk0 = kgen.GenPublicKeyNew(tc.sk0)
	tc.pk1 = kgen.GenPublicKeyNew(tc.sk1)

	tc.encryptorPk0 = ckks.NewEncryptor(tc.params, tc.pk0)
	tc.decryptorSk0 = ckks.NewDecryptor(tc.params, tc.sk0)
	tc.decryptorSk1 = ckks.NewDecryptor(tc.params, tc.sk1)

	return
}

func testE2SProtocol(tc *testContext, t *testing.T) {

	params := tc.params

	t.Run(testString("E2SProtocol", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		type Party struct {
			e2s            *E2SProtocol
			s2e            *S2EProtocol
			sk             *rlwe.SecretKey
			publicShareE2S *drlwe.CKSShare
			publicShareS2E *drlwe.CKSShare
			secretShare    *rlwe.AdditiveShareBigint
		}

		coeffs, _, ciphertext := newTestVectors(tc, tc.encryptorPk0, -1, 1)

		tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

		params := tc.params
		P := make([]Party, tc.NParties)
		for i := range P {
			P[i].e2s = NewE2SProtocol(params, 3.2)
			P[i].s2e = NewS2EProtocol(params, 3.2)
			P[i].sk = tc.sk0Shards[i]
			P[i].publicShareE2S = P[i].e2s.AllocateShare(minLevel)
			P[i].publicShareS2E = P[i].s2e.AllocateShare(params.MaxLevel())
			P[i].secretShare = NewAdditiveShareBigint(params, params.LogSlots())
		}

		for i, p := range P {
			// Enc(-M_i)
			p.e2s.GenShare(p.sk, logBound, params.LogSlots(), ciphertext, p.secretShare, p.publicShareE2S)

			if i > 0 {
				// Enc(sum(-M_i))
				p.e2s.AggregateShares(P[0].publicShareE2S, p.publicShareE2S, P[0].publicShareE2S)
			}
		}

		// sum(-M_i) + x
		P[0].e2s.GetShare(P[0].secretShare, P[0].publicShareE2S, params.LogSlots(), ciphertext, P[0].secretShare)

		// sum(-M_i) + x + sum(M_i) = x
		rec := NewAdditiveShareBigint(params, params.LogSlots())
		for _, p := range P {
			a := rec.Value
			b := p.secretShare.Value

			for i := range a {
				a[i].Add(a[i], b[i])
			}
		}

		pt := ckks.NewPlaintext(params, ciphertext.Level())
		pt.IsNTT = false
		pt.Scale = ciphertext.Scale
		tc.ringQ.AtLevel(pt.Level()).SetCoefficientsBigint(rec.Value, pt.Value)

		verifyTestVectors(tc, nil, coeffs, pt, t)

		crp := P[0].s2e.SampleCRP(params.MaxLevel(), tc.crs)

		for i, p := range P {
			p.s2e.GenShare(p.sk, crp, params.LogSlots(), p.secretShare, p.publicShareS2E)
			if i > 0 {
				p.s2e.AggregateShares(P[0].publicShareS2E, p.publicShareS2E, P[0].publicShareS2E)
			}
		}

		ctRec := ckks.NewCiphertext(params, 1, params.MaxLevel())
		ctRec.Scale = params.DefaultScale()
		P[0].s2e.GetEncryption(P[0].publicShareS2E, crp, ctRec)

		verifyTestVectors(tc, tc.decryptorSk0, coeffs, ctRec, t)

	})
}

func testRefresh(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	decryptorSk0 := tc.decryptorSk0
	params := tc.params

	t.Run(testString("Refresh", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		type Party struct {
			*RefreshProtocol
			s     *rlwe.SecretKey
			share *RefreshShare
		}

		levelIn := minLevel
		levelOut := params.MaxLevel()

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			if i == 0 {
				p.RefreshProtocol = NewRefreshProtocol(params, logBound, 3.2)
			} else {
				p.RefreshProtocol = RefreshParties[0].RefreshProtocol.ShallowCopy()
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare(levelIn, levelOut)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		for _, scale := range []float64{params.DefaultScale().Float64(), params.DefaultScale().Float64() * 128} {
			t.Run(fmt.Sprintf("atScale=%f", scale), func(t *testing.T) {
				coeffs, _, ciphertext := newTestVectorsAtScale(tc, encryptorPk0, -1, 1, rlwe.NewScale(scale))

				// Brings ciphertext to minLevel + 1
				tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

				crp := P0.SampleCRP(levelOut, tc.crs)

				for i, p := range RefreshParties {

					p.GenShare(p.s, logBound, params.LogSlots(), ciphertext, crp, p.share)

					if i > 0 {
						P0.AggregateShares(p.share, P0.share, P0.share)
					}
				}

				P0.Finalize(ciphertext, params.LogSlots(), crp, P0.share, ciphertext)

				verifyTestVectors(tc, decryptorSk0, coeffs, ciphertext, t)
			})
		}

	})
}

func testRefreshAndTransform(tc *testContext, t *testing.T) {

	var err error
	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	params := tc.params
	decryptorSk0 := tc.decryptorSk0

	t.Run(testString("RefreshAndTransform", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		type Party struct {
			*MaskedTransformProtocol
			s     *rlwe.SecretKey
			share *MaskedTransformShare
		}

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, -1, 1)

		// Drops the ciphertext to the minimum level that ensures correctness and 128-bit security
		tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

		levelIn := minLevel
		levelOut := params.MaxLevel()

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)

			if i == 0 {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(params, params, logBound, 3.2); err != nil {
					t.Log(err)
					t.Fail()
				}
			} else {
				p.MaskedTransformProtocol = RefreshParties[0].MaskedTransformProtocol.ShallowCopy()
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare(levelIn, levelOut)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]
		crp := P0.SampleCRP(levelOut, tc.crs)

		transform := &MaskedTransformFunc{
			Decode: true,
			Func: func(coeffs []*ring.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], ring.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], ring.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		for i, p := range RefreshParties {
			p.GenShare(p.s, p.s, logBound, params.LogSlots(), ciphertext, crp, transform, p.share)

			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.Transform(ciphertext, tc.params.LogSlots(), transform, crp, P0.share, ciphertext)

		for i := range coeffs {
			coeffs[i] = complex(real(coeffs[i])*0.9238795325112867, imag(coeffs[i])*0.7071067811865476)
		}

		verifyTestVectors(tc, decryptorSk0, coeffs, ciphertext, t)
	})
}

func testRefreshAndTransformSwitchParams(tc *testContext, t *testing.T) {

	var err error

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	params := tc.params

	t.Run(testString("RefreshAndTransformAndSwitchParams", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		type Party struct {
			*MaskedTransformProtocol
			sIn   *rlwe.SecretKey
			sOut  *rlwe.SecretKey
			share *MaskedTransformShare
		}

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, -1, 1)

		// Drops the ciphertext to the minimum level that ensures correctness and 128-bit security
		tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

		levelIn := minLevel

		// Target parameters
		var paramsOut ckks.Parameters
		paramsOut, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:     params.LogN() + 1,
			LogQ:     []int{54, 49, 49, 49, 49, 49, 49},
			LogP:     []int{52, 52},
			RingType: params.RingType(),
			LogSlots: params.MaxLogSlots() + 1,
			LogScale: 49,
		})

		require.Nil(t, err)

		levelOut := paramsOut.MaxLevel()

		RefreshParties := make([]*Party, tc.NParties)

		kgenParamsOut := rlwe.NewKeyGenerator(paramsOut.Parameters)
		skIdealOut := rlwe.NewSecretKey(paramsOut.Parameters)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)

			if i == 0 {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(params, paramsOut, logBound, 3.2); err != nil {
					t.Log(err)
					t.Fail()
				}
			} else {
				p.MaskedTransformProtocol = RefreshParties[0].MaskedTransformProtocol.ShallowCopy()
			}

			p.sIn = sk0Shards[i]

			p.sOut = kgenParamsOut.GenSecretKeyNew() // New shared secret key in target parameters
			paramsOut.RingQ().Add(skIdealOut.Value.Q, p.sOut.Value.Q, skIdealOut.Value.Q)

			p.share = p.AllocateShare(levelIn, levelOut)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]
		crp := P0.SampleCRP(levelOut, tc.crs)

		transform := &MaskedTransformFunc{
			Decode: true,
			Func: func(coeffs []*ring.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], ring.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], ring.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		for i, p := range RefreshParties {
			p.GenShare(p.sIn, p.sOut, logBound, params.LogSlots(), ciphertext, crp, transform, p.share)

			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.Transform(ciphertext, tc.params.LogSlots(), transform, crp, P0.share, ciphertext)

		for i := range coeffs {
			coeffs[i] = complex(real(coeffs[i])*0.9238795325112867, imag(coeffs[i])*0.7071067811865476)
		}

		precStats := ckks.GetPrecisionStats(paramsOut, ckks.NewEncoder(paramsOut), nil, coeffs, ckks.NewDecryptor(paramsOut, skIdealOut).DecryptNew(ciphertext), params.LogSlots(), 0)

		if *printPrecisionStats {
			t.Log(precStats.String())
		}

		require.GreaterOrEqual(t, precStats.MeanPrecision.Real, minPrec)
		require.GreaterOrEqual(t, precStats.MeanPrecision.Imag, minPrec)
	})
}

func testMarshalling(tc *testContext, t *testing.T) {
	params := tc.params

	t.Run(testString("Marshalling/Refresh", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		ciphertext := ckks.NewCiphertext(params, 1, minLevel)
		ciphertext.Scale = params.DefaultScale()
		tc.uniformSampler.AtLevel(minLevel).Read(ciphertext.Value[0])
		tc.uniformSampler.AtLevel(minLevel).Read(ciphertext.Value[1])

		// Testing refresh shares
		refreshproto := NewRefreshProtocol(tc.params, logBound, 3.2)
		refreshshare := refreshproto.AllocateShare(ciphertext.Level(), params.MaxLevel())

		crp := refreshproto.SampleCRP(params.MaxLevel(), tc.crs)

		refreshproto.GenShare(tc.sk0, logBound, params.LogSlots(), ciphertext, crp, refreshshare)

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

func newTestVectors(testContext *testContext, encryptor rlwe.Encryptor, a, b complex128) (values []complex128, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {
	return newTestVectorsAtScale(testContext, encryptor, a, b, testContext.params.DefaultScale())
}

func newTestVectorsAtScale(testContext *testContext, encryptor rlwe.Encryptor, a, b complex128, scale rlwe.Scale) (values []complex128, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {

	params := testContext.params

	logSlots := params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	for i := 0; i < 1<<logSlots; i++ {
		values[i] = complex(sampling.RandFloat64(real(a), real(b)), sampling.RandFloat64(imag(a), imag(b)))
	}

	plaintext = testContext.encoder.EncodeNew(values, params.MaxLevel(), scale, params.LogSlots())

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(tc *testContext, decryptor rlwe.Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	precStats := ckks.GetPrecisionStats(tc.params, tc.encoder, decryptor, valuesWant, element, tc.params.LogSlots(), 0)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, precStats.MeanPrecision.Real, minPrec)
	require.GreaterOrEqual(t, precStats.MeanPrecision.Imag, minPrec)
}
