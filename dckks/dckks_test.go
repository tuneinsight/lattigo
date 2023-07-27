package dckks

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func GetTestName(opname string, parties int, params ckks.Parameters) string {
	return fmt.Sprintf("%s/RingType=%s/logN=%d/logQP=%d/Qi=%d/Pi=%d/LogDefaultScale=%d/Parties=%d",
		opname,
		params.RingType(),
		params.LogN(),
		int(math.Round(params.LogQP())),
		params.QCount(),
		params.PCount(),
		int(math.Log2(params.DefaultScale().Float64())),
		parties)
}

type testContext struct {
	params   ckks.Parameters
	NParties int

	ringQ *ring.Ring
	ringP *ring.Ring

	encoder   *ckks.Encoder
	evaluator *ckks.Evaluator

	encryptorPk0 *rlwe.Encryptor
	decryptorSk0 *rlwe.Decryptor
	decryptorSk1 *rlwe.Decryptor

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
	default:
		testParams = testParamsLiteral
	}

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		for _, paramsLiteral := range testParams {

			paramsLiteral.RingType = ringType

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
				testEncToShareProtocol,
				testRefresh,
				testRefreshAndTransform,
				testRefreshAndTransformSwitchParams,
			} {
				testSet(tc, t)
				runtime.GC()
			}
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
	if tc.pk0, err = kgen.GenPublicKeyNew(tc.sk0); err != nil {
		return
	}

	if tc.pk1, err = kgen.GenPublicKeyNew(tc.sk1); err != nil {
		return
	}

	if tc.encryptorPk0, err = ckks.NewEncryptor(tc.params, tc.pk0); err != nil {
		return
	}

	if tc.decryptorSk0, err = ckks.NewDecryptor(tc.params, tc.sk0); err != nil {
		return
	}

	if tc.decryptorSk1, err = ckks.NewDecryptor(tc.params, tc.sk1); err != nil {
		return
	}

	return
}

func testEncToShareProtocol(tc *testContext, t *testing.T) {

	params := tc.params

	t.Run(GetTestName("EncToShareProtocol", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		type Party struct {
			e2s            EncToShareProtocol
			s2e            ShareToEncProtocol
			sk             *rlwe.SecretKey
			publicShareE2S drlwe.KeySwitchShare
			publicShareS2E drlwe.KeySwitchShare
			secretShare    drlwe.AdditiveShareBigint
		}

		coeffs, _, ciphertext := newTestVectors(tc, tc.encryptorPk0, -1, 1)

		tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

		params := tc.params
		P := make([]Party, tc.NParties)
		var err error
		for i := range P {

			P[i].e2s, err = NewEncToShareProtocol(params, params.Xe())
			require.NoError(t, err)

			P[i].s2e, err = NewShareToEncProtocol(params, params.Xe())
			require.NoError(t, err)

			P[i].sk = tc.sk0Shards[i]
			P[i].publicShareE2S = P[i].e2s.AllocateShare(minLevel)
			P[i].publicShareS2E = P[i].s2e.AllocateShare(params.MaxLevel())
			P[i].secretShare = NewAdditiveShare(params, ciphertext.LogSlots())
		}

		for i, p := range P {
			// Enc(-M_i)
			p.e2s.GenShare(p.sk, logBound, ciphertext, &p.secretShare, &p.publicShareE2S)

			if i > 0 {
				// Enc(sum(-M_i))
				p.e2s.AggregateShares(P[0].publicShareE2S, p.publicShareE2S, &P[0].publicShareE2S)
			}
		}

		// sum(-M_i) + x
		P[0].e2s.GetShare(&P[0].secretShare, P[0].publicShareE2S, ciphertext, &P[0].secretShare)

		// sum(-M_i) + x + sum(M_i) = x
		rec := NewAdditiveShare(params, ciphertext.LogSlots())
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
			p.s2e.GenShare(p.sk, crp, ciphertext.MetaData, p.secretShare, &p.publicShareS2E)
			if i > 0 {
				p.s2e.AggregateShares(P[0].publicShareS2E, p.publicShareS2E, &P[0].publicShareS2E)
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

	t.Run(GetTestName("Refresh", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		type Party struct {
			RefreshProtocol
			s     *rlwe.SecretKey
			share drlwe.RefreshShare
		}

		levelIn := minLevel
		levelOut := params.MaxLevel()

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			var err error
			if i == 0 {
				p.RefreshProtocol, err = NewRefreshProtocol(params, logBound, params.Xe())
				require.NoError(t, err)
			} else {
				p.RefreshProtocol = RefreshParties[0].RefreshProtocol.ShallowCopy()
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare(levelIn, levelOut)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		for _, scale := range []float64{params.DefaultScale().Float64(), params.DefaultScale().Float64() * 128} {
			t.Run(fmt.Sprintf("AtScale=%d", int(math.Round(math.Log2(scale)))), func(t *testing.T) {
				coeffs, _, ciphertext := newTestVectorsAtScale(tc, encryptorPk0, -1, 1, rlwe.NewScale(scale))

				// Brings ciphertext to minLevel + 1
				tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

				crp := P0.SampleCRP(levelOut, tc.crs)

				for i, p := range RefreshParties {

					p.GenShare(p.s, logBound, ciphertext, crp, &p.share)

					if i > 0 {
						P0.AggregateShares(&p.share, &P0.share, &P0.share)
					}
				}

				P0.Finalize(ciphertext, crp, P0.share, ciphertext)

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

	t.Run(GetTestName("RefreshAndTransform", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		type Party struct {
			MaskedTransformProtocol
			s     *rlwe.SecretKey
			share drlwe.RefreshShare
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
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(params, params, logBound, params.Xe()); err != nil {
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
			Func: func(coeffs []*bignum.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], bignum.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], bignum.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		for i, p := range RefreshParties {
			p.GenShare(p.s, p.s, logBound, ciphertext, crp, transform, &p.share)

			if i > 0 {
				P0.AggregateShares(&p.share, &P0.share, &P0.share)
			}
		}

		P0.Transform(ciphertext, transform, crp, P0.share, ciphertext)

		for i := range coeffs {
			coeffs[i][0].Mul(coeffs[i][0], bignum.NewFloat(0.9238795325112867, logBound))
			coeffs[i][1].Mul(coeffs[i][1], bignum.NewFloat(0.7071067811865476, logBound))
		}

		verifyTestVectors(tc, decryptorSk0, coeffs, ciphertext, t)
	})
}

func testRefreshAndTransformSwitchParams(tc *testContext, t *testing.T) {

	var err error

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	params := tc.params

	t.Run(GetTestName("RefreshAndTransformAndSwitchParams", tc.NParties, params), func(t *testing.T) {

		var minLevel int
		var logBound uint
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, params.DefaultScale(), tc.NParties, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
			t.Skip("Not enough levels to ensure correctness and 128 security")
		}

		type Party struct {
			MaskedTransformProtocol
			sIn   *rlwe.SecretKey
			sOut  *rlwe.SecretKey
			share drlwe.RefreshShare
		}

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, -1, 1)

		// Drops the ciphertext to the minimum level that ensures correctness and 128-bit security
		tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

		levelIn := minLevel

		// Target parameters
		var paramsOut ckks.Parameters
		paramsOut, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:            params.LogN() + 1,
			LogQ:            []int{54, 49, 49, 49, 49, 49, 49},
			LogP:            []int{52, 52},
			RingType:        params.RingType(),
			LogDefaultScale: 49,
		})

		require.Nil(t, err)

		levelOut := paramsOut.MaxLevel()

		RefreshParties := make([]*Party, tc.NParties)

		kgenParamsOut := rlwe.NewKeyGenerator(paramsOut.Parameters)
		skIdealOut := rlwe.NewSecretKey(paramsOut.Parameters)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)

			if i == 0 {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(params, paramsOut, logBound, params.Xe()); err != nil {
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
			Func: func(coeffs []*bignum.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], bignum.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], bignum.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		for i, p := range RefreshParties {
			p.GenShare(p.sIn, p.sOut, logBound, ciphertext, crp, transform, &p.share)

			if i > 0 {
				P0.AggregateShares(&p.share, &P0.share, &P0.share)
			}
		}

		P0.Transform(ciphertext, transform, crp, P0.share, ciphertext)

		for i := range coeffs {
			coeffs[i][0].Mul(coeffs[i][0], bignum.NewFloat(0.9238795325112867, logBound))
			coeffs[i][1].Mul(coeffs[i][1], bignum.NewFloat(0.7071067811865476, logBound))
		}

		dec, err := ckks.NewDecryptor(paramsOut, skIdealOut)
		require.NoError(t, err)

		precStats := ckks.GetPrecisionStats(paramsOut, ckks.NewEncoder(paramsOut), nil, coeffs, dec.DecryptNew(ciphertext), nil, false)

		if *printPrecisionStats {
			t.Log(precStats.String())
		}

		rf64, _ := precStats.MeanPrecision.Real.Float64()
		if64, _ := precStats.MeanPrecision.Imag.Float64()

		minPrec := math.Log2(paramsOut.DefaultScale().Float64())
		switch params.RingType() {
		case ring.Standard:
			minPrec -= float64(paramsOut.LogN()) + 2
		case ring.ConjugateInvariant:
			minPrec -= float64(paramsOut.LogN()) + 2.5
		}
		if minPrec < 0 {
			minPrec = 0
		}

		require.GreaterOrEqual(t, rf64, minPrec)
		require.GreaterOrEqual(t, if64, minPrec)
	})
}

func newTestVectors(tc *testContext, encryptor *rlwe.Encryptor, a, b complex128) (values []*bignum.Complex, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {
	return newTestVectorsAtScale(tc, encryptor, a, b, tc.params.DefaultScale())
}

func newTestVectorsAtScale(tc *testContext, encryptor *rlwe.Encryptor, a, b complex128, scale rlwe.Scale) (values []*bignum.Complex, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	prec := tc.encoder.Prec()

	pt = ckks.NewPlaintext(tc.params, tc.params.MaxLevel())
	pt.Scale = scale

	values = make([]*bignum.Complex, pt.Slots())

	switch tc.params.RingType() {
	case ring.Standard:
		for i := range values {
			values[i] = &bignum.Complex{
				bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
				bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
			}
		}
	case ring.ConjugateInvariant:
		for i := range values {
			values[i] = &bignum.Complex{
				bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
				new(big.Float),
			}
		}
	default:
		panic("invalid ring type")
	}

	tc.encoder.Encode(values, pt)

	if encryptor != nil {
		var err error
		ct, err = encryptor.EncryptNew(pt)
		if err != nil {
			panic(err)
		}
	}

	return values, pt, ct
}

func verifyTestVectors(tc *testContext, decryptor *rlwe.Decryptor, valuesWant, valuesHave interface{}, t *testing.T) {

	precStats := ckks.GetPrecisionStats(tc.params, tc.encoder, decryptor, valuesWant, valuesHave, nil, false)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64, _ := precStats.MeanPrecision.Real.Float64()
	if64, _ := precStats.MeanPrecision.Imag.Float64()

	minPrec := math.Log2(tc.params.DefaultScale().Float64()) - float64(tc.params.LogN()+2)
	if minPrec < 0 {
		minPrec = 0
	}

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}
