package mpckks

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
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

	crs            multiparty.CRS
	uniformSampler *ring.UniformSampler
}

func TestMULTIPARTYFloat(t *testing.T) {

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

	kgen := rlwe.NewKeyGenerator(tc.params)

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
	tc.encryptorPk0 = rlwe.NewEncryptor(tc.params, tc.pk0)
	tc.decryptorSk0 = rlwe.NewDecryptor(tc.params, tc.sk0)
	tc.decryptorSk1 = rlwe.NewDecryptor(tc.params, tc.sk1)

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
			publicShareE2S multiparty.KeySwitchShare
			publicShareS2E multiparty.KeySwitchShare
			secretShare    multiparty.AdditiveShareBigint
		}

		params := tc.params

		coeffs, _, ciphertext := newTestVectors(tc, tc.encryptorPk0, -1, 1, params.LogMaxSlots())

		tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

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

		ckks.VerifyTestVectors(params, tc.encoder, nil, coeffs, pt, params.LogDefaultScale(), 0, *printPrecisionStats, t)

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

		ckks.VerifyTestVectors(params, tc.encoder, tc.decryptorSk0, coeffs, ctRec, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testRefresh(tc *testContext, t *testing.T) {

	paramsIn := tc.params

	// To get the precision of the linear transformations
	_, logBound, _ := GetMinimumLevelForRefresh(128, paramsIn.DefaultScale(), tc.NParties, paramsIn.Q())

	t.Run(GetTestName("N->N/Transform=nil", tc.NParties, paramsIn), func(t *testing.T) {
		testRefreshParameterized(tc, paramsIn, tc.sk0Shards, nil, t)
	})

	t.Run(GetTestName("N->2N/Transform=nil", tc.NParties, paramsIn), func(t *testing.T) {

		var paramsOut ckks.Parameters
		var err error
		paramsOut, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:            paramsIn.LogN() + 1,
			LogQ:            []int{54, 54, 54, 49, 49, 49, 49, 49, 49},
			LogP:            []int{52, 52},
			RingType:        paramsIn.RingType(),
			LogDefaultScale: paramsIn.LogDefaultScale(),
		})

		require.NoError(t, err)

		kgenOut := rlwe.NewKeyGenerator(paramsOut)

		skOut := make([]*rlwe.SecretKey, tc.NParties)
		for i := range skOut {
			skOut[i] = kgenOut.GenSecretKeyNew()
		}

		testRefreshParameterized(tc, paramsOut, skOut, nil, t)
	})

	t.Run(GetTestName("2N->N/Transform=nil", tc.NParties, tc.params), func(t *testing.T) {

		var paramsOut ckks.Parameters
		var err error
		paramsOut, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:            paramsIn.LogN() - 1,
			LogQ:            []int{54, 54, 54, 49, 49, 49, 49, 49, 49},
			LogP:            []int{52, 52},
			RingType:        paramsIn.RingType(),
			LogDefaultScale: paramsIn.LogDefaultScale(),
		})

		require.NoError(t, err)

		kgenOut := rlwe.NewKeyGenerator(paramsOut)

		skOut := make([]*rlwe.SecretKey, tc.NParties)
		for i := range skOut {
			skOut[i] = kgenOut.GenSecretKeyNew()
		}

		testRefreshParameterized(tc, paramsOut, skOut, nil, t)
	})

	t.Run(GetTestName("N->N/Transform=true", tc.NParties, paramsIn), func(t *testing.T) {

		transform := &MaskedLinearTransformationFunc{
			Decode: true,
			Func: func(coeffs []*bignum.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], bignum.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], bignum.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		testRefreshParameterized(tc, paramsIn, tc.sk0Shards, transform, t)
	})

	t.Run(GetTestName("N->2N/Transform=true", tc.NParties, paramsIn), func(t *testing.T) {

		var paramsOut ckks.Parameters
		var err error
		paramsOut, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:            paramsIn.LogN() + 1,
			LogQ:            []int{54, 54, 54, 49, 49, 49, 49, 49, 49},
			LogP:            []int{52, 52},
			RingType:        paramsIn.RingType(),
			LogDefaultScale: paramsIn.LogDefaultScale(),
		})

		require.NoError(t, err)

		kgenOut := rlwe.NewKeyGenerator(paramsOut)

		skOut := make([]*rlwe.SecretKey, tc.NParties)
		for i := range skOut {
			skOut[i] = kgenOut.GenSecretKeyNew()
		}

		transform := &MaskedLinearTransformationFunc{
			Decode: true,
			Func: func(coeffs []*bignum.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], bignum.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], bignum.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		testRefreshParameterized(tc, paramsOut, skOut, transform, t)
	})

	t.Run(GetTestName("2N->N/Transform=true", tc.NParties, tc.params), func(t *testing.T) {

		var paramsOut ckks.Parameters
		var err error
		paramsOut, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
			LogN:            paramsIn.LogN() - 1,
			LogQ:            []int{54, 54, 54, 49, 49, 49, 49, 49, 49},
			LogP:            []int{52, 52},
			RingType:        paramsIn.RingType(),
			LogDefaultScale: paramsIn.LogDefaultScale(),
		})

		require.NoError(t, err)

		kgenOut := rlwe.NewKeyGenerator(paramsOut)

		skOut := make([]*rlwe.SecretKey, tc.NParties)
		for i := range skOut {
			skOut[i] = kgenOut.GenSecretKeyNew()
		}

		transform := &MaskedLinearTransformationFunc{
			Decode: true,
			Func: func(coeffs []*bignum.Complex) {
				for i := range coeffs {
					coeffs[i][0].Mul(coeffs[i][0], bignum.NewFloat(0.9238795325112867, logBound))
					coeffs[i][1].Mul(coeffs[i][1], bignum.NewFloat(0.7071067811865476, logBound))
				}
			},
			Encode: true,
		}

		testRefreshParameterized(tc, paramsOut, skOut, transform, t)
	})
}

func testRefreshParameterized(tc *testContext, paramsOut ckks.Parameters, skOut []*rlwe.SecretKey, transform *MaskedLinearTransformationFunc, t *testing.T) {

	var err error

	paramsIn := tc.params

	encIn := tc.encryptorPk0

	skIdealOut := rlwe.NewSecretKey(paramsOut)
	for i := 0; i < tc.NParties; i++ {
		paramsOut.RingQ().Add(skIdealOut.Value.Q, skOut[i].Value.Q, skIdealOut.Value.Q)
	}

	var minLevel int
	var logBound uint
	var ok bool
	if minLevel, logBound, ok = GetMinimumLevelForRefresh(128, paramsIn.DefaultScale(), tc.NParties, paramsIn.Q()); ok != true || minLevel+1 > paramsIn.MaxLevel() {
		t.Skip("Not enough levels to ensure correctness and 128 security")
	}

	type Party struct {
		MaskedLinearTransformationProtocol
		sIn   *rlwe.SecretKey
		sOut  *rlwe.SecretKey
		share multiparty.RefreshShare
	}

	coeffs, _, ciphertext := newTestVectors(tc, encIn, -1, 1, utils.Min(paramsIn.LogMaxSlots(), paramsOut.LogMaxSlots()))

	// Drops the ciphertext to the minimum level that ensures correctness and 128-bit security
	tc.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel-1)

	levelIn := minLevel

	require.Nil(t, err)

	levelOut := paramsOut.MaxLevel()

	RefreshParties := make([]*Party, tc.NParties)

	for i := 0; i < tc.NParties; i++ {
		p := new(Party)

		if i == 0 {
			if p.MaskedLinearTransformationProtocol, err = NewMaskedLinearTransformationProtocol(paramsIn, paramsOut, logBound, paramsIn.Xe()); err != nil {
				t.Log(err)
				t.Fail()
			}
		} else {
			p.MaskedLinearTransformationProtocol = RefreshParties[0].MaskedLinearTransformationProtocol.ShallowCopy()
		}

		p.sIn = tc.sk0Shards[i]
		p.sOut = skOut[i]

		p.share = p.AllocateShare(levelIn, levelOut)
		RefreshParties[i] = p
	}

	P0 := RefreshParties[0]
	crp := P0.SampleCRP(levelOut, tc.crs)

	for i, p := range RefreshParties {
		p.GenShare(p.sIn, p.sOut, logBound, ciphertext, crp, transform, &p.share)

		if i > 0 {
			P0.AggregateShares(&p.share, &P0.share, &P0.share)
		}
	}

	P0.Transform(ciphertext, transform, crp, P0.share, ciphertext)

	// Applies transform in plaintext

	if transform != nil {
		transform.Func(coeffs)
	}

	ckks.VerifyTestVectors(paramsOut, ckks.NewEncoder(paramsOut), rlwe.NewDecryptor(paramsOut, skIdealOut), coeffs, ciphertext, paramsOut.LogDefaultScale(), 0, *printPrecisionStats, t)
}

func newTestVectors(tc *testContext, encryptor *rlwe.Encryptor, a, b complex128, logSlots int) (values []*bignum.Complex, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {
	return newTestVectorsAtScale(tc, encryptor, a, b, tc.params.DefaultScale(), logSlots)
}

func newTestVectorsAtScale(tc *testContext, encryptor *rlwe.Encryptor, a, b complex128, scale rlwe.Scale, logSlots int) (values []*bignum.Complex, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	prec := tc.encoder.Prec()

	pt = ckks.NewPlaintext(tc.params, tc.params.MaxLevel())
	pt.Scale = scale
	pt.LogDimensions.Cols = logSlots

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

	if err := tc.encoder.Encode(values, pt); err != nil {
		panic(err)
	}

	if encryptor != nil {
		var err error
		if ct, err = encryptor.EncryptNew(pt); err != nil {
			panic(err)
		}
	}

	return values, pt, ct
}
