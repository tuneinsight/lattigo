package mpbgv

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"runtime"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func GetTestName(opname string, p bgv.Parameters, parties int) string {
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/logP=%d/LogSlots=%dx%d/logT=%d/Qi=%d/Pi=%d/parties=%d",
		opname,
		p.LogN(),
		int(math.Round(p.LogQ())),
		int(math.Round(p.LogP())),
		p.LogMaxDimensions().Rows,
		p.LogMaxDimensions().Cols,
		int(math.Round(p.LogT())),
		p.QCount(),
		p.PCount(),
		parties)
}

type testContext struct {
	params bgv.Parameters

	// Number of parties
	NParties int

	// Polynomial degree
	n int

	// Polynomial contexts
	ringT *ring.Ring
	ringQ *ring.Ring
	ringP *ring.Ring

	encoder *bgv.Encoder

	sk0Shards []*rlwe.SecretKey
	sk0       *rlwe.SecretKey

	sk1       *rlwe.SecretKey
	sk1Shards []*rlwe.SecretKey

	pk0 *rlwe.PublicKey
	pk1 *rlwe.PublicKey

	encryptorPk0 *rlwe.Encryptor
	decryptorSk0 *rlwe.Decryptor
	decryptorSk1 *rlwe.Decryptor
	evaluator    *bgv.Evaluator

	crs            multiparty.CRS
	uniformSampler *ring.UniformSampler
}

func TestInteger(t *testing.T) {

	var err error

	paramsLiterals := testParams

	if *flagParamString != "" {
		var jsonParams bgv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []bgv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals {

		for _, plaintextModulus := range testPlaintextModulus[:] {

			p.PlaintextModulus = plaintextModulus

			var params bgv.Parameters
			if params, err = bgv.NewParametersFromLiteral(p); err != nil {
				t.Fatal(err)
			}

			nParties := 3

			var tc *testContext
			if tc, err = gentestContext(nParties, params); err != nil {
				t.Fatal(err)
			}
			for _, testSet := range []func(tc *testContext, t *testing.T){
				testEncToShares,
				testRefresh,
				testRefreshAndPermutation,
				testRefreshAndTransformSwitchParams,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}
}

func gentestContext(nParties int, params bgv.Parameters) (tc *testContext, err error) {

	tc = new(testContext)

	tc.params = params

	tc.NParties = nParties

	tc.n = params.N()

	tc.ringT = params.RingT()
	tc.ringQ = params.RingQ()
	tc.ringP = params.RingP()

	prng, _ := sampling.NewKeyedPRNG([]byte{'t', 'e', 's', 't'})
	tc.crs = prng
	tc.uniformSampler = ring.NewUniformSampler(prng, params.RingQ())

	tc.encoder = bgv.NewEncoder(tc.params)
	tc.evaluator = bgv.NewEvaluator(tc.params, nil)

	kgen := rlwe.NewKeyGenerator(tc.params)

	// SecretKeys
	tc.sk0Shards = make([]*rlwe.SecretKey, nParties)
	tc.sk1Shards = make([]*rlwe.SecretKey, nParties)

	tc.sk0 = rlwe.NewSecretKey(tc.params.Parameters)
	tc.sk1 = rlwe.NewSecretKey(tc.params.Parameters)

	ringQP := params.RingQP()
	for j := 0; j < nParties; j++ {
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

func testEncToShares(tc *testContext, t *testing.T) {

	coeffs, _, ciphertext := newTestVectors(tc, tc.encryptorPk0, t)

	type Party struct {
		e2s         EncToShareProtocol
		s2e         ShareToEncProtocol
		sk          *rlwe.SecretKey
		publicShare multiparty.KeySwitchShare
		secretShare multiparty.AdditiveShare
	}

	params := tc.params
	P := make([]Party, tc.NParties)

	var err error
	for i := range P {
		if i == 0 {
			P[i].e2s, err = NewEncToShareProtocol(params, params.Xe())
			require.NoError(t, err)
			P[i].s2e, err = NewShareToEncProtocol(params, params.Xe())
			require.NoError(t, err)
		} else {
			P[i].e2s = P[0].e2s.ShallowCopy()
			P[i].s2e = P[0].s2e.ShallowCopy()
		}

		P[i].sk = tc.sk0Shards[i]
		P[i].publicShare = P[i].e2s.AllocateShare(ciphertext.Level())
		P[i].secretShare = NewAdditiveShare(params)
	}

	// The EncToShare protocol is run in all tests, as a setup to the ShareToEnc test.
	for i, p := range P {
		p.e2s.GenShare(p.sk, ciphertext, &p.secretShare, &p.publicShare)
		if i > 0 {
			p.e2s.AggregateShares(P[0].publicShare, p.publicShare, &P[0].publicShare)
		}
	}

	P[0].e2s.GetShare(&P[0].secretShare, P[0].publicShare, ciphertext, &P[0].secretShare)

	t.Run(GetTestName("EncToShareProtocol", tc.params, tc.NParties), func(t *testing.T) {

		rec := NewAdditiveShare(params)
		for _, p := range P {
			tc.ringT.Add(rec.Value, p.secretShare.Value, rec.Value)
		}

		ptRt := tc.params.RingT().NewPoly()
		ptRt.Copy(rec.Value)
		values := make([]uint64, len(coeffs))

		tc.encoder.DecodeRingT(ptRt, ciphertext.Scale, values)

		assert.True(t, slices.Equal(coeffs, values))
	})

	crp := P[0].e2s.SampleCRP(params.MaxLevel(), tc.crs)

	t.Run(GetTestName("ShareToEncProtocol", tc.params, tc.NParties), func(t *testing.T) {

		for i, p := range P {
			p.s2e.GenShare(p.sk, crp, p.secretShare, &p.publicShare)
			if i > 0 {
				p.s2e.AggregateShares(P[0].publicShare, p.publicShare, &P[0].publicShare)
			}
		}

		ctRec := bgv.NewCiphertext(tc.params, 1, tc.params.MaxLevel())
		*ctRec.MetaData = *ciphertext.MetaData
		P[0].s2e.GetEncryption(P[0].publicShare, crp, ctRec)

		verifyTestVectors(tc, tc.decryptorSk0, coeffs, ctRec, t)
	})
}

func testRefresh(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	encoder := tc.encoder
	decryptorSk0 := tc.decryptorSk0

	minLevel := 0
	maxLevel := tc.params.MaxLevel()

	t.Run(GetTestName("Refresh", tc.params, tc.NParties), func(t *testing.T) {

		type Party struct {
			RefreshProtocol
			s     *rlwe.SecretKey
			share multiparty.RefreshShare
		}

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			if i == 0 {
				var err error
				p.RefreshProtocol, err = NewRefreshProtocol(tc.params, tc.params.Xe())
				require.NoError(t, err)
			} else {
				p.RefreshProtocol = RefreshParties[0].RefreshProtocol.ShallowCopy()
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare(minLevel, maxLevel)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crp := P0.SampleCRP(maxLevel, tc.crs)

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)
		ciphertext.Resize(ciphertext.Degree(), minLevel)

		for i, p := range RefreshParties {
			p.GenShare(p.s, ciphertext, crp, &p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, &P0.share)
			}

		}

		P0.Finalize(ciphertext, crp, P0.share, ciphertext)

		//Decrypts and compare
		require.True(t, ciphertext.Level() == maxLevel)
		have := make([]uint64, tc.params.MaxSlots())
		encoder.Decode(decryptorSk0.DecryptNew(ciphertext), have)
		require.True(t, slices.Equal(coeffs, have))
	})
}

func testRefreshAndPermutation(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	encoder := tc.encoder
	decryptorSk0 := tc.decryptorSk0
	var err error

	minLevel := 0
	maxLevel := tc.params.MaxLevel()

	t.Run(GetTestName("RefreshAndPermutation", tc.params, tc.NParties), func(t *testing.T) {

		type Party struct {
			MaskedTransformProtocol
			s     *rlwe.SecretKey
			share multiparty.RefreshShare
		}

		RefreshParties := make([]*Party, tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			if i == 0 {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(tc.params, tc.params, tc.params.Xe()); err != nil {
					t.Fatal(err)
				}
			} else {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(tc.params, tc.params, tc.params.Xe()); err != nil {
					t.Fatal(err)
				}
			}

			p.s = sk0Shards[i]
			p.share = p.AllocateShare(minLevel, maxLevel)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crp := P0.SampleCRP(maxLevel, tc.crs)

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)
		ciphertext.Resize(ciphertext.Degree(), minLevel)

		permutation := make([]uint64, len(coeffs))
		N := uint64(len(coeffs))
		prng, _ := sampling.NewPRNG()
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

		maskedTransform := &MaskedTransformFunc{
			Decode: true,
			Func:   permute,
			Encode: true,
		}

		for i, p := range RefreshParties {
			p.GenShare(p.s, p.s, ciphertext, crp, maskedTransform, &p.share)
			if i > 0 {
				P0.AggregateShares(P0.share, p.share, &P0.share)
			}
		}

		P0.Transform(ciphertext, maskedTransform, crp, P0.share, ciphertext)

		coeffsPermute := make([]uint64, len(coeffs))
		for i := range coeffsPermute {
			coeffsPermute[i] = coeffs[permutation[i]]
		}

		coeffsHave := make([]uint64, tc.params.MaxSlots())
		encoder.Decode(decryptorSk0.DecryptNew(ciphertext), coeffsHave)

		//Decrypts and compares
		require.True(t, ciphertext.Level() == maxLevel)
		require.True(t, slices.Equal(coeffsPermute, coeffsHave))
	})
}

func testRefreshAndTransformSwitchParams(tc *testContext, t *testing.T) {

	encryptorPk0 := tc.encryptorPk0
	sk0Shards := tc.sk0Shards
	paramsIn := tc.params

	t.Run(GetTestName("RefreshAndTransformSwitchparams", tc.params, tc.NParties), func(t *testing.T) {

		var paramsOut bgv.Parameters
		var err error
		paramsOut, err = bgv.NewParametersFromLiteral(bgv.ParametersLiteral{
			LogN:             paramsIn.LogN(),
			LogQ:             []int{54, 49, 49, 49},
			LogP:             []int{52, 52},
			PlaintextModulus: paramsIn.PlaintextModulus(),
		})

		minLevel := 0
		maxLevel := paramsOut.MaxLevel()

		require.Nil(t, err)

		type Party struct {
			MaskedTransformProtocol
			sIn   *rlwe.SecretKey
			sOut  *rlwe.SecretKey
			share multiparty.RefreshShare
		}

		RefreshParties := make([]*Party, tc.NParties)
		kgenParamsOut := rlwe.NewKeyGenerator(paramsOut.Parameters)
		skIdealOut := rlwe.NewSecretKey(paramsOut.Parameters)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			if i == 0 {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(paramsIn, paramsOut, paramsIn.Xe()); err != nil {
					t.Fatal(err)
				}
			} else {
				if p.MaskedTransformProtocol, err = NewMaskedTransformProtocol(paramsIn, paramsOut, paramsIn.Xe()); err != nil {
					t.Fatal(err)
				}
			}

			p.sIn = sk0Shards[i]

			p.sOut = kgenParamsOut.GenSecretKeyNew() // New shared secret key in target parameters
			paramsOut.RingQ().Add(skIdealOut.Value.Q, p.sOut.Value.Q, skIdealOut.Value.Q)

			p.share = p.AllocateShare(minLevel, maxLevel)

			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crp := P0.SampleCRP(paramsOut.MaxLevel(), tc.crs)

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		permutation := make([]uint64, len(coeffs))
		N := uint64(len(coeffs))
		prng, _ := sampling.NewPRNG()
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
			p.GenShare(p.sIn, p.sOut, ciphertext, crp, transform, &p.share)
			if i > 0 {
				P0.AggregateShares(P0.share, p.share, &P0.share)
			}
		}

		P0.Transform(ciphertext, transform, crp, P0.share, ciphertext)

		transform.Func(coeffs)

		coeffsHave := make([]uint64, tc.params.MaxSlots())
		dec := rlwe.NewDecryptor(paramsOut.Parameters, skIdealOut)
		bgv.NewEncoder(paramsOut).Decode(dec.DecryptNew(ciphertext), coeffsHave)

		//Decrypts and compares
		require.True(t, ciphertext.Level() == maxLevel)
		require.True(t, slices.Equal(coeffs, coeffsHave))
	})
}

func newTestVectors(tc *testContext, encryptor *rlwe.Encryptor, t *testing.T) (coeffs []uint64, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {

	prng, err := sampling.NewPRNG()
	require.NoError(t, err)
	uniformSampler := ring.NewUniformSampler(prng, tc.ringT)
	coeffsPol := uniformSampler.ReadNew()

	for i := range coeffsPol.Coeffs[0] {
		coeffsPol.Coeffs[0][i] = uint64(1)
	}

	plaintext = bgv.NewPlaintext(tc.params, tc.params.MaxLevel())
	plaintext.Scale = tc.params.NewScale(2)
	require.NoError(t, tc.encoder.Encode(coeffsPol.Coeffs[0], plaintext))
	ciphertext, err = encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(tc *testContext, decryptor *rlwe.Decryptor, coeffs []uint64, ciphertext *rlwe.Ciphertext, t *testing.T) {
	have := make([]uint64, tc.params.MaxSlots())
	tc.encoder.Decode(decryptor.DecryptNew(ciphertext), have)
	require.True(t, slices.Equal(coeffs, have))
}
