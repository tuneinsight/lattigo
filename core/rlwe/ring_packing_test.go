package rlwe

import (
	"math"
	"math/big"
	"math/bits"
	"math/rand"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

const (
	LogNLarge = 10
	LogNSmall = 8
)

func TestRingPacking(t *testing.T) {

	var err error

	paramsLit := ParametersLiteral{
		LogN:    LogNLarge,
		LogQ:    []int{60},
		LogP:    []int{60},
		NTTFlag: true,
	}

	for _, NTTFlag := range []bool{true, false} {

		paramsLit.NTTFlag = NTTFlag

		var params Parameters
		if params, err = NewParametersFromLiteral(paramsLit); err != nil {
			t.Fatal(err)
		}

		tc, err := NewTestContext(params)
		require.NoError(t, err)

		for _, testSet := range []func(tc *TestContext, t *testing.T){
			testRingPacking,
		} {
			testSet(tc, t)
			runtime.GC()
		}
	}
}

func testRingPacking(tc *TestContext, t *testing.T) {

	params := tc.params
	sk := tc.sk
	enc := tc.enc
	dec := tc.dec
	level := params.MaxLevel()

	evkParams := EvaluationKeyParameters{
		LevelQ: utils.Pointy(params.MaxLevelQ()),
		LevelP: utils.Pointy(params.MaxLevelP()),
	}

	evkRP := RingPackingEvaluationKey{}

	ski, err := evkRP.GenRingSwitchingKeys(params, sk, LogNSmall, evkParams)
	require.NoError(t, err)

	evkRP.GenRepackEvaluationKeys(evkRP.Parameters[LogNSmall], ski[LogNSmall], evkParams)
	evkRP.GenRepackEvaluationKeys(evkRP.Parameters[params.LogN()], ski[params.LogN()], evkParams)
	evkRP.GenExtractEvaluationKeys(evkRP.Parameters[LogNSmall], ski[LogNSmall], evkParams)

	eval := NewRingPackingEvaluator(&evkRP)

	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), 0, "Split"), func(t *testing.T) {

		pt := genPlaintextNTT(params, level, 1<<40)
		ct, err := enc.EncryptNew(pt)
		require.NoError(t, err)

		ctEvenNHalf, ctOddNHalf, err := eval.SplitNew(ct)

		if eval.MaxLogN() == eval.MinLogN() {
			require.Error(t, err)
			t.Skip("eval.MaxLogN() = eval.MinLogN()")
		} else {
			require.NoError(t, err)

			paramsNHalf := eval.Parameters[ctEvenNHalf.LogN()].GetRLWEParameters()
			r := paramsNHalf.RingQ().AtLevel(ct.Level())

			decNHalf := NewDecryptor(paramsNHalf, ski[paramsNHalf.LogN()])

			ptEve := decNHalf.DecryptNew(ctEvenNHalf)
			ptOdd := decNHalf.DecryptNew(ctOddNHalf)

			if ptEve.IsNTT {
				r.INTT(ptEve.Value, ptEve.Value)
			}

			if ptOdd.IsNTT {
				r.INTT(ptOdd.Value, ptOdd.Value)
			}

			if pt.IsNTT {
				params.RingQ().AtLevel(ct.Level()).INTT(pt.Value, pt.Value)
			}

			for i := 0; i < level+1; i++ {

				Q := r.SubRings[i].Modulus
				ref := pt.Value.Coeffs[i]
				eve := ptEve.Value.Coeffs[i]
				odd := ptOdd.Value.Coeffs[i]

				for j := 0; j < paramsNHalf.N(); j++ {
					eve[j] = ring.CRed(eve[j]+Q-ref[j*2+0], Q)
					odd[j] = ring.CRed(odd[j]+Q-ref[j*2+1], Q)
				}
			}

			require.GreaterOrEqual(t, float64(paramsNHalf.LogN()+1), r.Log2OfStandardDeviation(ptEve.Value))
			require.GreaterOrEqual(t, float64(paramsNHalf.LogN()+1), r.Log2OfStandardDeviation(ptOdd.Value))
		}
	})

	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), 0, "Merge"), func(t *testing.T) {

		if eval.MaxLogN() == eval.MinLogN() {
			t.Skip("eval.MaxLogN() = eval.MinLogN()")
		}

		paramsNHalf := *eval.Parameters[params.LogN()-1].GetRLWEParameters()
		encNHalf := NewEncryptor(paramsNHalf, ski[paramsNHalf.LogN()])

		ptEve := genPlaintextNTT(paramsNHalf, level, 1<<40)
		ptOdd := genPlaintextNTT(paramsNHalf, level, 1<<40)

		ctEve, err := encNHalf.EncryptNew(ptEve)
		require.NoError(t, err)

		ctOdd, err := encNHalf.EncryptNew(ptOdd)
		require.NoError(t, err)

		ct, err := eval.MergeNew(ctEve, ctOdd)
		require.NoError(t, err)

		pt := dec.DecryptNew(ct)

		if pt.IsNTT {
			params.RingQ().AtLevel(level).INTT(pt.Value, pt.Value)
		}

		if ptEve.IsNTT {
			paramsNHalf.RingQ().AtLevel(level).INTT(ptEve.Value, ptEve.Value)
		}

		if ptOdd.IsNTT {
			paramsNHalf.RingQ().AtLevel(level).INTT(ptOdd.Value, ptOdd.Value)
		}

		for i := 0; i < level+1; i++ {
			Q := params.RingQ().SubRings[i].Modulus
			ref := pt.Value.Coeffs[i]
			eve := ptEve.Value.Coeffs[i]
			odd := ptOdd.Value.Coeffs[i]
			for j := 0; j < paramsNHalf.N(); j++ {
				ref[2*j+0] = ring.CRed(ref[2*j+0]+Q-eve[j], Q)
				ref[2*j+1] = ring.CRed(ref[2*j+1]+Q-odd[j], Q)
			}
		}

		require.GreaterOrEqual(t, float64(params.LogN()+1), params.RingQ().AtLevel(level).Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), 0, "Extract/Naive=False"), func(t *testing.T) {

		if params.RingType() != ring.Standard {
			t.Skip("Expand not supported for ring.Type = ring.ConjugateInvariant")
		}

		ringQ := params.RingQ().AtLevel(level)

		pt := genPlaintextNTT(params, level, 1<<40)

		ct, err := enc.EncryptNew(pt)
		require.NoError(t, err)

		gap := 17
		logGap := bits.Len64(uint64(gap))
		idx := map[int]bool{}
		for i := 0; i < params.N()/gap; i++ {
			idx[i*gap] = true
		}

		ciphertexts, err := eval.Extract(ct, idx)
		require.NoError(t, err)

		// Checks that the number of returned ciphertexts is equal
		// to the size of the index and that each element in the
		// index list has a corresponding extracted ciphertext.
		require.Equal(t, len(ciphertexts), len(idx))
		for i := range idx {
			_, ok := ciphertexts[i]
			require.True(t, ok)
		}

		// Decrypts & Checks
		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		paramsSmallN := evkRP.Parameters[ciphertexts[0].LogN()].GetRLWEParameters()

		ptDec := NewPlaintext(paramsSmallN, level)

		ringQSmallN := paramsSmallN.RingQ().AtLevel(level)
		Q := ringQSmallN.ModuliChain()

		decSmallN := NewDecryptor(paramsSmallN, ski[paramsSmallN.LogN()])

		for i := range idx {

			require.Equal(t, ciphertexts[i].LogN(), paramsSmallN.LogN())

			decSmallN.Decrypt(ciphertexts[i], ptDec)

			if ptDec.IsNTT {
				ringQSmallN.INTT(ptDec.Value, ptDec.Value)
			}

			for j := 0; j < level+1; j++ {
				ptDec.Value.Coeffs[j][0] = ring.CRed(ptDec.Value.Coeffs[j][0]+Q[j]-pt.Value.Coeffs[j][i], Q[j])
			}

			// Logs the noise
			require.GreaterOrEqual(t, float64(params.LogN()+logGap+1), ringQSmallN.Log2OfStandardDeviation(ptDec.Value))
		}
	})

	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), 0, "Extract/Naive=True"), func(t *testing.T) {

		if params.RingType() != ring.Standard {
			t.Skip("Expand not supported for ring.Type = ring.ConjugateInvariant")
		}

		ringQ := params.RingQ().AtLevel(level)

		pt := genPlaintextNTT(params, level, 1<<40)

		ct, err := enc.EncryptNew(pt)
		require.NoError(t, err)

		// Generates some extraction index map that contains
		// elements which are both not power and where the
		// smallest gap is not a power of two (to test the
		// worst case)
		gap := 17
		idx := map[int]bool{}
		for i := 0; i < params.N()/gap; i++ {
			idx[i*gap] = true
		}

		// Extract & returns a map containing the extracted RLWE ciphertexts.
		ciphertexts, err := eval.ExtractNaive(ct, idx)
		require.NoError(t, err)

		// Checks that the number of returned ciphertexts is equal
		// to the size of the index and that each element in the
		// index list has a corresponding extracted ciphertext.
		require.Equal(t, len(ciphertexts), len(idx))
		for i := range idx {
			_, ok := ciphertexts[i]
			require.True(t, ok)
		}

		// Decrypts & Checks
		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		paramsSmallN := evkRP.Parameters[ciphertexts[0].LogN()].GetRLWEParameters()

		ptDec := NewPlaintext(paramsSmallN, level)

		ringQSmallN := paramsSmallN.RingQ().AtLevel(level)
		Q := ringQSmallN.ModuliChain()

		decSmallN := NewDecryptor(paramsSmallN, ski[paramsSmallN.LogN()])

		for i := range idx {

			require.Equal(t, ciphertexts[i].LogN(), paramsSmallN.LogN())

			decSmallN.Decrypt(ciphertexts[i], ptDec)

			if ptDec.IsNTT {
				ringQSmallN.INTT(ptDec.Value, ptDec.Value)
			}

			for j := 0; j < level+1; j++ {
				ptDec.Value.Coeffs[j][0] = ring.CRed(ptDec.Value.Coeffs[j][0]+Q[j]-pt.Value.Coeffs[j][i], Q[j])
			}

			// Logs the noise
			coeffs := []*big.Int{new(big.Int)}
			ringQSmallN.PolyToBigintCentered(ptDec.Value, ringQ.N(), coeffs)
			noise := math.Log2(math.Abs(float64(coeffs[0].Int64())))

			require.GreaterOrEqual(t, float64(params.LogN()), noise)
		}
	})

	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), 0, "Repack"), func(t *testing.T) {

		if params.RingType() != ring.Standard {
			t.Skip("Pack not supported for ring.Type = ring.ConjugateInvariant")
		}

		pt := NewPlaintext(params, level)
		ringQ := tc.params.RingQ().AtLevel(level)

		ptPacked := genPlaintextNTT(params, level, 1<<40)
		ciphertexts := make(map[int]*Ciphertext)

		// Generates ciphertexts where the i-th ciphertext
		// having as constant coefficients the i-th coefficient
		// of the plaintext.
		// Generates a list of ciphertexts indexed by non-power-of-two
		// and where the smallest gap is not a power of two to test
		// the worst case.
		XInvNTT := GenXPow2NTT(ringQ, 1, true)[0]
		gap := 3
		for i := 0; i < params.N(); i++ {

			if i%gap == 0 {
				if ciphertexts[i], err = enc.EncryptNew(ptPacked); err != nil {
					t.Fatal(err)
				}
			}

			ringQ.MulCoeffsMontgomery(ptPacked.Value, XInvNTT, ptPacked.Value)
		}

		// Resets plaintext as it has been modified by being sequentially multiplied with X^-1
		ptPacked = genPlaintextNTT(params, level, 1<<40)

		// Repacks the ciphertexts
		ct, err := eval.Repack(ciphertexts)
		require.NoError(t, err)

		// Decrypts & Checks
		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		if ptPacked.IsNTT {
			ringQ.INTT(ptPacked.Value, ptPacked.Value)
		}

		for i := 0; i < level+1; i++ {
			Q := ringQ.SubRings[i].Modulus
			have := pt.Value.Coeffs[i]
			ref := ptPacked.Value.Coeffs[i]
			for j := 0; j < params.N(); j += gap {
				have[j] = ring.CRed(have[j]+Q-ref[j], Q)
			}
		}

		// Logs the noise
		require.GreaterOrEqual(t, float64(params.LogN()+5), ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), 0, "Extract[naive=false]->Permute->Repack[naive=true]"), func(t *testing.T) {
		testExtractPermuteRepack(params, level, enc, dec, eval, false, true, t)
	})

	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), 0, "Extract[naive=true]->Permute->Repack[naive=false]"), func(t *testing.T) {
		testExtractPermuteRepack(params, level, enc, dec, eval, true, false, t)
	})
}

func testExtractPermuteRepack(params Parameters, level int, enc *Encryptor, dec *Decryptor, eval *RingPackingEvaluator, ExtractNaive, RepackNaive bool, t *testing.T) {
	if params.RingType() != ring.Standard {
		t.Skip("Expand not supported for ring.Type = ring.ConjugateInvariant")
	}

	ringQ := params.RingQ().AtLevel(level)

	N := params.N()

	pt := genPlaintextNTT(params, level, 1<<40)

	ct, err := enc.EncryptNew(pt)
	require.NoError(t, err)

	// Ensures that ct is encrypted at the max
	// defined ring degree
	require.Equal(t, ct.LogN(), eval.MaxLogN())

	// Generates a random index selection
	// of size N/2 (to test that omitted
	// elements output zero coefficients)
	r := rand.New(rand.NewSource(0))
	list := make([]int, params.N())
	for i := range list {
		list[i] = i
	}
	r.Shuffle(len(list), func(i, j int) { list[i], list[j] = list[j], list[i] })

	idx := map[int]bool{}
	for _, i := range list[:params.N()>>1] {
		idx[i] = true
	}

	// Extract the coefficients at the given index
	var cts map[int]*Ciphertext
	if ExtractNaive {
		cts, err = eval.ExtractNaive(ct, idx)
	} else {
		cts, err = eval.Extract(ct, idx)
	}

	require.NoError(t, err)

	// Checks that the output ciphertext match the smallest
	// defined ring degree
	for i := range cts {
		require.Equal(t, cts[i].LogN(), eval.MinLogN())
	}

	// Defines a new mapping
	permute := func(x int) (y int) {
		return ((x + N/2) & (N - 1))
	}

	// Applies the mapping
	ctsPermute := map[int]*Ciphertext{}
	for i := range cts {
		ctsPermute[permute(i)] = cts[i]
	}

	// Repacks with the new permutation
	if RepackNaive {
		ct, err = eval.RepackNaive(ctsPermute)
	} else {
		ct, err = eval.Repack(ctsPermute)
	}
	require.NoError(t, err)

	// Decrypts & Checks
	ptHave := dec.DecryptNew(ct)

	if pt.IsNTT {
		ringQ.INTT(pt.Value, pt.Value)
	}

	if ptHave.IsNTT {
		ringQ.INTT(ptHave.Value, ptHave.Value)
	}

	for i := 0; i < level+1; i++ {
		Q := ringQ.SubRings[i].Modulus
		have := ptHave.Value.Coeffs[i]
		ref := pt.Value.Coeffs[i]
		for k0 := range idx {
			k1 := permute(k0)
			have[k1] = ring.CRed(have[k1]+Q-ref[k0], Q)
		}
	}

	// Logs the noise
	require.GreaterOrEqual(t, float64(params.LogN()+5), ringQ.Log2OfStandardDeviation(ptHave.Value))
}

func genPlaintextNTT(params Parameters, level, max int) (pt *Plaintext) {

	N := params.N()

	step := float64(max) / float64(N)

	pt = NewPlaintext(params, level)

	for i := 0; i < level+1; i++ {
		c := pt.Value.Coeffs[i]
		for j := 0; j < N; j++ {
			c[j] = uint64(float64(j) * step)
		}
	}

	params.RingQ().AtLevel(level).NTT(pt.Value, pt.Value)
	pt.IsNTT = true

	return
}
