package rlwe

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

// TestParams is a set of test parameters for the correctness of the rlwe pacakge.
var TestParams = []ParametersLiteral{TestPN10QP27, TestPN11QP54, TestPN12QP109, TestPN13QP218, TestPN14QP438, TestPN15QP880, TestPN16QP240, TestPN17QP360}

func testString(params Parameters, opname string) string {
	return fmt.Sprintf("%s/logN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

func TestRLWE(t *testing.T) {

	var err error

	defaultParams := TestParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = TestParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams[:] {
		var params Parameters
		if params, err = NewParametersFromLiteral(defaultParam); err != nil {
			t.Fatal(err)
		}

		kgen := NewKeyGenerator(params)

		for _, testSet := range []func(kgen KeyGenerator, t *testing.T){
			testGenKeyPair,
			testSwitchKeyGen,
			testEncryptor,
			testDecryptor,
			testKeySwitcher,
			testKeySwitchDimension,
			testMergeRLWE,
			testExpandRLWE,
			testMarshaller,
		} {
			testSet(kgen, t)
			runtime.GC()
		}
	}
}

// Returns the ceil(log2) of the sum of the absolute value of all the coefficients
func log2OfInnerSum(level int, ringQ *ring.Ring, poly *ring.Poly) (logSum int) {
	sumRNS := make([]uint64, level+1)
	var sum uint64
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		qiHalf := qi >> 1
		coeffs := poly.Coeffs[i]
		sum = 0

		for j := 0; j < ringQ.N; j++ {

			v := coeffs[j]

			if v >= qiHalf {
				sum = ring.CRed(sum+qi-v, qi)
			} else {
				sum = ring.CRed(sum+v, qi)
			}
		}

		sumRNS[i] = sum
	}

	var smallNorm = true
	for i := 1; i < level+1; i++ {
		smallNorm = smallNorm && (sumRNS[0] == sumRNS[i])
	}

	if !smallNorm {
		var crtReconstruction *big.Int

		sumBigInt := ring.NewUint(0)
		QiB := new(big.Int)
		tmp := new(big.Int)
		modulusBigint := ringQ.ModulusAtLevel[level]

		for i := 0; i < level+1; i++ {
			QiB.SetUint64(ringQ.Modulus[i])
			crtReconstruction = new(big.Int)
			crtReconstruction.Quo(modulusBigint, QiB)
			tmp.ModInverse(crtReconstruction, QiB)
			tmp.Mod(tmp, QiB)
			crtReconstruction.Mul(crtReconstruction, tmp)

			sumBigInt.Add(sumBigInt, tmp.Mul(ring.NewUint(sumRNS[i]), crtReconstruction))
		}

		sumBigInt.Mod(sumBigInt, modulusBigint)

		logSum = sumBigInt.BitLen()
	} else {
		logSum = bits.Len64(sumRNS[0])
	}

	return
}

func testGenKeyPair(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	// Checks that the secret-key has exactly params.h non-zero coefficients
	t.Run(testString(params, "SK"), func(t *testing.T) {

		skInvNTT := NewSecretKey(params)

		if params.PCount() > 0 {
			params.RingP().InvNTTLvl(sk.Value.P.Level(), sk.Value.P, skInvNTT.Value.P)
			for i := range skInvNTT.Value.P.Coeffs {
				var zeros int
				for j := range skInvNTT.Value.P.Coeffs[i] {
					if skInvNTT.Value.P.Coeffs[i][j] == 0 {
						zeros++
					}
				}
				require.Equal(t, params.ringP.N, zeros+params.h)
			}
		}

		params.RingQ().InvNTTLvl(sk.Value.Q.Level(), sk.Value.Q, skInvNTT.Value.Q)
		for i := range skInvNTT.Value.Q.Coeffs {
			var zeros int
			for j := range skInvNTT.Value.Q.Coeffs[i] {
				if skInvNTT.Value.Q.Coeffs[i][j] == 0 {
					zeros++
				}
			}
			require.Equal(t, params.ringQ.N, zeros+params.h)
		}

	})

	// Checks that sum([-as + e, a] + [as])) <= N * 6 * sigma
	t.Run(testString(params, "PK"), func(t *testing.T) {

		if params.PCount() > 0 {

			params.RingQP().MulCoeffsMontgomeryAndAddLvl(sk.Value.Q.Level(), sk.Value.P.Level(), sk.Value, pk.Value[1], pk.Value[0])
			params.RingQP().InvNTTLvl(sk.Value.Q.Level(), sk.Value.P.Level(), pk.Value[0], pk.Value[0])

			log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].Q.Level(), params.RingQ(), pk.Value[0].Q))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].P.Level(), params.RingP(), pk.Value[0].P))
		} else {
			params.RingQ().MulCoeffsMontgomeryAndAdd(sk.Value.Q, pk.Value[1].Q, pk.Value[0].Q)
			params.RingQ().InvNTT(pk.Value[0].Q, pk.Value[0].Q)

			log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].Q.Level(), params.RingQ(), pk.Value[0].Q))
		}

	})

}

func testSwitchKeyGen(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	// Checks that switching keys are en encryption under the output key
	// of the RNS decomposition of the input key by
	// 1) Decrypting the RNS decomposed input key
	// 2) Reconstructing the key
	// 3) Checking that the difference with the input key has a small norm
	t.Run(testString(params, "SWKGen"), func(t *testing.T) {

		ringQ := params.RingQ()
		ringP := params.RingP()
		ringQP := params.RingQP()
		skIn := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()
		levelQ, levelP := params.QCount()-1, params.PCount()-1
		decompBIT := params.DecompBIT(levelQ, levelP)

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		swk := NewSwitchingKey(params, params.QCount()-1, params.PCount()-1)
		kgen.(*keyGenerator).genSwitchingKey(skIn.Value.Q, skOut.Value, swk)

		// Decrypts
		// [-asIn + w*P*sOut + e, a] + [asIn]
		for i := range swk.Value {
			for _, el := range swk.Value[i] {
				ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, el[1], skOut.Value, el[0])
			}
		}

		// Sums all basis together (equivalent to multiplying with CRT decomposition of 1)
		// sum([1]_w * [RNS*BIT*P*sOut + e]) = BIT*P*sOut + sum(e)
		for i := range swk.Value { // RNS decomp
			if i > 0 {
				for j := range swk.Value[i] { // BIT decomp
					ringQP.AddLvl(levelQ, levelP, swk.Value[0][j][0], swk.Value[i][j][0], swk.Value[0][j][0])
				}
			}
		}

		// sOut * P * BIT
		if levelP != -1 {
			ringQ.MulScalarBigint(skIn.Value.Q, ringP.ModulusAtLevel[levelP], skIn.Value.Q)
		}

		log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()*len(swk.Value)))
		for i := 0; i < decompBIT; i++ {

			// P*s^i + sum(e) - P*s^i = sum(e)
			ringQ.Sub(swk.Value[0][i][0].Q, skIn.Value.Q, swk.Value[0][i][0].Q)

			// Checks that the error is below the bound
			// Worst error bound is N * floor(6*sigma) * #Keys

			ringQP.InvNTTLvl(levelQ, levelP, swk.Value[0][i][0], swk.Value[0][i][0])
			ringQP.InvMFormLvl(levelQ, levelP, swk.Value[0][i][0], swk.Value[0][i][0])

			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(levelQ, ringQ, swk.Value[0][i][0].Q))

			if levelP != -1 {
				require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(levelP, ringP, swk.Value[0][i][0].P))
			}

			// sOut * P * BIT
			ringQ.MulScalar(skIn.Value.Q, 1<<params.logbase2, skIn.Value.Q)
		}
	})
}

func testEncryptor(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	ringQ := params.RingQ()

	t.Run(testString(params, "Encrypt/Pk/MaxLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Pk/MinLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Pk/ShallowCopy"), func(t *testing.T) {
		enc1 := NewEncryptor(params, pk)
		enc2 := enc1.ShallowCopy()
		pkEnc1, pkEnc2 := enc1.(*pkEncryptor), enc2.(*pkEncryptor)
		require.True(t, pkEnc1.params.Equals(pkEnc2.params))
		require.True(t, pkEnc1.pk == pkEnc2.pk)
		require.False(t, (pkEnc1.basisextender == pkEnc2.basisextender) && (pkEnc1.basisextender != nil) && (pkEnc2.basisextender != nil))
		require.False(t, pkEnc1.encryptorBuffers == pkEnc2.encryptorBuffers)
		require.False(t, pkEnc1.ternarySampler == pkEnc2.ternarySampler)
		require.False(t, pkEnc1.gaussianSampler == pkEnc2.gaussianSampler)
	})

	t.Run(testString(params, "Encrypt/Sk/MaxLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, sk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Sk/MinLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, sk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Sk/ShallowCopy"), func(t *testing.T) {
		enc1 := NewEncryptor(params, sk)
		enc2 := enc1.ShallowCopy()
		skEnc1, skEnc2 := enc1.(*skEncryptor), enc2.(*skEncryptor)
		require.True(t, skEnc1.params.Equals(skEnc2.params))
		require.True(t, skEnc1.sk == skEnc2.sk)
		require.False(t, (skEnc1.basisextender == skEnc2.basisextender) && (skEnc1.basisextender != nil) && (skEnc2.basisextender != nil))
		require.False(t, skEnc1.encryptorBuffers == skEnc2.encryptorBuffers)
		require.False(t, skEnc1.ternarySampler == skEnc2.ternarySampler)
		require.False(t, skEnc1.gaussianSampler == skEnc2.gaussianSampler)
	})

	t.Run(testString(params, "Encrypt/WithKey/Sk->Sk"), func(t *testing.T) {
		sk2 := kgen.GenSecretKey()
		enc1 := NewEncryptor(params, sk)
		enc2 := enc1.WithKey(sk2)
		skEnc1, skEnc2 := enc1.(*skEncryptor), enc2.(*skEncryptor)
		require.True(t, skEnc1.params.Equals(skEnc2.params))
		require.True(t, skEnc1.sk.Value.Equals(sk.Value))
		require.True(t, skEnc2.sk.Value.Equals(sk2.Value))
		require.True(t, skEnc1.basisextender == skEnc2.basisextender)
		require.True(t, skEnc1.encryptorBuffers == skEnc2.encryptorBuffers)
		require.True(t, skEnc1.ternarySampler == skEnc2.ternarySampler)
		require.True(t, skEnc1.gaussianSampler == skEnc2.gaussianSampler)
	})
}

func testDecryptor(kgen KeyGenerator, t *testing.T) {
	params := kgen.(*keyGenerator).params
	sk := kgen.GenSecretKey()
	ringQ := params.RingQ()
	encryptor := NewEncryptor(params, sk)
	decryptor := NewDecryptor(params, sk)

	t.Run(testString(params, "Decrypt/MaxLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})

	t.Run(testString(params, "Encrypt/MinLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})
}

func testKeySwitcher(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	t.Run(testString(params, "KeySwitch"), func(t *testing.T) {

		sk := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()
		eval := NewEvaluator(params, nil)

		ringQ := params.RingQ()

		levelQ := params.MaxLevel()

		plaintext := NewPlaintext(params, levelQ)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, sk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)

		// Test that Dec(KS(Enc(ct, sk), skOut), skOut) has a small norm
		swk := kgen.GenSwitchingKey(sk, skOut)
		eval.GadgetProduct(ciphertext.Value[1].Level(), ciphertext.Value[1], swk.Ciphertext, eval.BuffQP[1].Q, eval.BuffQP[2].Q)
		ringQ.Add(ciphertext.Value[0], eval.BuffQP[1].Q, ciphertext.Value[0])
		ring.CopyValues(eval.BuffQP[2].Q, ciphertext.Value[1])
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], skOut.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 11+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})
}

func testKeySwitchDimension(kgen KeyGenerator, t *testing.T) {

	paramsLargeDim := kgen.(*keyGenerator).params

	t.Run(testString(paramsLargeDim, "KeySwitchDimension"), func(t *testing.T) {

		var Q []uint64
		if len(paramsLargeDim.Q()) > 1 {
			Q = paramsLargeDim.Q()[:1]
		} else {
			Q = paramsLargeDim.Q()
		}

		var P []uint64
		if len(paramsLargeDim.P()) != 0 {
			P = paramsLargeDim.P()[:1]
		} else {
			P = []uint64{}
		}

		paramsSmallDim, err := NewParametersFromLiteral(ParametersLiteral{
			LogN:     paramsLargeDim.LogN() - 1,
			Q:        Q,
			P:        P,
			Sigma:    DefaultSigma,
			RingType: paramsLargeDim.RingType(),
		})

		assert.Nil(t, err)

		t.Run("LargeToSmall", func(t *testing.T) {

			ringQLargeDim := paramsLargeDim.RingQ()
			ringQSmallDim := paramsSmallDim.RingQ()

			kgenLargeDim := NewKeyGenerator(paramsLargeDim)
			skLargeDim := kgenLargeDim.GenSecretKey()
			kgenSmallDim := NewKeyGenerator(paramsSmallDim)
			skSmallDim := kgenSmallDim.GenSecretKey()

			swk := kgenLargeDim.GenSwitchingKey(skLargeDim, skSmallDim)

			plaintext := NewPlaintext(paramsLargeDim, paramsLargeDim.MaxLevel())
			plaintext.Value.IsNTT = true
			encryptor := NewEncryptor(paramsLargeDim, skLargeDim)
			ctLargeDim := NewCiphertextNTT(paramsLargeDim, 1, plaintext.Level())
			encryptor.Encrypt(plaintext, ctLargeDim)

			eval := NewEvaluator(paramsLargeDim, nil)
			eval.GadgetProduct(paramsSmallDim.MaxLevel(), ctLargeDim.Value[1], swk.Ciphertext, eval.BuffQP[1].Q, eval.BuffQP[2].Q)
			ringQLargeDim.AddLvl(paramsSmallDim.MaxLevel(), ctLargeDim.Value[0], eval.BuffQP[1].Q, ctLargeDim.Value[0])
			ring.CopyValues(eval.BuffQP[2].Q, ctLargeDim.Value[1])

			//Extracts Coefficients
			ctSmallDim := NewCiphertextNTT(paramsSmallDim, 1, paramsSmallDim.MaxLevel())

			SwitchCiphertextRingDegreeNTT(ctLargeDim, ringQSmallDim, ringQLargeDim, ctSmallDim)

			// Decrypts with smaller dimension key
			ringQSmallDim.MulCoeffsMontgomeryAndAddLvl(ctSmallDim.Level(), ctSmallDim.Value[1], skSmallDim.Value.Q, ctSmallDim.Value[0])
			ringQSmallDim.InvNTTLvl(ctSmallDim.Level(), ctSmallDim.Value[0], ctSmallDim.Value[0])

			require.GreaterOrEqual(t, 11+paramsSmallDim.LogN(), log2OfInnerSum(ctSmallDim.Level(), ringQSmallDim, ctSmallDim.Value[0]))
		})

		t.Run("SmallToLarge", func(t *testing.T) {

			ringQLargeDim := paramsLargeDim.RingQ()

			kgenLargeDim := NewKeyGenerator(paramsLargeDim)
			skLargeDim := kgenLargeDim.GenSecretKey()
			kgenSmallDim := NewKeyGenerator(paramsSmallDim)
			skSmallDim := kgenSmallDim.GenSecretKey()

			swk := kgenLargeDim.GenSwitchingKey(skSmallDim, skLargeDim)

			plaintext := NewPlaintext(paramsSmallDim, paramsSmallDim.MaxLevel())
			plaintext.Value.IsNTT = true
			encryptor := NewEncryptor(paramsSmallDim, skSmallDim)
			ctSmallDim := NewCiphertextNTT(paramsSmallDim, 1, plaintext.Level())
			encryptor.Encrypt(plaintext, ctSmallDim)

			//Extracts Coefficients
			ctLargeDim := NewCiphertextNTT(paramsLargeDim, 1, plaintext.Level())

			SwitchCiphertextRingDegreeNTT(ctSmallDim, nil, nil, ctLargeDim)

			eval := NewEvaluator(paramsLargeDim, nil)
			eval.GadgetProduct(ctLargeDim.Value[1].Level(), ctLargeDim.Value[1], swk.Ciphertext, eval.BuffQP[1].Q, eval.BuffQP[2].Q)
			ringQLargeDim.Add(ctLargeDim.Value[0], eval.BuffQP[1].Q, ctLargeDim.Value[0])
			ring.CopyValues(eval.BuffQP[2].Q, ctLargeDim.Value[1])

			// Decrypts with smaller dimension key
			ringQLargeDim.MulCoeffsMontgomeryAndAddLvl(ctLargeDim.Level(), ctLargeDim.Value[1], skLargeDim.Value.Q, ctLargeDim.Value[0])
			ringQLargeDim.InvNTTLvl(ctLargeDim.Level(), ctLargeDim.Value[0], ctLargeDim.Value[0])

			require.GreaterOrEqual(t, 11+paramsSmallDim.LogN(), log2OfInnerSum(ctLargeDim.Level(), ringQLargeDim, ctLargeDim.Value[0]))
		})
	})
}

func testMergeRLWE(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	t.Run(testString(params, "MergeRLWE"), func(t *testing.T) {

		kgen := NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encryptor := NewEncryptor(params, sk)
		decryptor := NewDecryptor(params, sk)
		pt := NewPlaintext(params, params.MaxLevel())

		for i := 0; i < pt.Level()+1; i++ {
			for j := 0; j < params.N(); j++ {
				pt.Value.Coeffs[i][j] = (1 << 30) + uint64(j)*(1<<20)
			}
		}

		ciphertexts := make(map[int]*Ciphertext)
		slotIndex := make(map[int]bool)
		for i := 0; i < params.N(); i += params.N() / 16 {
			ciphertexts[i] = NewCiphertextNTT(params, 1, params.MaxLevel())
			encryptor.Encrypt(pt, ciphertexts[i])
			slotIndex[i] = true
		}

		// Rotation Keys
		galEls := params.GaloisElementsForMergeRLWE()
		rtks := kgen.GenRotationKeys(galEls, sk)

		eval := NewEvaluator(params, &EvaluationKey{Rtks: rtks})

		ciphertext := eval.MergeRLWE(ciphertexts)

		decryptor.Decrypt(ciphertext, pt)

		bound := uint64(params.N() * params.N())

		for i := 0; i < pt.Level()+1; i++ {

			Q := params.RingQ().Modulus[i]
			QHalf := Q >> 1

			for i, c := range pt.Value.Coeffs[i] {

				if c >= QHalf {
					c = Q - c
				}

				if _, ok := slotIndex[i]; !ok {
					if c > bound {
						t.Fatal(i, c)
					}
				}
			}
		}
	})
}

func testExpandRLWE(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	t.Run(testString(params, "ExpandRLWE"), func(t *testing.T) {

		kgen := NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encryptor := NewEncryptor(params, sk)
		decryptor := NewDecryptor(params, sk)
		pt := NewPlaintext(params, params.MaxLevel())

		logN := 3

		values := make([]uint64, params.N())

		scale := 1 << 22

		for i := 0; i < 1<<logN; i++ {
			values[i] = uint64(scale * i)
		}

		for i := 0; i < pt.Level()+1; i++ {
			copy(pt.Value.Coeffs[i], values)
		}

		ctIn := NewCiphertextNTT(params, 1, params.MaxLevel())
		encryptor.Encrypt(pt, ctIn)

		// Rotation Keys
		galEls := params.GaloisElementForExpandRLWE(logN)

		rtks := kgen.GenRotationKeys(galEls, sk)

		eval := NewEvaluator(params, &EvaluationKey{Rtks: rtks})

		ciphertexts := eval.ExpandRLWE(ctIn, logN)

		bound := uint64(params.N() * params.N())

		for i := range ciphertexts {

			decryptor.Decrypt(ciphertexts[i], pt)

			for j := 0; j < pt.Level()+1; j++ {

				Q := params.RingQ().Modulus[j]
				QHalf := Q >> 1

				for k, c := range pt.Value.Coeffs[j] {

					if c >= QHalf {
						c = Q - c
					}

					if k != 0 {
						require.Greater(t, bound, c)
					} else {
						require.InDelta(t, 0, math.Abs(float64(values[i])-float64(c))/float64(scale), 0.01)
					}
				}
			}
		}
	})
}

func testMarshaller(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	t.Run(testString(params, "Marshaller/Parameters/Binary"), func(t *testing.T) {
		bytes, err := params.MarshalBinary()
		assert.Nil(t, err)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		assert.Nil(t, err)
		assert.Equal(t, params, p)
		assert.Equal(t, params.RingQ(), p.RingQ())
	})

	t.Run(testString(params, "Marshaller/Parameters/JSON"), func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(params)
		assert.Nil(t, err)
		assert.NotNil(t, data)

		// checks that Parameters can be unmarshalled without error
		var rlweParams Parameters
		err = json.Unmarshal(data, &rlweParams)
		assert.Nil(t, err)
		assert.True(t, params.Equals(rlweParams))
	})

	t.Run(testString(params, "Marshaller/Ciphertext"), func(t *testing.T) {

		prng, _ := utils.NewPRNG()

		for degree := 0; degree < 4; degree++ {
			t.Run(fmt.Sprintf("degree=%d", degree), func(t *testing.T) {
				ciphertextWant := NewCiphertextRandom(prng, params, degree, params.MaxLevel())

				marshalledCiphertext, err := ciphertextWant.MarshalBinary()
				require.NoError(t, err)

				ciphertextTest := new(Ciphertext)
				require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

				require.Equal(t, ciphertextWant.Degree(), ciphertextTest.Degree())
				require.Equal(t, ciphertextWant.Level(), ciphertextTest.Level())

				for i := range ciphertextWant.Value {
					require.True(t, params.RingQ().EqualLvl(ciphertextWant.Level(), ciphertextWant.Value[i], ciphertextTest.Value[i]))
				}
			})
		}
	})

	t.Run(testString(params, "Marshaller/Sk"), func(t *testing.T) {

		marshalledSk, err := sk.MarshalBinary()
		require.NoError(t, err)

		skTest := new(SecretKey)
		err = skTest.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, sk.Value.Equals(skTest.Value))
	})

	t.Run(testString(params, "Marshaller/Pk"), func(t *testing.T) {

		marshalledPk, err := pk.MarshalBinary()
		require.NoError(t, err)

		pkTest := new(PublicKey)
		err = pkTest.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		require.True(t, pk.Equals(pkTest))
	})

	t.Run(testString(params, "Marshaller/EvaluationKey"), func(t *testing.T) {

		evalKey := kgen.GenRelinearizationKey(sk, 3)
		data, err := evalKey.MarshalBinary()
		require.NoError(t, err)

		resEvalKey := new(RelinearizationKey)
		err = resEvalKey.UnmarshalBinary(data)
		require.NoError(t, err)

		require.True(t, evalKey.Equals(resEvalKey))
	})

	t.Run(testString(params, "Marshaller/SwitchingKey"), func(t *testing.T) {

		skOut := kgen.GenSecretKey()

		switchingKey := kgen.GenSwitchingKey(sk, skOut)
		data, err := switchingKey.MarshalBinary()
		require.NoError(t, err)

		resSwitchingKey := new(SwitchingKey)
		err = resSwitchingKey.UnmarshalBinary(data)
		require.NoError(t, err)

		require.True(t, switchingKey.Equals(resSwitchingKey))
	})

	t.Run(testString(params, "Marshaller/RotationKey"), func(t *testing.T) {

		rots := []int{1, -1, 63, -63}
		galEls := []uint64{}
		if params.RingType() == ring.Standard {
			galEls = append(galEls, params.GaloisElementForRowRotation())
		}

		for _, n := range rots {
			galEls = append(galEls, params.GaloisElementForColumnRotationBy(n))
		}

		rotationKey := kgen.GenRotationKeys(galEls, sk)

		data, err := rotationKey.MarshalBinary()
		require.NoError(t, err)

		resRotationKey := new(RotationKeySet)
		err = resRotationKey.UnmarshalBinary(data)
		require.NoError(t, err)

		rotationKey.Equals(resRotationKey)
	})
}
