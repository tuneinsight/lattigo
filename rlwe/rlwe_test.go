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
	"time"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

// TestParams is a set of test parameters for the correctness of the rlwe pacakge.
<<<<<<< HEAD
var TestParams = []ParametersLiteral{TestPN11QP54, TestPN12QP109, TestPN13QP218, TestPN14QP438, TestPN15QP880}
=======
var TestParams = []ParametersLiteral{TestPN12QP109, TestPN13QP218, TestPN14QP438, TestPN15QP880, TestPN16QP240, TestPN17QP360}
>>>>>>> dev_rckks

func testString(params Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

func TestRLWE(t *testing.T) {
	defaultParams := TestParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = TestParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams[:1] {
		params, err := NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		kgen := NewKeyGenerator(params)

		for _, testSet := range []func(kgen KeyGenerator, t *testing.T){
			testRGSW,
			testGenKeyPair,
			testSwitchKeyGen,
			testEncryptor,
			testDecryptor,
			testKeySwitcher,
			testKeySwitchDimension,
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
		var qi uint64
		var crtReconstruction *big.Int

		sumBigInt := ring.NewUint(0)
		QiB := new(big.Int)
		tmp := new(big.Int)
		modulusBigint := ring.NewUint(1)

		for i := 0; i < level+1; i++ {

			qi = ringQ.Modulus[i]
			QiB.SetUint64(qi)

			modulusBigint.Mul(modulusBigint, QiB)

			crtReconstruction = new(big.Int)
			crtReconstruction.Quo(ringQ.ModulusBigint, QiB)
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

func testRGSW(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	t.Run(testString(params, "RGSW/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		sk := kgen.GenSecretKey()
		ringQ := params.RingQ()

		encryptor := NewEncryptor(params, sk)

		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = false
		plaintext.Value.Coeffs[0][1] = 1 << 20

		ct := NewCiphertextNTT(params, 1, params.MaxLevel())
		encryptor.Encrypt(plaintext, ct)

		plaintextRGSW := NewPlaintext(params, params.MaxLevel())
		plaintextRGSW.Value.IsNTT = false
		plaintextRGSW.Value.Coeffs[0][0] = 1

		rgswciphertext := make([]*RGSWCiphertext, 1<<12)
		for i := 0; i < 1; i++ {
			fmt.Println(i)
			rgswciphertext[i] = NewCiphertextRGSWNTT(params, params.MaxLevel())
			encryptor.(*skEncryptor).EncryptRGSW(plaintextRGSW, rgswciphertext[i])

		}

		ks := NewKeySwitcher(params)

		ctOut := NewCiphertextNTT(params, 1, params.MaxLevel())

		now := time.Now()
		for i := 0; i < 1<<10; i++ {
			ks.MulRGSW(ct, rgswciphertext[0], ctOut)
		}
		fmt.Printf("Done: %s", time.Since(now))

		ringQ.MulCoeffsMontgomeryAndAddLvl(ctOut.Level(), ctOut.Value[1], sk.Value.Q, ctOut.Value[0])
		ringQ.InvNTTLvl(ctOut.Level(), ctOut.Value[0], ctOut.Value[0])
		fmt.Println(ctOut.Value[0].Coeffs[0][:4])
	})

}

func testGenKeyPair(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	// Checks that sum([-as + e, a] + [as])) <= N * 6 * sigma
	t.Run(testString(params, "PKGen/"), func(t *testing.T) {
		sk, pk := kgen.GenKeyPair()

		// [-as + e] + [as]
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
	t.Run(testString(params, "SWKGen/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		ringQ := params.RingQ()
		ringP := params.RingP()
		ringQP := params.RingQP()
		skIn := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()
		levelQ, levelP := params.QCount()-1, params.PCount()-1

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		swk := NewSwitchingKey(params, params.QCount()-1, params.PCount()-1)
		kgen.(*keyGenerator).genSwitchingKey(skIn.Value.Q, skOut.Value, swk)

		// Decrypts
		// [-asIn + w*P*sOut + e, a] + [asIn]
		for j := range swk.Value {
			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, swk.Value[j][1], skOut.Value, swk.Value[j][0])
		}

		// Sums all basis together (equivalent to multiplying with CRT decomposition of 1)
		// sum([1]_w * [w*P*sOut + e]) = P*sOut + sum(e)
		for j := range swk.Value {
			if j > 0 {
				ringQP.AddLvl(levelQ, levelP, swk.Value[0][0], swk.Value[j][0], swk.Value[0][0])
			}
		}

		// sOut * P
		ringQ.MulScalarBigint(skIn.Value.Q, ringP.ModulusBigint, skIn.Value.Q)

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(swk.Value[0][0].Q, skIn.Value.Q, swk.Value[0][0].Q)

		// Checks that the error is below the bound
		// Worst error bound is N * floor(6*sigma) * #Keys

		ringQP.InvNTTLvl(levelQ, levelP, swk.Value[0][0], swk.Value[0][0])
		ringQP.InvMFormLvl(levelQ, levelP, swk.Value[0][0], swk.Value[0][0])

		log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()*len(swk.Value)))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringQ.Modulus)-1, ringQ, swk.Value[0][0].Q))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringP.Modulus)-1, ringP, swk.Value[0][0].P))

	})
}

func testEncryptor(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	ringQ := params.RingQ()

	t.Run(testString(params, "Encrypt/Pk/MaxLevel/"), func(t *testing.T) {
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

	t.Run(testString(params, "Encrypt/Pk/MinLevel/"), func(t *testing.T) {
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

	t.Run(testString(params, "Encrypt/Sk/MaxLevel/"), func(t *testing.T) {
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

	t.Run(testString(params, "Encrypt/Sk/MinLevel/"), func(t *testing.T) {
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
}

func testDecryptor(kgen KeyGenerator, t *testing.T) {
	params := kgen.(*keyGenerator).params
	sk := kgen.GenSecretKey()
	ringQ := params.RingQ()
	encryptor := NewEncryptor(params, sk)
	decryptor := NewDecryptor(params, sk)

	t.Run(testString(params, "Decrypt/MaxLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})

	t.Run(testString(params, "Encrypt/MinLevel/"), func(t *testing.T) {
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

	t.Run(testString(params, "KeySwitch/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		sk := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()
		ks := NewKeySwitcher(params)

		ringQ := params.RingQ()
		ringP := params.RingP()

		levelQ := params.MaxLevel()
		alpha := params.PCount()
		levelP := alpha - 1

		QBig := ring.NewUint(1)
		for i := range ringQ.Modulus[:levelQ+1] {
			QBig.Mul(QBig, ring.NewUint(ringQ.Modulus[i]))
		}

		PBig := ring.NewUint(1)
		for i := range ringP.Modulus[:levelP+1] {
			PBig.Mul(PBig, ring.NewUint(ringP.Modulus[i]))
		}

		plaintext := NewPlaintext(params, levelQ)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, sk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)

		// Tests that a random polynomial decomposed is equal to its
		// reconstruction mod each RNS
		t.Run(testString(params, "DecomposeNTT/"), func(t *testing.T) {

			c2InvNTT := ringQ.NewPolyLvl(ciphertext.Level())
			ringQ.InvNTT(ciphertext.Value[1], c2InvNTT)

			coeffsBigintHaveQ := make([]*big.Int, ringQ.N)
			coeffsBigintHaveP := make([]*big.Int, ringQ.N)
			coeffsBigintRef := make([]*big.Int, ringQ.N)
			coeffsBigintWant := make([]*big.Int, ringQ.N)

			for i := range coeffsBigintRef {
				coeffsBigintHaveQ[i] = new(big.Int)
				coeffsBigintHaveP[i] = new(big.Int)
				coeffsBigintRef[i] = new(big.Int)
				coeffsBigintWant[i] = new(big.Int)
			}

			ringQ.PolyToBigintCenteredLvl(ciphertext.Level(), c2InvNTT, coeffsBigintRef)

			tmpQ := ringQ.NewPolyLvl(ciphertext.Level())
			tmpP := ringP.NewPolyLvl(levelP)

			for i := 0; i < len(ks.PoolDecompQP); i++ {

				ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, ciphertext.Value[1], c2InvNTT, ks.PoolDecompQP[i].Q, ks.PoolDecompQP[i].P)

				// Compute q_alpha_i in bigInt
				qalphai := ring.NewInt(1)

				for j := 0; j < alpha; j++ {
					idx := i*alpha + j
					if idx > levelQ {
						break
					}
					qalphai.Mul(qalphai, ring.NewUint(ringQ.Modulus[idx]))
				}

				ringQ.ReduceLvl(levelQ, ks.PoolDecompQP[i].Q, ks.PoolDecompQP[i].Q)
				ringP.ReduceLvl(levelP, ks.PoolDecompQP[i].P, ks.PoolDecompQP[i].P)

				ringQ.InvNTTLvl(levelQ, ks.PoolDecompQP[i].Q, tmpQ)
				ringP.InvNTTLvl(levelP, ks.PoolDecompQP[i].P, tmpP)

				ringQ.PolyToBigintCenteredLvl(levelQ, tmpQ, coeffsBigintHaveQ)
				ringP.PolyToBigintCenteredLvl(levelP, tmpP, coeffsBigintHaveP)

				// Checks that Reconstruct(NTT(c2 mod Q)) mod q_alpha_i == Reconstruct(NTT(Decomp(c2 mod Q, q_alpha-i) mod QP))
				for i := range coeffsBigintWant[:1] {

					coeffsBigintWant[i].Mod(coeffsBigintRef[i], qalphai)
					coeffsBigintWant[i].Mod(coeffsBigintWant[i], QBig)
					coeffsBigintHaveQ[i].Mod(coeffsBigintHaveQ[i], QBig)
					require.Equal(t, coeffsBigintHaveQ[i].Cmp(coeffsBigintWant[i]), 0)

					coeffsBigintWant[i].Mod(coeffsBigintRef[i], qalphai)
					coeffsBigintWant[i].Mod(coeffsBigintWant[i], PBig)
					coeffsBigintHaveP[i].Mod(coeffsBigintHaveP[i], PBig)
					require.Equal(t, coeffsBigintHaveP[i].Cmp(coeffsBigintWant[i]), 0)

				}
			}
		})

		// Test that Dec(KS(Enc(ct, sk), skOut), skOut) has a small norm
		t.Run(testString(params, "KeySwitch/Standard/"), func(t *testing.T) {
			swk := kgen.GenSwitchingKey(sk, skOut)
			ks.SwitchKeysInPlace(ciphertext.Value[1].Level(), ciphertext.Value[1], swk, ks.Pool[1].Q, ks.Pool[2].Q)
			ringQ.Add(ciphertext.Value[0], ks.Pool[1].Q, ciphertext.Value[0])
			ring.CopyValues(ks.Pool[2].Q, ciphertext.Value[1])
			ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], skOut.Value.Q, ciphertext.Value[0])
			ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
			require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
		})
	})
}

func testKeySwitchDimension(kgen KeyGenerator, t *testing.T) {

	paramsLargeDim := kgen.(*keyGenerator).params
<<<<<<< HEAD
=======
	paramsSmallDim, _ := NewParametersFromLiteral(ParametersLiteral{
		LogN:     paramsLargeDim.LogN() - 1,
		Q:        paramsLargeDim.Q()[:1],
		P:        paramsLargeDim.P()[:1],
		Sigma:    DefaultSigma,
		RingType: paramsLargeDim.RingType(),
	})
>>>>>>> dev_rckks

	t.Run(testString(paramsLargeDim, "KeySwitchDimension/"), func(t *testing.T) {

		if paramsLargeDim.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		paramsSmallDim, _ := NewParametersFromLiteral(ParametersLiteral{
			LogN:     paramsLargeDim.LogN() - 1,
			Q:        paramsLargeDim.Q()[:1],
			P:        paramsLargeDim.P()[:1],
			Sigma:    DefaultSigma,
			RingType: paramsLargeDim.RingType(),
		})

		t.Run("LargeToSmall/", func(t *testing.T) {

			ringQLargeDim := paramsLargeDim.RingQ()
			ringQSmallDim := paramsSmallDim.RingQ()

			kgenLargeDim := NewKeyGenerator(paramsLargeDim)
			skLargeDim := kgenLargeDim.GenSecretKey()
			kgenSmallDim := NewKeyGenerator(paramsSmallDim)
			skSmallDim := kgenSmallDim.GenSecretKey()

			swk := kgenLargeDim.GenSwitchingKey(skLargeDim, skSmallDim)

<<<<<<< HEAD
			plaintext := NewPlaintext(paramsLargeDim, paramsLargeDim.MaxLevel())
			plaintext.Value.IsNTT = true
			encryptor := NewEncryptor(paramsLargeDim, skLargeDim)
			ctLargeDim := NewCiphertextNTT(paramsLargeDim, 1, plaintext.Level())
			encryptor.Encrypt(plaintext, ctLargeDim)
=======
		SwitchCiphertextRingDegreeNTT(ctLargeDim, ringQSmallDim, ringQLargeDim, ctSmallDim)
>>>>>>> dev_rckks

			ks := NewKeySwitcher(paramsLargeDim)
			ks.SwitchKeysInPlace(paramsSmallDim.MaxLevel(), ctLargeDim.Value[1], swk, ks.Pool[1].Q, ks.Pool[2].Q)
			ringQLargeDim.AddLvl(paramsSmallDim.MaxLevel(), ctLargeDim.Value[0], ks.Pool[1].Q, ctLargeDim.Value[0])
			ring.CopyValues(ks.Pool[2].Q, ctLargeDim.Value[1])

			//Extracts Coefficients
			ctSmallDim := NewCiphertextNTT(paramsSmallDim, 1, paramsSmallDim.MaxLevel())

			SwitchCiphertextRingDegreeNTT(ctLargeDim, ringQSmallDim, ringQLargeDim, ctSmallDim)

			// Decrypts with smaller dimension key
			ringQSmallDim.MulCoeffsMontgomeryAndAddLvl(ctSmallDim.Level(), ctSmallDim.Value[1], skSmallDim.Value.Q, ctSmallDim.Value[0])
			ringQSmallDim.InvNTTLvl(ctSmallDim.Level(), ctSmallDim.Value[0], ctSmallDim.Value[0])

			require.GreaterOrEqual(t, 10+paramsSmallDim.LogN(), log2OfInnerSum(ctSmallDim.Level(), ringQSmallDim, ctSmallDim.Value[0]))
		})

		t.Run("SmallToLarge/", func(t *testing.T) {

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

<<<<<<< HEAD
			//Extracts Coefficients
			ctLargeDim := NewCiphertextNTT(paramsLargeDim, 1, plaintext.Level())
=======
		SwitchCiphertextRingDegreeNTT(ctSmallDim, nil, nil, ctLargeDim)
>>>>>>> dev_rckks

			SwitchCiphertextRingDegreeNTT(ctSmallDim, nil, nil, ctLargeDim)

			ks := NewKeySwitcher(paramsLargeDim)
			ks.SwitchKeysInPlace(ctLargeDim.Value[1].Level(), ctLargeDim.Value[1], swk, ks.Pool[1].Q, ks.Pool[2].Q)
			ringQLargeDim.Add(ctLargeDim.Value[0], ks.Pool[1].Q, ctLargeDim.Value[0])
			ring.CopyValues(ks.Pool[2].Q, ctLargeDim.Value[1])

			// Decrypts with smaller dimension key
			ringQLargeDim.MulCoeffsMontgomeryAndAddLvl(ctLargeDim.Level(), ctLargeDim.Value[1], skLargeDim.Value.Q, ctLargeDim.Value[0])
			ringQLargeDim.InvNTTLvl(ctLargeDim.Level(), ctLargeDim.Value[0], ctLargeDim.Value[0])

			require.GreaterOrEqual(t, 10+paramsSmallDim.LogN(), log2OfInnerSum(ctLargeDim.Level(), ringQLargeDim, ctLargeDim.Value[0]))
		})
	})
}

func testMarshaller(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	t.Run("Marshaller/Parameters/Binary", func(t *testing.T) {
		bytes, err := params.MarshalBinary()
		assert.Nil(t, err)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		assert.Nil(t, err)
		assert.Equal(t, params, p)
		assert.Equal(t, params.RingQ(), p.RingQ())
	})

	t.Run("Marshaller/Parameters/JSON", func(t *testing.T) {
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

	t.Run(testString(params, "Marshaller/Ciphertext/"), func(t *testing.T) {

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

	t.Run(testString(params, "Marshaller/Sk/"), func(t *testing.T) {

		marshalledSk, err := sk.MarshalBinary()
		require.NoError(t, err)

		skTest := new(SecretKey)
		err = skTest.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, sk.Value.Equals(skTest.Value))
	})

	t.Run(testString(params, "Marshaller/Pk/"), func(t *testing.T) {

		marshalledPk, err := pk.MarshalBinary()
		require.NoError(t, err)

		pkTest := new(PublicKey)
		err = pkTest.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		require.True(t, pk.Equals(pkTest))
	})

	t.Run(testString(params, "Marshaller/EvaluationKey/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		evalKey := kgen.GenRelinearizationKey(sk, 3)
		data, err := evalKey.MarshalBinary()
		require.NoError(t, err)

		resEvalKey := new(RelinearizationKey)
		err = resEvalKey.UnmarshalBinary(data)
		require.NoError(t, err)

		require.True(t, evalKey.Equals(resEvalKey))
	})

	t.Run(testString(params, "Marshaller/SwitchingKey/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		skOut := kgen.GenSecretKey()

		switchingKey := kgen.GenSwitchingKey(sk, skOut)
		data, err := switchingKey.MarshalBinary()
		require.NoError(t, err)

		resSwitchingKey := new(SwitchingKey)
		err = resSwitchingKey.UnmarshalBinary(data)
		require.NoError(t, err)

		require.True(t, switchingKey.Equals(resSwitchingKey))
	})

	t.Run(testString(params, "Marshaller/RotationKey/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		rots := []int{1, -1, 63, -63}
		galEls := []uint64{}
<<<<<<< HEAD
		if params.RingType() == RingStandard {
=======
		if params.RingType() == ring.Standard {
>>>>>>> dev_rckks
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
