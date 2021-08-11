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

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

// TestParams is a set of test parameters for the correctness of the rlwe pacakge.
var TestParams = []ParametersLiteral{TestPN12QP109, TestPN13QP218, TestPN14QP438, TestPN15QP880}

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

	for _, defaultParam := range defaultParams[:] {
		params, err := NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		kgen := NewKeyGenerator(params)

		for _, testSet := range []func(kgen KeyGenerator, t *testing.T){
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

func testGenKeyPair(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	// Checks that sum([-as + e, a] + [as])) <= N * 6 * sigma
	t.Run(testString(params, "PKGen/"), func(t *testing.T) {

		ringQ := params.RingQ()
		ringP := params.RingP()

		sk, pk := kgen.GenKeyPair()

		// [-as + e] + [as]
		ringQ.MulCoeffsMontgomeryAndAdd(sk.Value[0], pk.Value[1][0], pk.Value[0][0])
		ringQ.InvNTT(pk.Value[0][0], pk.Value[0][0])

		log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0][0].Level(), ringQ, pk.Value[0][0]))

		ringP.MulCoeffsMontgomeryAndAdd(sk.Value[1], pk.Value[1][1], pk.Value[0][1])
		ringP.InvNTT(pk.Value[0][1], pk.Value[0][1])
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0][1].Level(), ringP, pk.Value[0][1]))
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

		ringQ := params.RingQ()
		ringP := params.RingP()
		skIn := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		swk := NewSwitchingKey(params, params.QCount()-1, params.PCount()-1)
		kgen.(*keyGenerator).genSwitchingKey(skIn.Value[0], skOut.Value, swk)

		// Decrypts
		// [-asIn + w*P*sOut + e, a] + [asIn]
		for j := range swk.Value {
			ringQ.MulCoeffsMontgomeryAndAdd(swk.Value[j][1][0], skOut.Value[0], swk.Value[j][0][0])
			ringP.MulCoeffsMontgomeryAndAdd(swk.Value[j][1][1], skOut.Value[1], swk.Value[j][0][1])
		}

		polyQ := swk.Value[0][0][0]
		polyP := swk.Value[0][0][1]

		// Sums all basis together (equivalent to multiplying with CRT decomposition of 1)
		// sum([1]_w * [w*P*sOut + e]) = P*sOut + sum(e)
		for j := range swk.Value {
			if j > 0 {
				ringQ.Add(polyQ, swk.Value[j][0][0], polyQ)
				ringP.Add(polyP, swk.Value[j][0][1], polyP)
			}
		}

		// sOut * P
		ringQ.MulScalarBigint(skIn.Value[0], ringP.ModulusBigint, skIn.Value[0])

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(polyQ, skIn.Value[0], polyQ)

		// Checks that the error is below the bound
		// Worst error bound is N * floor(6*sigma) * #Keys
		ringQ.InvNTT(polyQ, polyQ)
		ringQ.InvMForm(polyQ, polyQ)
		ringP.InvNTT(polyP, polyP)
		ringP.InvMForm(polyP, polyP)

		log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()*len(swk.Value)))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringQ.Modulus)-1, ringQ, polyQ))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringP.Modulus)-1, ringP, polyP))

	})
}

func testEncryptor(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	ringQ := params.RingQ()

	t.Run(testString(params, "Encrypt/Pk/Fast/MaxLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewFastEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 12+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Pk/Fast/MinLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewFastEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 12+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Pk/Slow/MaxLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Pk/Slow/MinLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
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
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
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
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
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

		for i := 0; i < len(ks.PoolDecompQ); i++ {

			ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, ciphertext.Value[1], c2InvNTT, ks.PoolDecompQ[i], ks.PoolDecompP[i])

			// Compute q_alpha_i in bigInt
			qalphai := ring.NewInt(1)

			for j := 0; j < alpha; j++ {
				idx := i*alpha + j
				if idx > levelQ {
					break
				}
				qalphai.Mul(qalphai, ring.NewUint(ringQ.Modulus[idx]))
			}

			ringQ.ReduceLvl(levelQ, ks.PoolDecompQ[i], ks.PoolDecompQ[i])
			ringP.ReduceLvl(levelP, ks.PoolDecompP[i], ks.PoolDecompP[i])

			ringQ.InvNTTLvl(levelQ, ks.PoolDecompQ[i], tmpQ)
			ringP.InvNTTLvl(levelP, ks.PoolDecompP[i], tmpP)

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
		ks.SwitchKeysInPlace(ciphertext.Value[1].Level(), ciphertext.Value[1], swk, ks.PoolQ[1], ks.PoolQ[2])
		ringQ.Add(ciphertext.Value[0], ks.PoolQ[1], ciphertext.Value[0])
		ring.CopyValues(ks.PoolQ[2], ciphertext.Value[1])
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], skOut.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})
}

func testKeySwitchDimension(kgen KeyGenerator, t *testing.T) {

	paramsLargeDim := kgen.(*keyGenerator).params
	paramsSmallDim, _ := NewParametersFromLiteral(ParametersLiteral{
		LogN:  paramsLargeDim.LogN() - 1,
		Q:     paramsLargeDim.Q()[:1],
		P:     paramsLargeDim.P()[:1],
		Sigma: DefaultSigma,
	})

	t.Run(testString(paramsLargeDim, "KeySwitchDimension/LargeToSmall/"), func(t *testing.T) {

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

		ks := NewKeySwitcher(paramsLargeDim)
		ks.SwitchKeysInPlace(paramsSmallDim.MaxLevel(), ctLargeDim.Value[1], swk, ks.PoolQ[1], ks.PoolQ[2])
		ringQLargeDim.AddLvl(paramsSmallDim.MaxLevel(), ctLargeDim.Value[0], ks.PoolQ[1], ctLargeDim.Value[0])
		ring.CopyValues(ks.PoolQ[2], ctLargeDim.Value[1])

		//Extracts Coefficients
		ctSmallDim := NewCiphertextNTT(paramsSmallDim, 1, paramsSmallDim.MaxLevel())
		for i := range ctSmallDim.Value {
			ringQLargeDim.InvNTT(ctLargeDim.Value[i], ctLargeDim.Value[i])
			for j := range ctSmallDim.Value[i].Coeffs {
				tmp0 := ctSmallDim.Value[i].Coeffs[j]
				tmp1 := ctLargeDim.Value[i].Coeffs[j]
				gap := paramsLargeDim.N() / paramsSmallDim.N()
				for w := 0; w < paramsSmallDim.N(); w++ {
					tmp0[w] = tmp1[w*gap]
				}
			}
			ringQSmallDim.NTTLvl(ctSmallDim.Level(), ctSmallDim.Value[i], ctSmallDim.Value[i])
		}

		// Decrypts with smaller dimension key
		ringQSmallDim.MulCoeffsMontgomeryAndAddLvl(ctSmallDim.Level(), ctSmallDim.Value[1], skSmallDim.Value[0], ctSmallDim.Value[0])
		ringQSmallDim.InvNTTLvl(ctSmallDim.Level(), ctSmallDim.Value[0], ctSmallDim.Value[0])

		require.GreaterOrEqual(t, 10+paramsSmallDim.LogN(), log2OfInnerSum(ctSmallDim.Level(), ringQSmallDim, ctSmallDim.Value[0]))
	})

	t.Run(testString(paramsLargeDim, "KeySwitchDimension/SmallToLarge/"), func(t *testing.T) {

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

		ring.MapSmallDimensionToLargerDimensionNTT(ctSmallDim.Value[0], ctLargeDim.Value[0])
		ring.MapSmallDimensionToLargerDimensionNTT(ctSmallDim.Value[1], ctLargeDim.Value[1])

		ks := NewKeySwitcher(paramsLargeDim)
		ks.SwitchKeysInPlace(ctLargeDim.Value[1].Level(), ctLargeDim.Value[1], swk, ks.PoolQ[1], ks.PoolQ[2])
		ringQLargeDim.Add(ctLargeDim.Value[0], ks.PoolQ[1], ctLargeDim.Value[0])
		ring.CopyValues(ks.PoolQ[2], ctLargeDim.Value[1])

		// Decrypts with smaller dimension key
		ringQLargeDim.MulCoeffsMontgomeryAndAddLvl(ctLargeDim.Level(), ctLargeDim.Value[1], skLargeDim.Value[0], ctLargeDim.Value[0])
		ringQLargeDim.InvNTTLvl(ctLargeDim.Level(), ctLargeDim.Value[0], ctLargeDim.Value[0])

		require.GreaterOrEqual(t, 10+paramsSmallDim.LogN(), log2OfInnerSum(ctLargeDim.Level(), ringQLargeDim, ctLargeDim.Value[0]))
	})
}

func testMarshaller(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params
	ringQ := params.RingQ()
	ringP := params.RingP()

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

	t.Run(testString(params, "Marshaller/Sk/"), func(t *testing.T) {

		marshalledSk, err := sk.MarshalBinary()
		require.NoError(t, err)

		skTest := new(SecretKey)
		err = skTest.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, ringQ.Equal(sk.Value[0], skTest.Value[0]))
		require.True(t, ringP.Equal(sk.Value[1], skTest.Value[1]))

	})

	t.Run(testString(params, "Marshaller/Pk/"), func(t *testing.T) {

		marshalledPk, err := pk.MarshalBinary()
		require.NoError(t, err)

		pkTest := new(PublicKey)
		err = pkTest.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		for k := range pk.Value {
			require.Truef(t, ringQ.Equal(pk.Value[k][0], pkTest.Value[k][0]), "Marshal PublicKey element [%d][0]", k)
			require.Truef(t, ringP.Equal(pk.Value[k][1], pkTest.Value[k][1]), "Marshal PublicKey element [%d][1]", k)
		}
	})

	t.Run(testString(params, "Marshaller/EvaluationKey/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		evalKey := kgen.GenRelinearizationKey(sk, 3)
		data, err := evalKey.MarshalBinary()
		require.NoError(t, err)

		resEvalKey := new(RelinearizationKey)
		err = resEvalKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := evalKey.Keys[0].Value
		evakeyTest := resEvalKey.Keys[0].Value

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.Truef(t, ringQ.Equal(evakeyWant[j][k][0], evakeyTest[j][k][0]), "Marshal EvaluationKey element [%d][%d][0]", j, k)
				require.Truef(t, ringP.Equal(evakeyWant[j][k][1], evakeyTest[j][k][1]), "Marshal EvaluationKey element [%d][%d][1]", j, k)
			}
		}
	})

	t.Run(testString(params, "Marshaller/SwitchingKey/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		skOut := kgen.GenSecretKey()

		switchingKey := kgen.GenSwitchingKey(sk, skOut)
		data, err := switchingKey.MarshalBinary()
		require.NoError(t, err)

		resSwitchingKey := new(SwitchingKey)
		err = resSwitchingKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := switchingKey.Value
		evakeyTest := resSwitchingKey.Value

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.True(t, ringQ.Equal(evakeyWant[j][k][0], evakeyTest[j][k][0]), "Marshal SwitchingKey element [%d][%d][0]", j, k)
				require.True(t, ringP.Equal(evakeyWant[j][k][1], evakeyTest[j][k][1]), "Marshal SwitchingKey element [%d][%d][1]", j, k)
			}
		}
	})

	t.Run(testString(params, "Marshaller/RotationKey/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		rots := []int{1, -1, 63, -63}
		galEls := []uint64{params.GaloisElementForRowRotation()}
		for _, n := range rots {
			galEls = append(galEls, params.GaloisElementForColumnRotationBy(n))
		}

		rotationKey := kgen.GenRotationKeys(galEls, sk)

		data, err := rotationKey.MarshalBinary()
		require.NoError(t, err)

		resRotationKey := new(RotationKeySet)
		err = resRotationKey.UnmarshalBinary(data)
		require.NoError(t, err)

		for _, galEl := range galEls {

			evakeyWant := rotationKey.Keys[galEl].Value
			evakeyTest := resRotationKey.Keys[galEl].Value

			for j := range evakeyWant {
				for k := range evakeyWant[j] {
					require.Truef(t, ringQ.Equal(evakeyWant[j][k][0], evakeyTest[j][k][0]), "Marshal RotationKey RotateLeft %d element [%d][%d][0]", galEl, j, k)
					require.Truef(t, ringP.Equal(evakeyWant[j][k][1], evakeyTest[j][k][1]), "Marshal RotationKey RotateLeft %d element [%d][%d][1]", galEl, j, k)
				}
			}
		}
	})
}
