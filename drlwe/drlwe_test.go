package drlwe

import (
	"encoding/json"
	"flag"
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/stretchr/testify/require"
	"math"
	"math/big"
	"math/bits"
	"runtime"
	"testing"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func testString(params rlwe.Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

// TestParams is a set of test parameters for the correctness of the rlwe pacakge.
var TestParams = []rlwe.ParametersLiteral{rlwe.TestPN12QP109, rlwe.TestPN13QP218, rlwe.TestPN14QP438, rlwe.TestPN15QP880}

type testContext struct {
	params                 rlwe.Parameters
	kgen                   rlwe.KeyGenerator
	sk0, sk1, sk2, skIdeal *rlwe.SecretKey
	crpGenerator           UniformSampler
}

func newTestContext(params rlwe.Parameters) testContext {

	kgen := rlwe.NewKeyGenerator(params)
	sk0 := kgen.GenSecretKey()
	sk1 := kgen.GenSecretKey()
	sk2 := kgen.GenSecretKey()
	skIdeal := sk0.CopyNew()
	params.RingQ().Add(skIdeal.Value[0], sk1.Value[0], skIdeal.Value[0])
	params.RingP().Add(skIdeal.Value[1], sk1.Value[1], skIdeal.Value[1])
	params.RingQ().Add(skIdeal.Value[0], sk2.Value[0], skIdeal.Value[0])
	params.RingP().Add(skIdeal.Value[1], sk2.Value[1], skIdeal.Value[1])

	crpGenerator, _ := NewUniformSampler([]byte{}, params)

	return testContext{params, kgen, sk0, sk1, sk2, skIdeal, crpGenerator}
}

func TestDRLWE(t *testing.T) {
	defaultParams := TestParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = TestParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}

	if *flagParamString != "" {
		var jsonParams rlwe.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []rlwe.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams {
		params, err := rlwe.NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		textCtx := newTestContext(params)

		for _, testSet := range []func(textCtx testContext, t *testing.T){
			testPublicKeyGen,
			testKeySwitching,
			testPublicKeySwitching,
			testRelinKeyGen,
			testRotKeyGen,
			testMarshalling,
		} {
			testSet(textCtx, t)
			runtime.GC()
		}
	}
}

func testPublicKeyGen(testCtx testContext, t *testing.T) {

	params := testCtx.params
	ringQ := params.RingQ()
	ringP := params.RingP()

	t.Run(testString(params, "PublicKeyGen/"), func(t *testing.T) {

		CKGProtocol := NewCKGProtocol(params)

		share0 := CKGProtocol.AllocateShares()
		share1 := CKGProtocol.AllocateShares()
		share2 := CKGProtocol.AllocateShares()

		crp := testCtx.crpGenerator.ReadForCPKNew()

		CKGProtocol.GenShare(testCtx.sk0, crp, share0)
		CKGProtocol.GenShare(testCtx.sk1, crp, share1)
		CKGProtocol.GenShare(testCtx.sk2, crp, share2)

		CKGProtocol.AggregateShares(share0, share1, share0)
		CKGProtocol.AggregateShares(share0, share2, share0)

		pk := rlwe.NewPublicKey(params)
		CKGProtocol.GenPublicKey(share0, crp, pk)

		// [-as + e] + [as]
		ringQ.MulCoeffsMontgomeryAndAdd(testCtx.skIdeal.Value[0], pk.Value[1][0], pk.Value[0][0])
		ringP.MulCoeffsMontgomeryAndAdd(testCtx.skIdeal.Value[1], pk.Value[1][1], pk.Value[0][1])
		ringQ.InvNTT(pk.Value[0][0], pk.Value[0][0])
		ringP.InvNTT(pk.Value[0][1], pk.Value[0][1])

		log2Bound := bits.Len64(3 * uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0][0].Level(), ringQ, pk.Value[0][0]))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0][1].Level(), ringP, pk.Value[0][1]))
	})
}

func testKeySwitching(testCtx testContext, t *testing.T) {

	params := testCtx.params
	ringQ := params.RingQ()
	t.Run(testString(params, "KeySwitching/"), func(t *testing.T) {

		sk0Out := testCtx.kgen.GenSecretKey()
		sk1Out := testCtx.kgen.GenSecretKey()
		sk2Out := testCtx.kgen.GenSecretKey()

		skOutIdeal := sk0Out.Value[0].CopyNew()
		params.RingQ().Add(skOutIdeal, sk1Out.Value[0], skOutIdeal)
		params.RingQ().Add(skOutIdeal, sk2Out.Value[0], skOutIdeal)

		ciphertext := &rlwe.Ciphertext{Value: []*ring.Poly{ringQ.NewPoly(), testCtx.crpGenerator.ReadQNew()}}
		ringQ.MulCoeffsMontgomeryAndSub(ciphertext.Value[1], testCtx.skIdeal.Value[0], ciphertext.Value[0])
		ciphertext.Value[0].IsNTT = true
		ciphertext.Value[1].IsNTT = true

		CKSProtocol := NewCKSProtocol(params, rlwe.DefaultSigma)

		share0 := CKSProtocol.AllocateShare(ciphertext.Level())
		share1 := CKSProtocol.AllocateShare(ciphertext.Level())
		share2 := CKSProtocol.AllocateShare(ciphertext.Level())

		CKSProtocol.GenShare(testCtx.sk0, sk0Out, ciphertext, share0)
		CKSProtocol.GenShare(testCtx.sk1, sk1Out, ciphertext, share1)
		CKSProtocol.GenShare(testCtx.sk2, sk2Out, ciphertext, share2)

		CKSProtocol.AggregateShares(share0, share1, share0)
		CKSProtocol.AggregateShares(share0, share2, share0)

		ksCiphertext := &rlwe.Ciphertext{Value: []*ring.Poly{params.RingQ().NewPoly(), params.RingQ().NewPoly()}}

		CKSProtocol.KeySwitch(share0, ciphertext, ksCiphertext)

		// [-as + e] + [as]
		ringQ.MulCoeffsMontgomeryAndAdd(ksCiphertext.Value[1], skOutIdeal, ksCiphertext.Value[0])
		ringQ.InvNTT(ksCiphertext.Value[0], ksCiphertext.Value[0])
		log2Bound := bits.Len64(3 * uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(ksCiphertext.Value[0].Level(), ringQ, ksCiphertext.Value[0]))

	})
}

func testPublicKeySwitching(testCtx testContext, t *testing.T) {

	params := testCtx.params
	ringQ := params.RingQ()

	t.Run(testString(params, "PublicKeySwitching/"), func(t *testing.T) {

		skOut, pkOut := testCtx.kgen.GenKeyPair()

		ciphertext := &rlwe.Ciphertext{Value: []*ring.Poly{ringQ.NewPoly(), testCtx.crpGenerator.ReadQNew()}}
		ringQ.MulCoeffsMontgomeryAndSub(ciphertext.Value[1], testCtx.skIdeal.Value[0], ciphertext.Value[0])
		ciphertext.Value[0].IsNTT = true
		ciphertext.Value[1].IsNTT = true

		PCKSProtocol := NewPCKSProtocol(params, rlwe.DefaultSigma)

		share0 := PCKSProtocol.AllocateShare(ciphertext.Level())
		share1 := PCKSProtocol.AllocateShare(ciphertext.Level())
		share2 := PCKSProtocol.AllocateShare(ciphertext.Level())

		PCKSProtocol.GenShare(testCtx.sk0, pkOut, ciphertext, share0)
		PCKSProtocol.GenShare(testCtx.sk1, pkOut, ciphertext, share1)
		PCKSProtocol.GenShare(testCtx.sk2, pkOut, ciphertext, share2)

		PCKSProtocol.AggregateShares(share0, share1, share0)
		PCKSProtocol.AggregateShares(share0, share2, share0)

		ksCiphertext := &rlwe.Ciphertext{Value: []*ring.Poly{params.RingQ().NewPoly(), params.RingQ().NewPoly()}}

		PCKSProtocol.KeySwitch(share0, ciphertext, ksCiphertext)

		// [-as + e] + [as]
		ringQ.MulCoeffsMontgomeryAndAdd(ksCiphertext.Value[1], skOut.Value[0], ksCiphertext.Value[0])
		ringQ.InvNTT(ksCiphertext.Value[0], ksCiphertext.Value[0])
		log2Bound := bits.Len64(3 * uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound+5, log2OfInnerSum(ksCiphertext.Value[0].Level(), ringQ, ksCiphertext.Value[0]))

	})
}

func testRelinKeyGen(testCtx testContext, t *testing.T) {
	params := testCtx.params
	ringQ := params.RingQ()
	ringP := params.RingP()

	t.Run(testString(params, "RelinKeyGen/"), func(t *testing.T) {

		RKGProtocol := NewRKGProtocol(params, rlwe.DefaultSigma)

		ephSk0, share10, share20 := RKGProtocol.AllocateShares()
		ephSk1, share11, share21 := RKGProtocol.AllocateShares()
		ephSk2, share12, share22 := RKGProtocol.AllocateShares()

		crp := testCtx.crpGenerator.ReadForRKGNew()

		RKGProtocol.GenShareRoundOne(testCtx.sk0, crp, ephSk0, share10)
		RKGProtocol.GenShareRoundOne(testCtx.sk1, crp, ephSk1, share11)
		RKGProtocol.GenShareRoundOne(testCtx.sk2, crp, ephSk2, share12)

		RKGProtocol.AggregateShares(share10, share11, share10)
		RKGProtocol.AggregateShares(share10, share12, share10)

		RKGProtocol.GenShareRoundTwo(ephSk0, testCtx.sk0, share10, crp, share20)
		RKGProtocol.GenShareRoundTwo(ephSk1, testCtx.sk1, share10, crp, share21)
		RKGProtocol.GenShareRoundTwo(ephSk2, testCtx.sk2, share10, crp, share22)

		RKGProtocol.AggregateShares(share20, share21, share20)
		RKGProtocol.AggregateShares(share20, share22, share20)

		rlk := rlwe.NewRelinKey(params, 2)
		RKGProtocol.GenRelinearizationKey(share10, share20, rlk)

		skIn := testCtx.skIdeal.CopyNew()
		skOut := testCtx.skIdeal.CopyNew()
		ringQ.MulCoeffsMontgomery(skIn.Value[0], skIn.Value[0], skIn.Value[0])
		ringP.MulCoeffsMontgomery(skIn.Value[1], skIn.Value[1], skIn.Value[1])

		swk := rlk.Keys[0]

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
		ringP.InvNTT(polyP, polyP)
		ringQ.InvMForm(polyQ, polyQ)
		ringP.InvMForm(polyP, polyP)

		// Worst bound of inner sum
		// N*#Keys*(N * #Parties * floor(sigma*6) + #Parties * floor(sigma*6) + N * #Parties  +  #Parties * floor(6*sigma))
		log2Bound := bits.Len64(uint64(params.N() * len(swk.Value) * (params.N()*3*int(math.Floor(rlwe.DefaultSigma*6)) + 2*3*int(math.Floor(rlwe.DefaultSigma*6)) + params.N()*3)))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringQ.Modulus)-1, ringQ, polyQ))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringP.Modulus)-1, ringP, polyP))
	})
}

func testRotKeyGen(testCtx testContext, t *testing.T) {

	params := testCtx.params
	ringQ := params.RingQ()
	ringP := params.RingP()

	t.Run(testString(params, "RotKeyGen/"), func(t *testing.T) {

		crp := testCtx.crpGenerator.ReadForRTGNew()

		RTGProtocol := NewRTGProtocol(params)

		share0 := RTGProtocol.AllocateShares()
		share1 := RTGProtocol.AllocateShares()
		share2 := RTGProtocol.AllocateShares()

		galEl := params.GaloisElementForRowRotation()

		RTGProtocol.GenShare(testCtx.sk0, galEl, crp, share0)
		RTGProtocol.GenShare(testCtx.sk1, galEl, crp, share1)
		RTGProtocol.GenShare(testCtx.sk2, galEl, crp, share2)

		RTGProtocol.Aggregate(share0, share1, share0)
		RTGProtocol.Aggregate(share0, share2, share0)

		rotKeySet := rlwe.NewRotationKeySet(params, []uint64{galEl})
		RTGProtocol.GenRotationKey(share0, crp, rotKeySet.Keys[galEl])

		skIn := testCtx.skIdeal.CopyNew()
		skOut := testCtx.skIdeal.CopyNew()
		galElInv := ring.ModExp(galEl, int(4*params.N()-1), uint64(4*params.N()))
		ring.PermuteNTT(testCtx.skIdeal.Value[0], galElInv, skOut.Value[0])
		ring.PermuteNTT(testCtx.skIdeal.Value[1], galElInv, skOut.Value[1])

		swk := rotKeySet.Keys[galEl]

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

		// Worst bound of inner sum
		// N*#Keys*(N * #Parties * floor(sigma*6) + #Parties * floor(sigma*6) + N * #Parties  +  #Parties * floor(6*sigma))
		log2Bound := bits.Len64(3 * uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringQ.Modulus)-1, ringQ, polyQ))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringP.Modulus)-1, ringP, polyP))
	})
}

func testMarshalling(testCtx testContext, t *testing.T) {

	crs := testCtx.crpGenerator.ReadForCPKNew()
	params := testCtx.params

	ciphertext := &rlwe.Ciphertext{Value: []*ring.Poly{testCtx.crpGenerator.ReadQNew(), testCtx.crpGenerator.ReadQNew()}}

	t.Run(testString(params, "Marshalling/CPK/"), func(t *testing.T) {
		keygenProtocol := NewCKGProtocol(testCtx.params)
		KeyGenShareBefore := keygenProtocol.AllocateShares()
		keygenProtocol.GenShare(testCtx.sk0, crs, KeyGenShareBefore)
		//now we marshall it
		data, err := KeyGenShareBefore.MarshalBinary()

		if err != nil {
			t.Error("Could not marshal the CKGShare : ", err)
		}

		KeyGenShareAfter := new(CKGShare)
		err = KeyGenShareAfter.UnmarshalBinary(data)
		if err != nil {
			t.Error("Could not unmarshal the CKGShare : ", err)
		}

		//comparing the results
		require.Equal(t, KeyGenShareBefore.Value[0].Degree(), KeyGenShareAfter.Value[0].Degree())
		require.Equal(t, KeyGenShareBefore.Value[1].Degree(), KeyGenShareAfter.Value[1].Degree())
		require.Equal(t, KeyGenShareBefore.Value[0].LenModuli(), KeyGenShareAfter.Value[0].LenModuli())
		require.Equal(t, KeyGenShareBefore.Value[1].LenModuli(), KeyGenShareAfter.Value[1].LenModuli())

		require.Equal(t, KeyGenShareAfter.Value[0].Coeffs[:], KeyGenShareBefore.Value[0].Coeffs[:])
		require.Equal(t, KeyGenShareAfter.Value[1].Coeffs[:], KeyGenShareBefore.Value[1].Coeffs[:])
	})

	t.Run(testString(params, "Marshalling/PCKS/"), func(t *testing.T) {
		//Check marshalling for the PCKS

		KeySwitchProtocol := NewPCKSProtocol(testCtx.params, testCtx.params.Sigma())
		SwitchShare := KeySwitchProtocol.AllocateShare(ciphertext.Level())
		_, pkOut := testCtx.kgen.GenKeyPair()
		KeySwitchProtocol.GenShare(testCtx.sk0, pkOut, ciphertext, SwitchShare)

		data, err := SwitchShare.MarshalBinary()
		require.NoError(t, err)

		SwitchShareReceiver := new(PCKSShare)
		err = SwitchShareReceiver.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := 0; i < 2; i++ {
			//compare the shares.
			ringBefore := SwitchShare.Value[i]
			ringAfter := SwitchShareReceiver.Value[i]
			require.Equal(t, ringBefore.Degree(), ringAfter.Degree())
			require.Equal(t, ringAfter.Coeffs, ringBefore.Coeffs)
		}
	})

	t.Run(testString(params, "Marshalling/CKS/"), func(t *testing.T) {

		//Now for CKSShare ~ its similar to PKSShare
		cksp := NewCKSProtocol(testCtx.params, testCtx.params.Sigma())
		cksshare := cksp.AllocateShare(ciphertext.Level())
		cksp.GenShare(testCtx.sk0, testCtx.sk1, ciphertext, cksshare)

		data, err := cksshare.MarshalBinary()
		require.NoError(t, err)
		cksshareAfter := new(CKSShare)
		err = cksshareAfter.UnmarshalBinary(data)
		require.NoError(t, err)

		//now compare both shares.

		require.Equal(t, cksshare.Value.Degree(), cksshareAfter.Value.Degree())
		require.Equal(t, cksshare.Value.LenModuli(), cksshareAfter.Value.LenModuli())

		require.Equal(t, cksshare.Value.Coeffs, cksshareAfter.Value.Coeffs)
	})

	t.Run(testString(params, "Marshalling/RKG/"), func(t *testing.T) {

		//check RTGShare

		RKGProtocol := NewRKGProtocol(params, rlwe.DefaultSigma)

		ephSk0, share10, _ := RKGProtocol.AllocateShares()

		crp := testCtx.crpGenerator.ReadForRKGNew()

		RKGProtocol.GenShareRoundOne(testCtx.sk0, crp, ephSk0, share10)

		data, err := share10.MarshalBinary()
		require.NoError(t, err)

		rkgShare := new(RKGShare)
		err = rkgShare.UnmarshalBinary(data)
		require.NoError(t, err)

		require.Equal(t, len(rkgShare.Value), len(share10.Value))
		for i, val := range share10.Value {
			require.Equal(t, len(rkgShare.Value[i][0][0].Coeffs), len(val[0][0].Coeffs))
			require.Equal(t, len(rkgShare.Value[i][0][1].Coeffs), len(val[0][1].Coeffs))
			require.Equal(t, rkgShare.Value[i][0][0].Coeffs, val[0][0].Coeffs)
			require.Equal(t, rkgShare.Value[i][0][1].Coeffs, val[0][1].Coeffs)

			require.Equal(t, len(rkgShare.Value[i][1][0].Coeffs), len(val[1][0].Coeffs))
			require.Equal(t, len(rkgShare.Value[i][1][1].Coeffs), len(val[1][1].Coeffs))
			require.Equal(t, rkgShare.Value[i][1][0].Coeffs, val[1][0].Coeffs)
			require.Equal(t, rkgShare.Value[i][1][1].Coeffs, val[1][1].Coeffs)

		}
	})

	t.Run(testString(params, "Marshalling/RTG/"), func(t *testing.T) {

		//check RTGShare

		crp := testCtx.crpGenerator.ReadForRTGNew()

		galEl := testCtx.params.GaloisElementForColumnRotationBy(64)

		RTGProtocol := NewRTGProtocol(testCtx.params)
		rtgShare := RTGProtocol.AllocateShares()
		RTGProtocol.GenShare(testCtx.sk1, galEl, crp, rtgShare)

		data, err := rtgShare.MarshalBinary()
		require.NoError(t, err)

		resRTGShare := new(RTGShare)
		err = resRTGShare.UnmarshalBinary(data)
		require.NoError(t, err)

		require.Equal(t, len(resRTGShare.Value), len(rtgShare.Value))

		for i, val := range rtgShare.Value {
			require.Equal(t, len(resRTGShare.Value[i][0].Coeffs), len(val[0].Coeffs))
			require.Equal(t, resRTGShare.Value[i][0].Coeffs, val[0].Coeffs)
			require.Equal(t, len(resRTGShare.Value[i][1].Coeffs), len(val[1].Coeffs))
			require.Equal(t, resRTGShare.Value[i][1].Coeffs, val[1].Coeffs)

		}
	})
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
