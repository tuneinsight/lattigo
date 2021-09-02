package drlwe

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
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/require"
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
	uniformSampler         *ring.UniformSampler
	crs                    utils.PRNG
}

func newTestContext(params rlwe.Parameters) testContext {

	kgen := rlwe.NewKeyGenerator(params)
	sk0 := kgen.GenSecretKey()
	sk1 := kgen.GenSecretKey()
	sk2 := kgen.GenSecretKey()
	skIdeal := sk0.CopyNew()
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	params.RingQP().AddLvl(levelQ, levelP, skIdeal.Value, sk1.Value, skIdeal.Value)
	params.RingQP().AddLvl(levelQ, levelP, skIdeal.Value, sk2.Value, skIdeal.Value)

	prng, _ := utils.NewKeyedPRNG([]byte{'t', 'e', 's', 't'})
	unifSampler := ring.NewUniformSampler(prng, params.RingQ())

	return testContext{params, kgen, sk0, sk1, sk2, skIdeal, unifSampler, prng}
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
	ringQP := params.RingQP()
	levelQ, levelP := params.QCount()-1, params.PCount()-1

	t.Run(testString(params, "PublicKeyGen/"), func(t *testing.T) {

		ckg := NewCKGProtocol(params)

		share0 := ckg.AllocateShares()
		share1 := ckg.AllocateShares()
		share2 := ckg.AllocateShares()

		crp := ckg.SampleCRP(testCtx.crs)

		ckg.GenShare(testCtx.sk0, crp, share0)
		ckg.GenShare(testCtx.sk1, crp, share1)
		ckg.GenShare(testCtx.sk2, crp, share2)

		ckg.AggregateShares(share0, share1, share0)
		ckg.AggregateShares(share0, share2, share0)

		pk := rlwe.NewPublicKey(params)
		ckg.GenPublicKey(share0, crp, pk)

		// [-as + e] + [as]
		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, testCtx.skIdeal.Value, pk.Value[1], pk.Value[0])
		ringQP.InvNTTLvl(levelQ, levelP, pk.Value[0], pk.Value[0])

		log2Bound := bits.Len64(3 * uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].Q.Level(), ringQ, pk.Value[0].Q))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].P.Level(), ringP, pk.Value[0].P))
	})
}

func testKeySwitching(testCtx testContext, t *testing.T) {

	params := testCtx.params
	ringQ := params.RingQ()
	ringQP := params.RingQP()
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	t.Run(testString(params, "KeySwitching/"), func(t *testing.T) {

		sk0Out := testCtx.kgen.GenSecretKey()
		sk1Out := testCtx.kgen.GenSecretKey()
		sk2Out := testCtx.kgen.GenSecretKey()

		skOutIdeal := sk0Out.CopyNew()
		ringQP.AddLvl(levelQ, levelP, skOutIdeal.Value, sk1Out.Value, skOutIdeal.Value)
		ringQP.AddLvl(levelQ, levelP, skOutIdeal.Value, sk2Out.Value, skOutIdeal.Value)

		ciphertext := &rlwe.Ciphertext{Value: []*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}}
		testCtx.uniformSampler.Read(ciphertext.Value[1])
		ringQ.MulCoeffsMontgomeryAndSub(ciphertext.Value[1], testCtx.skIdeal.Value.Q, ciphertext.Value[0])
		ciphertext.Value[0].IsNTT = true
		ciphertext.Value[1].IsNTT = true

		cks := NewCKSProtocol(params, rlwe.DefaultSigma)

		share0 := cks.AllocateShare(ciphertext.Level())
		share1 := cks.AllocateShare(ciphertext.Level())
		share2 := cks.AllocateShare(ciphertext.Level())

		cks.GenShare(testCtx.sk0, sk0Out, ciphertext, share0)
		cks.GenShare(testCtx.sk1, sk1Out, ciphertext, share1)
		cks.GenShare(testCtx.sk2, sk2Out, ciphertext, share2)

		cks.AggregateShares(share0, share1, share0)
		cks.AggregateShares(share0, share2, share0)

		ksCiphertext := &rlwe.Ciphertext{Value: []*ring.Poly{params.RingQ().NewPoly(), params.RingQ().NewPoly()}}

		cks.KeySwitch(share0, ciphertext, ksCiphertext)

		// [-as + e] + [as]
		ringQ.MulCoeffsMontgomeryAndAdd(ksCiphertext.Value[1], skOutIdeal.Value.Q, ksCiphertext.Value[0])
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

		ciphertext := &rlwe.Ciphertext{Value: []*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}}
		testCtx.uniformSampler.Read(ciphertext.Value[1])
		ringQ.MulCoeffsMontgomeryAndSub(ciphertext.Value[1], testCtx.skIdeal.Value.Q, ciphertext.Value[0])
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
		ringQ.MulCoeffsMontgomeryAndAdd(ksCiphertext.Value[1], skOut.Value.Q, ksCiphertext.Value[0])
		ringQ.InvNTT(ksCiphertext.Value[0], ksCiphertext.Value[0])
		log2Bound := bits.Len64(3 * uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound+5, log2OfInnerSum(ksCiphertext.Value[0].Level(), ringQ, ksCiphertext.Value[0]))

	})
}

func testRelinKeyGen(testCtx testContext, t *testing.T) {
	params := testCtx.params
	ringQ := params.RingQ()
	ringP := params.RingP()
	ringQP := params.RingQP()
	levelQ, levelP := params.QCount()-1, params.PCount()-1

	t.Run(testString(params, "RelinKeyGen/"), func(t *testing.T) {

		rkg := NewRKGProtocol(params, rlwe.DefaultSigma)

		ephSk0, share10, share20 := rkg.AllocateShares()
		ephSk1, share11, share21 := rkg.AllocateShares()
		ephSk2, share12, share22 := rkg.AllocateShares()

		crp := rkg.SampleCRP(testCtx.crs)

		rkg.GenShareRoundOne(testCtx.sk0, crp, ephSk0, share10)
		rkg.GenShareRoundOne(testCtx.sk1, crp, ephSk1, share11)
		rkg.GenShareRoundOne(testCtx.sk2, crp, ephSk2, share12)

		rkg.AggregateShares(share10, share11, share10)
		rkg.AggregateShares(share10, share12, share10)

		rkg.GenShareRoundTwo(ephSk0, testCtx.sk0, share10, crp, share20)
		rkg.GenShareRoundTwo(ephSk1, testCtx.sk1, share10, crp, share21)
		rkg.GenShareRoundTwo(ephSk2, testCtx.sk2, share10, crp, share22)

		rkg.AggregateShares(share20, share21, share20)
		rkg.AggregateShares(share20, share22, share20)

		rlk := rlwe.NewRelinKey(params, 2)
		rkg.GenRelinearizationKey(share10, share20, rlk)

		skIn := testCtx.skIdeal.CopyNew()
		skOut := testCtx.skIdeal.CopyNew()
		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, skIn.Value, skIn.Value, skIn.Value)

		swk := rlk.Keys[0]

		// Decrypts
		// [-asIn + w*P*sOut + e, a] + [asIn]
		for j := range swk.Value {
			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, swk.Value[j][1], skOut.Value, swk.Value[j][0])
		}

		poly := swk.Value[0][0]

		// Sums all basis together (equivalent to multiplying with CRT decomposition of 1)
		// sum([1]_w * [w*P*sOut + e]) = P*sOut + sum(e)
		for j := range swk.Value {
			if j > 0 {
				ringQP.AddLvl(levelQ, levelP, poly, swk.Value[j][0], poly)
			}
		}

		// sOut * P
		ringQ.MulScalarBigint(skIn.Value.Q, ringP.ModulusBigint, skIn.Value.Q)

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(swk.Value[0][0].Q, skIn.Value.Q, swk.Value[0][0].Q)

		// Checks that the error is below the bound
		// Worst error bound is N * floor(6*sigma) * #Keys
		ringQP.InvNTTLvl(levelQ, levelP, poly, poly)
		ringQP.InvMFormLvl(levelQ, levelP, poly, poly)

		// Worst bound of inner sum
		// N*#Keys*(N * #Parties * floor(sigma*6) + #Parties * floor(sigma*6) + N * #Parties  +  #Parties * floor(6*sigma))
		log2Bound := bits.Len64(uint64(params.N() * len(swk.Value) * (params.N()*3*int(math.Floor(rlwe.DefaultSigma*6)) + 2*3*int(math.Floor(rlwe.DefaultSigma*6)) + params.N()*3)))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringQ.Modulus)-1, ringQ, swk.Value[0][0].Q))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringP.Modulus)-1, ringP, swk.Value[0][0].P))
	})
}

func testRotKeyGen(testCtx testContext, t *testing.T) {

	params := testCtx.params
	ringQ := params.RingQ()
	ringP := params.RingP()
	ringQP := params.RingQP()
	levelQ, levelP := params.QCount()-1, params.PCount()-1

	t.Run(testString(params, "RotKeyGen/"), func(t *testing.T) {

		rtg := NewRTGProtocol(params)

		share0 := rtg.AllocateShares()
		share1 := rtg.AllocateShares()
		share2 := rtg.AllocateShares()

		crp := rtg.SampleCRP(testCtx.crs)

		galEl := params.GaloisElementForRowRotation()

		rtg.GenShare(testCtx.sk0, galEl, crp, share0)
		rtg.GenShare(testCtx.sk1, galEl, crp, share1)
		rtg.GenShare(testCtx.sk2, galEl, crp, share2)

		rtg.Aggregate(share0, share1, share0)
		rtg.Aggregate(share0, share2, share0)

		rotKeySet := rlwe.NewRotationKeySet(params, []uint64{galEl})
		rtg.GenRotationKey(share0, crp, rotKeySet.Keys[galEl])

		skIn := testCtx.skIdeal.CopyNew()
		skOut := testCtx.skIdeal.CopyNew()
		galElInv := ring.ModExp(galEl, uint64(2*params.N()-1), uint64(2*params.N()))
		ring.PermuteNTT(testCtx.skIdeal.Value.Q, galElInv, skOut.Value.Q)
		ring.PermuteNTT(testCtx.skIdeal.Value.P, galElInv, skOut.Value.P)

		swk := rotKeySet.Keys[galEl]

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

		// Worst bound of inner sum
		// N*#Keys*(N * #Parties * floor(sigma*6) + #Parties * floor(sigma*6) + N * #Parties  +  #Parties * floor(6*sigma))
		log2Bound := bits.Len64(3 * uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringQ.Modulus)-1, ringQ, swk.Value[0][0].Q))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringP.Modulus)-1, ringP, swk.Value[0][0].P))
	})
}

func testMarshalling(testCtx testContext, t *testing.T) {

	params := testCtx.params

	ciphertext := &rlwe.Ciphertext{Value: []*ring.Poly{params.RingQ().NewPoly(), params.RingQ().NewPoly()}}
	testCtx.uniformSampler.Read(ciphertext.Value[0])
	testCtx.uniformSampler.Read(ciphertext.Value[1])

	t.Run(testString(params, "Marshalling/CKG/"), func(t *testing.T) {
		ckg := NewCKGProtocol(testCtx.params)
		KeyGenShareBefore := ckg.AllocateShares()
		crs := ckg.SampleCRP(testCtx.crs)

		ckg.GenShare(testCtx.sk0, crs, KeyGenShareBefore)
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
		require.Equal(t, KeyGenShareBefore.Value.Q.Degree(), KeyGenShareAfter.Value.Q.Degree())
		require.Equal(t, KeyGenShareBefore.Value.P.Degree(), KeyGenShareAfter.Value.P.Degree())
		require.Equal(t, KeyGenShareBefore.Value.Q.LenModuli(), KeyGenShareAfter.Value.Q.LenModuli())
		require.Equal(t, KeyGenShareBefore.Value.P.LenModuli(), KeyGenShareAfter.Value.P.LenModuli())
		require.Equal(t, KeyGenShareAfter.Value.Q.Coeffs, KeyGenShareBefore.Value.Q.Coeffs)
		require.Equal(t, KeyGenShareAfter.Value.P.Coeffs, KeyGenShareBefore.Value.P.Coeffs)
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

		require.Equal(t, SwitchShare.Value[0].Degree(), SwitchShareReceiver.Value[0].Degree())
		require.Equal(t, SwitchShare.Value[1].Degree(), SwitchShareReceiver.Value[1].Degree())
		require.Equal(t, SwitchShare.Value[0].LenModuli(), SwitchShareReceiver.Value[0].LenModuli())
		require.Equal(t, SwitchShare.Value[1].LenModuli(), SwitchShareReceiver.Value[1].LenModuli())
		require.Equal(t, SwitchShare.Value[0].Coeffs, SwitchShareReceiver.Value[0].Coeffs)
		require.Equal(t, SwitchShare.Value[1].Coeffs, SwitchShareReceiver.Value[1].Coeffs)
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

		crp := RKGProtocol.SampleCRP(testCtx.crs)

		RKGProtocol.GenShareRoundOne(testCtx.sk0, crp, ephSk0, share10)

		data, err := share10.MarshalBinary()
		require.NoError(t, err)

		rkgShare := new(RKGShare)
		err = rkgShare.UnmarshalBinary(data)
		require.NoError(t, err)

		require.Equal(t, len(rkgShare.Value), len(share10.Value))
		for i, val := range share10.Value {
			require.Equal(t, len(rkgShare.Value[i][0].Q.Coeffs), len(val[0].Q.Coeffs))
			require.Equal(t, len(rkgShare.Value[i][0].P.Coeffs), len(val[0].P.Coeffs))
			require.Equal(t, rkgShare.Value[i][0].Q.Coeffs, val[0].Q.Coeffs)
			require.Equal(t, rkgShare.Value[i][0].P.Coeffs, val[0].P.Coeffs)

			require.Equal(t, len(rkgShare.Value[i][1].Q.Coeffs), len(val[1].Q.Coeffs))
			require.Equal(t, len(rkgShare.Value[i][1].P.Coeffs), len(val[1].P.Coeffs))
			require.Equal(t, rkgShare.Value[i][1].Q.Coeffs, val[1].Q.Coeffs)
			require.Equal(t, rkgShare.Value[i][1].P.Coeffs, val[1].P.Coeffs)
		}
	})

	t.Run(testString(params, "Marshalling/RTG/"), func(t *testing.T) {

		//check RTGShare

		galEl := testCtx.params.GaloisElementForColumnRotationBy(64)

		rtg := NewRTGProtocol(testCtx.params)
		rtgShare := rtg.AllocateShares()

		crp := rtg.SampleCRP(testCtx.crs)

		rtg.GenShare(testCtx.sk1, galEl, crp, rtgShare)

		data, err := rtgShare.MarshalBinary()
		require.NoError(t, err)

		resRTGShare := new(RTGShare)
		err = resRTGShare.UnmarshalBinary(data)
		require.NoError(t, err)

		require.Equal(t, len(resRTGShare.Value), len(rtgShare.Value))

		for i, val := range rtgShare.Value {
			require.Equal(t, len(resRTGShare.Value[i].Q.Coeffs), len(val.Q.Coeffs))
			require.Equal(t, resRTGShare.Value[i].Q.Coeffs, val.Q.Coeffs)
			require.Equal(t, len(resRTGShare.Value[i].P.Coeffs), len(val.P.Coeffs))
			require.Equal(t, resRTGShare.Value[i].P.Coeffs, val.P.Coeffs)
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
