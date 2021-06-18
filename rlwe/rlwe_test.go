package rlwe

import (
	"math/big"
	"encoding/json"
	"flag"
	"math"
	"runtime"
	"testing"
	//"github.com/ldsec/lattigo/v2/ring"
	//"github.com/ldsec/lattigo/v2/utils"
	//"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func TestRLWE(t *testing.T) {
	defaultParams := TestParams[:1] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = TestParams[:1] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = append(TestParams) //, DefaultPostQuantumParams...) // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams {
		params, err := NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		kgen := NewKeyGenerator(params)

		for _, testSet := range []func(kgen KeyGenerator, t *testing.T){
			testGenKeyPair,
			testSwitchKeyGen,
		} {
			testSet(kgen, t)
			runtime.GC()
		}
	}
}

func testGenKeyPair(kgen KeyGenerator, t *testing.T){


	t.Run("PKGen", func(t *testing.T){
		// Checks that sum([-as + e, a] + [as])) < N * 6 * sigma
		params := kgen.(*keyGenerator).params

		ringQP := params.RingQP()

		sk, pk := kgen.GenKeyPair()

		// [-as + e] + [as]
		ringQP.MulCoeffsMontgomeryAndAdd(sk.Value, pk.Value[1], pk.Value[0])
		ringQP.InvNTT(pk.Value[0], pk.Value[0])
		coeffsBigint := make([]*big.Int, ringQP.N)
		for i := range coeffsBigint{
			coeffsBigint[i] = new(big.Int)
		}
		// Reconstruct the poly mod QP
		ringQP.PolyToBigintCenteredLvl(pk.Value[0].Level(), pk.Value[0], coeffsBigint)

		// Sums all coefficients
		sum := new(big.Int).SetInt64(0)
		for i := range coeffsBigint{
			sum.Add(sum, coeffsBigint[i].Abs(coeffsBigint[i]))
		}

		// Compares with the bound (6*sigma)*N
		bound := new(big.Int).SetInt64(int64(math.Floor(DefaultSigma*6)))
		bound.Mul(bound, new(big.Int).SetInt64(int64(params.N())))
		sign := bound.Cmp(sum)

		require.Greater(t, sign, -1)
	})
}

func testSwitchKeyGen(kgen KeyGenerator, t *testing.T){

	t.Run("SWKGen", func(t *testing.T){

		params := kgen.(*keyGenerator).params

		ringQ := params.RingQ()
		skIn := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		swk := NewSwitchingKey(params)
		kgen.(*keyGenerator).newSwitchingKey(skIn.Value, skOut.Value, swk)

		// Decrypts
		// [-asIn + w*P*sOut + e, a] + [asIn]
		for j := range swk.Value{
			ringQ.MulCoeffsMontgomeryAndAdd(swk.Value[j][1], skOut.Value, swk.Value[j][0])
		}

		poly := swk.Value[0][0]

		// Sums all basis together (equivalent to multiplying with CRT decomposition of 1)
		// sum([1]_w * [w*P*sOut + e]) = P*sOut + sum(e)
		for j := range swk.Value{
			if j > 0{
				ringQ.Add(poly, swk.Value[j][0], poly)
			}
		}

		// sOut * P 
		ringQ.MulScalarBigint(skIn.Value, kgen.(*keyGenerator).pBigInt, skIn.Value)

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(poly, skIn.Value, poly)

		// Checks that the error is below the bound
		ringQ.InvNTT(poly, poly)
		ringQ.InvMForm(poly, poly)

		coeffsBigint := make([]*big.Int, ringQ.N)
		for i := range coeffsBigint{
			coeffsBigint[i] = new(big.Int)
		}
		ringQ.PolyToBigintCenteredLvl(params.QCount()-1, poly, coeffsBigint)

		// Sums all coefficients
		sum := new(big.Int).SetInt64(0)
		for i := range coeffsBigint{
			sum.Add(sum, coeffsBigint[i].Abs(coeffsBigint[i]))
		}

		// Compares with the bound
		bound := new(big.Int).SetInt64(int64(math.Floor(DefaultSigma*6)))
		bound.Mul(bound, new(big.Int).SetInt64(int64(params.N() * len(swk.Value))))
		sign := bound.Cmp(sum)

		require.Greater(t, sign, -1)

	})
}
