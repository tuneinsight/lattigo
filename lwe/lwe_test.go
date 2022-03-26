package lwe

import (
	"encoding/json"
	"flag"
	"fmt"
	"github.com/stretchr/testify/assert"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"math"
	"runtime"
	"testing"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

var TestParams = []rlwe.ParametersLiteral{rlwe.TestPN12QP109, rlwe.TestPN13QP218}

func testString(params rlwe.Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

func TestLWE(t *testing.T) {
	defaultParams := TestParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = TestParams[:1] // the short test suite runs for ring degree N=2^12, 2^13
	}

	if *flagParamString != "" {
		var jsonParams rlwe.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []rlwe.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams[1:] {

		params, err := rlwe.NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		for _, testSet := range []func(params rlwe.Parameters, t *testing.T){
			testBinFHE,
			testLUT,
			testRLWEToLWE,
			testLWEToRLWE,
		} {
			testSet(params, t)
			runtime.GC()
		}
	}
}

func sign(x float64) (y float64) {
	if x < 0 {
		return -1
	}

	if x == 0 {
		return 0
	}

	return 1
}
func testBinFHE(params rlwe.Parameters, t *testing.T) {
	var err error

	// N=1024, Q=0x7fff801 -> 2^131
	paramsLUT, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:     10,
		Q:        []uint64{0x7fff801},
		P:        []uint64{},
		Sigma:    rlwe.DefaultSigma,
		LogBase2: 7,
	})

	assert.Nil(t, err)

	// N=512, Q=0x3001 -> 2^135
	paramsLWE, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:  9,
		Q:     []uint64{0x3001},
		P:     []uint64{},
		Sigma: rlwe.DefaultSigma,
	})

	assert.Nil(t, err)

	paramsKS, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:     9,
		Q:        []uint64{0x7fff801},
		P:        []uint64{},
		Sigma:    rlwe.DefaultSigma,
		LogBase2: 0,
	})

	assert.Nil(t, err)

	t.Run(testString(paramsLUT, "BinFHE/"), func(t *testing.T) {

		ringLWE := paramsLWE.RingQ()
		ringLUT := paramsLUT.RingQ()

		scaleLWE := float64(paramsLWE.Q()[0]) / 4.0
		scaleLUT := float64(paramsLUT.Q()[0]) / 4.0

		slots := 1

		LUTPoly := InitGate(xorGate, ringLUT)

		lutPolyMap := make(map[int]*ring.Poly)
		for i := 0; i < slots; i++ {
			lutPolyMap[i] = LUTPoly
		}

		skLWE := rlwe.NewKeyGenerator(paramsLWE).GenSecretKey()
		encryptorLWE := rlwe.NewEncryptor(paramsLWE, skLWE)

		m0 := rlwe.NewPlaintext(paramsLWE, paramsLWE.MaxLevel())
		m0.Value.Coeffs[0][0] = uint64(1 * scaleLWE / 4.0)

		m1 := rlwe.NewPlaintext(paramsLWE, paramsLWE.MaxLevel())
		m1.Value.Coeffs[0][0] = uint64(1 * scaleLWE / 4.0)

		ctm0 := rlwe.NewCiphertextNTT(paramsLWE, 1, paramsLWE.MaxLevel())
		encryptorLWE.Encrypt(m0, ctm0)

		ctm1 := rlwe.NewCiphertextNTT(paramsLWE, 1, paramsLWE.MaxLevel())
		encryptorLWE.Encrypt(m1, ctm1)

		handler := NewHandler(paramsLUT, paramsLWE, nil)

		kgenLUT := rlwe.NewKeyGenerator(paramsLUT)
		skLUT := kgenLUT.GenSecretKey()
		LUTKEY := handler.GenLUTKey(skLUT, skLWE)

		skLWELarge := rlwe.NewSecretKey(paramsLWE)
		skLWELarge2 := rlwe.NewSecretKey(paramsLUT)
		ringLWE.InvNTT(skLWE.Value.Q, skLWELarge.Value.Q)
		ringLWE.InvMForm(skLWELarge.Value.Q, skLWELarge.Value.Q)

		for i := range skLWELarge.Value.Q.Coeffs[0] {
			c := skLWELarge.Value.Q.Coeffs[0][i]
			if c == paramsLWE.Q()[0]-1 {
				skLWELarge.Value.Q.Coeffs[0][i] = paramsLUT.Q()[0] - 1
				skLWELarge2.Value.Q.Coeffs[0][i*2] = paramsLUT.Q()[0] - 1
			} else if c != 0 {
				skLWELarge.Value.Q.Coeffs[0][i] = 1
				skLWELarge2.Value.Q.Coeffs[0][i*2] = 1
			}
		}

		paramsKS.RingQ().NTT(skLWELarge.Value.Q, skLWELarge.Value.Q)
		paramsKS.RingQ().MForm(skLWELarge.Value.Q, skLWELarge.Value.Q)

		paramsLUT.RingQ().NTT(skLWELarge2.Value.Q, skLWELarge2.Value.Q)
		paramsLUT.RingQ().MForm(skLWELarge2.Value.Q, skLWELarge2.Value.Q)

		skLUT2skLWE := kgenLUT.GenSwitchingKey(skLUT, skLWELarge)

		ringLWE.Add(ctm0.Value[0], ctm1.Value[0], ctm0.Value[0])
		ringLWE.Add(ctm0.Value[1], ctm1.Value[1], ctm0.Value[1])

		ctsLUT := handler.ExtractAndEvaluateLUT(ctm0, lutPolyMap, LUTKEY)

		tmp := rlwe.NewCiphertextNTT(paramsLUT, 1, paramsLUT.MaxLevel())
		handler.evalRLWE.SwitchKeysInPlace(0, ctsLUT[0].Value[1], skLUT2skLWE, handler.evalRLWE.Pool[1].Q, handler.evalRLWE.Pool[2].Q)
		ringLUT.AddLvl(0, ctsLUT[0].Value[0], handler.evalRLWE.Pool[1].Q, tmp.Value[0])
		ring.CopyValuesLvl(0, handler.evalRLWE.Pool[2].Q, tmp.Value[1])
		ctLWE := rlwe.NewCiphertextNTT(paramsKS, 1, paramsKS.MaxLevel())
		rlwe.SwitchCiphertextRingDegreeNTT(tmp, paramsKS.RingQ(), paramsLUT.RingQ(), ctLWE)

		for i := range ctLWE.Value {
			paramsKS.RingQ().InvNTT(ctLWE.Value[i], ctLWE.Value[i])

			Q := paramsKS.Q()[0]
			q := paramsLWE.Q()[0]
			ratio := float64(q) / float64(Q)

			for j := 0; j < paramsLWE.N(); j++ {

				c := ctLWE.Value[i].Coeffs[0][j]
				c = uint64(float64(c)*ratio + 0.5)
				ctLWE.Value[i].Coeffs[0][j] = c
			}

			paramsLWE.RingQ().NTT(ctLWE.Value[i], ctLWE.Value[i])
		}

		q := paramsLWE.Q()[0]
		qHalf := q >> 1

		decryptorLWE := rlwe.NewDecryptor(paramsLWE, skLWE)
		ptLWE := rlwe.NewPlaintext(paramsLWE, paramsLWE.MaxLevel())

		decryptorLWE.Decrypt(ctLWE, ptLWE)

		_ = scaleLUT

		c := ptLWE.Value.Coeffs[0][0]

		var a float64
		if c >= qHalf {
			a = -float64(q-c) / scaleLWE
		} else {
			a = float64(c) / scaleLWE
		}

		fmt.Println(math.Round(a * 4))
	})
}

func testLUT(params rlwe.Parameters, t *testing.T) {
	var err error

	// N=1024, Q=0x7fff801 -> 2^131
	paramsLUT, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:     8,
		Q:        []uint64{0x7fff801},
		P:        []uint64{},
		Sigma:    rlwe.DefaultSigma,
		LogBase2: 7,
	})

	assert.Nil(t, err)

	// N=512, Q=0x3001 -> 2^135
	paramsLWE, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:  7,
		Q:     []uint64{0x3001},
		P:     []uint64{},
		Sigma: rlwe.DefaultSigma,
	})

	assert.Nil(t, err)

	t.Run(testString(paramsLUT, "LUT/"), func(t *testing.T) {

		scaleLWE := float64(paramsLWE.Q()[0]) / 4.0
		scaleLUT := float64(paramsLUT.Q()[0]) / 4.0

		slots := 16

		LUTPoly := InitLUT(nandGate, scaleLUT, paramsLUT.RingQ(), -1, 1)

		lutPolyMap := make(map[int]*ring.Poly)
		for i := 0; i < slots; i++ {
			lutPolyMap[i] = LUTPoly
		}

		skLWE := rlwe.NewKeyGenerator(paramsLWE).GenSecretKey()
		encryptorLWE := rlwe.NewEncryptor(paramsLWE, skLWE)

		values := make([]float64, slots)
		for i := 0; i < slots; i++ {
			values[i] = -1 + float64(2*i)/float64(slots)
		}

		ptLWE := rlwe.NewPlaintext(paramsLWE, paramsLWE.MaxLevel())
		for i := range values {
			if values[i] < 0 {
				ptLWE.Value.Coeffs[0][i] = paramsLWE.Q()[0] - uint64(-values[i]*scaleLWE)
			} else {
				ptLWE.Value.Coeffs[0][i] = uint64(values[i] * scaleLWE)
			}
		}
		ctLWE := rlwe.NewCiphertextNTT(paramsLWE, 1, paramsLWE.MaxLevel())
		encryptorLWE.Encrypt(ptLWE, ctLWE)

		handler := NewHandler(paramsLUT, paramsLWE, nil)

		skLUT := rlwe.NewKeyGenerator(paramsLUT).GenSecretKey()
		LUTKEY := handler.GenLUTKey(skLUT, skLWE)

		ctsLUT := handler.ExtractAndEvaluateLUT(ctLWE, lutPolyMap, LUTKEY)

		q := paramsLUT.Q()[0]
		qHalf := q >> 1
		decryptorLUT := rlwe.NewDecryptor(paramsLUT, skLUT)
		ptLUT := rlwe.NewPlaintext(paramsLUT, paramsLUT.MaxLevel())
		for i := 0; i < slots; i++ {

			decryptorLUT.Decrypt(ctsLUT[i], ptLUT)

			c := ptLUT.Value.Coeffs[0][i]

			var a float64
			if c >= qHalf {
				a = -float64(q-c) / scaleLUT
			} else {
				a = float64(c) / scaleLUT
			}

			fmt.Printf("%7.4f - %7.4f - %7.4f\n", math.Round(a*32)/32, math.Round(a*8)/8, values[i])
		}
	})
}

func testRLWEToLWE(params rlwe.Parameters, t *testing.T) {
	t.Run(testString(params, "RLWEToLWE/"), func(t *testing.T) {
		kgen := rlwe.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encryptor := rlwe.NewEncryptor(params, sk)
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		ct := rlwe.NewCiphertextNTT(params, 1, params.MaxLevel())
		encryptor.Encrypt(pt, ct)

		skInvNTT := params.RingQ().NewPoly()
		params.RingQ().InvNTT(sk.Value.Q, skInvNTT)

		slotIndex := make(map[int]bool)
		for i := 0; i < params.N(); i++ {
			slotIndex[i] = true
		}

		LWE := RLWEToLWE(ct, params.RingQ(), slotIndex)

		for i := 0; i < params.RingQ().N; i++ {
			if math.Abs(DecryptLWE(LWE[i], params.RingQ(), skInvNTT)) > 19 {
				t.Error()
			}
		}
	})
}

func testLWEToRLWE(params rlwe.Parameters, t *testing.T) {
	t.Run(testString(params, "LWEToRLWE/"), func(t *testing.T) {
		kgen := rlwe.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encryptor := rlwe.NewEncryptor(params, sk)
		decryptor := rlwe.NewDecryptor(params, sk)
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		ct := rlwe.NewCiphertextNTT(params, 1, params.MaxLevel())
		encryptor.Encrypt(pt, ct)

		skInvNTT := params.RingQ().NewPoly()
		params.RingQ().InvNTT(sk.Value.Q, skInvNTT)

		slotIndex := make(map[int]bool)
		for i := 0; i < params.N(); i++ {
			slotIndex[i] = true
		}

		ctLWE := RLWEToLWE(ct, params.RingQ(), slotIndex)

		DecryptLWE(ctLWE[0], params.RingQ(), skInvNTT)

		handler := NewHandler(params, params, nil)

		ctRLWE := handler.LWEToRLWE(ctLWE)

		for i := 0; i < len(ctRLWE); i++ {
			decryptor.Decrypt(ctRLWE[i], pt)

			for j := 0; j < pt.Level()+1; j++ {

				c := pt.Value.Coeffs[j][0]

				if c >= params.RingQ().Modulus[j]>>1 {
					c = params.RingQ().Modulus[j] - c
				}

				if c > 19 {
					t.Fatal(i, j, c)
				}
			}
		}
	})
}
