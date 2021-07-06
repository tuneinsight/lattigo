package main

import(
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/utils"
)

func main(){
	var err error 
	ParametersLiteral := ckks.ParametersLiteral{
		LogN:     13,
		LogSlots: 12,
		Scale:    1 << 50,
		Sigma:    rlwe.DefaultSigma,
		Q: []uint64{
			0x10000000006e0001, // 60 Q0
			0x100000000060001, // 58 CtS
			0xfffffffff00001,  // 58 CtS
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
		},
	}

	SlotsToCoeffsParameters := ckks.SlotsToCoeffsParameters{
		LevelStart:  2,
		BSGSRatio:   16.0,
		BitReversed: false,
		ScalingFactor: [][]float64{
			{0x100000000060001},
			{0xfffffffff00001},
		},
	}

	H := 192

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil{
		panic(err)
	}

	ringQ := params.RingQ()

	encoder := ckks.NewEncoder(params)
	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeySparse(H)
	encryptor := ckks.NewEncryptor(params, sk)
	decryptor := ckks.NewDecryptor(params, sk)

	SlotsToCoeffsMatrix := SlotsToCoeffsParameters.GenSlotsToCoeffsMatrix(&params, params.LogSlots(), 1.0, encoder)

	// Gets the rotations indexes for SlotsToCoeffs
	rotations := SlotsToCoeffsParameters.RotationsForSlotsToCoeffs(&params, params.LogSlots())

	rotKey := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)

	eval := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotKey})

	values := make([]complex128, params.Slots())
	for i := range values {
		values[i] = complex(float64(i+1)/float64(params.Slots()), 1+float64(i+1)/float64(params.Slots()))
	}

	plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
	// Must encode with 2*Slots because a real vector is returned
	encoder.Encode(plaintext, values, utils.MinInt(params.LogSlots()+1, params.LogN()-1)) 
	ct := encryptor.EncryptNew(plaintext)

	// Homomorphic Decoding
	ct = ckks.SlotsToCoeffs(ct, nil, SlotsToCoeffsMatrix, eval)
 	
 	// Decrypt and print coefficient domain
	coeffsFloat := encoder.DecodeCoeffsPublic(decryptor.DecryptNew(ct), 0)
	valuesFloat := make([]complex128, params.Slots())
	gap := params.N() / (2 * params.Slots())
	for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {
		valuesFloat[i] = complex(coeffsFloat[idx], coeffsFloat[idx+(params.N()>>1)])
	}
	ckks.SliceBitReverseInPlaceComplex128(valuesFloat, params.Slots())
	for i := range valuesFloat{
		fmt.Println(i, valuesFloat[i])
	}

	// RLWE -> LWE Extraction
	lweReal, lweImag := ExtractLWESamplesBitReversed(ct, params)

	skInvNTT := sk.Value.CopyNew()
	ringQ.InvNTT(skInvNTT, skInvNTT)

	for i := 0; i < params.Slots(); i++{
		DecryptLWE(ringQ, lweReal[i], params.Scale(), skInvNTT)
		DecryptLWE(ringQ, lweImag[i], params.Scale(), skInvNTT)
	}
}

func DecryptLWE(ringQ *ring.Ring, lwe RNSLWESample, scale float64, skInvNTT *ring.Poly) (float64){

	ringQ.MulCoeffsMontgomeryLvl(0, lwe.a, skInvNTT, lwe.a)
	qi := ringQ.Modulus[0]
	tmp0 := lwe.a.Coeffs[0]
	tmp1 := lwe.b[0]
	for j := 0; j < ringQ.N; j++{
		tmp1 = ring.CRed(tmp1+tmp0[j], qi)
	}

	if tmp1 >= ringQ.Modulus[0]>>1{
		tmp1 = ringQ.Modulus[0]-tmp1
		return -float64(tmp1)/scale
	}else{
		return float64(tmp1)/scale
	}
}

type RNSLWESample struct{
	b []uint64
	a *ring.Poly
}

func ExtractLWESamplesBitReversed(ct *ckks.Ciphertext, params ckks.Parameters) (LWEReal, LWEImag []RNSLWESample){

	ringQ := params.RingQ()

	ringQ.InvNTTLvl(ct.Level(), ct.Value[0], ct.Value[0])
	ringQ.InvNTTLvl(ct.Level(), ct.Value[1], ct.Value[1])

	LWEReal = make([]RNSLWESample, params.Slots())
	LWEImag = make([]RNSLWESample, params.Slots())

	// Copies coefficients multiplied by X^{N-1} in reverse order : 
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	acc := ringQ.NewPolyLvl(ct.Level())
	for i, qi := range ringQ.Modulus[:ct.Level()+1]{
		tmp0 := acc.Coeffs[i]
		tmp1 := ct.Value[1].Coeffs[i]
		tmp0[0] = tmp1[0]
		for j := 1; j < ringQ.N; j++{
			tmp0[j] = qi-tmp1[ringQ.N-j]
		}
	}

	pol := ct.Value[0]

	gap := params.N() / (2 * params.Slots()) // Gap between plaintext coefficient if sparse packed

	// Real values
	for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {

		iRev := utils.BitReverse64(uint64(i), uint64(params.LogSlots())) 

		LWEReal[iRev].b = make([]uint64, len(pol.Coeffs))
		for j := range pol.Coeffs{
			LWEReal[iRev].b[j] = pol.Coeffs[j][idx]
		}

		LWEReal[iRev].a = acc.CopyNew()

		// Multiplies the accumulator by X^{N/(2*slots)}
		MulBySmallMonomial(ringQ, acc, gap)
	}

	// Imaginary values
	for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {

		iRev := utils.BitReverse64(uint64(i), uint64(params.LogSlots())) 

		LWEImag[iRev].b = make([]uint64, len(pol.Coeffs))
		for j := range pol.Coeffs{
			LWEImag[iRev].b[j] = pol.Coeffs[j][idx+(params.N()>>1)]
		}

		LWEImag[iRev].a = acc.CopyNew()

		// Multiplies the accumulator by X^{N/(2*slots)}
		MulBySmallMonomial(ringQ, acc, gap)
	}
	return
}

func MulBySmallMonomial(ringQ *ring.Ring, pol *ring.Poly, n int){
	for i, qi := range ringQ.Modulus[:pol.Level()+1]{
		pol.Coeffs[i] = append(pol.Coeffs[i][ringQ.N-n:], pol.Coeffs[i][:ringQ.N-n]...)
		tmp := pol.Coeffs[i]
		for j := 0; j < n; j++{
			tmp[j] = qi - tmp[j]
		}
	}
}