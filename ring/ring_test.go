package ring

import (
	"bufio"
	"fmt"
	"log"
	"math/bits"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"testing"
	"time"
)

var folder = "test_data/"

// Name of the test vectors files

var files_60 = []string{
	"test_pol_60____8_2",
	"test_pol_60___16_2",
	"test_pol_60___32_2",
	"test_pol_60___64_2",
	"test_pol_60__128_2",
	"test_pol_60__256_2",
	"test_pol_60__512_2",
}

// Name of the test vectors files

var filesNTT_60 = []string{
	"test_pol_NTT_60____8_2",
	"test_pol_NTT_60___16_2",
	"test_pol_NTT_60___32_2",
	"test_pol_NTT_60___64_2",
	"test_pol_NTT_60__128_2",
	"test_pol_NTT_60__256_2",
	"test_pol_NTT_60__512_2",
}

func Test_Polynomial(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	for i := uint64(0); i < 1; i++ {

		N := uint64(2 << (12 + i))
		T := uint64(65537)

		Qi := Qi60[uint64(len(Qi60))-2<<i:]
		Pi := Pi60[uint64(len(Pi60))-((2<<i)+1):]

		sigma := 3.19

		contextT := NewContext()
		contextT.SetParameters(N, []uint64{T})
		contextT.GenNTTParams()

		contextQ := NewContext()
		contextQ.SetParameters(N, Qi)
		contextQ.GenNTTParams()

		contextP := NewContext()
		contextP.SetParameters(N, Pi)
		contextP.GenNTTParams()

		contextQP := NewContext()
		contextQP.Merge(contextQ, contextP)

		// ok!
		test_GenerateNTTPrimes(N, Qi[0], t)

		// ok!
		test_ImportExportPolyString(contextQ, t)

		// ok!
		test_Marshaler(contextQ, t)

		// TODO : check that the coefficients are within the bound
		test_GaussianPoly(sigma, contextQ, t)

		// ok!
		test_BRed(contextQ, t)

		// ok!
		test_MRed(contextQ, t)

		// ok!
		test_Shift(contextQ, t)

		// ok!
		test_GaloisShift(contextQ, t)

		// ok!
		test_MForm(contextQ, t)

		// ok!
		test_MulPoly(contextQ, t)

		// ok!
		test_MulPoly_Montgomery(contextQ, t)

		// ok!
		test_ExtendBasis(contextQ, contextP, contextQP, t)

		// ok!
		test_SimpleScaling(T, contextQ, contextP, t)

		// ok!
		test_ComplexScaling(T, contextQ, contextP, contextQP, t)

		test_MultByMonomial(contextQ, t)
	}
}

// Parses a file and return a slice whose elements are each line of the file
func getParamsFromString(filename string) []string {

	file, _ := os.Open(filename)

	defer file.Close()

	var lines []string

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}

	return lines
}

// creates a context from a slice of string of N and Modulies
func constructContextFromString(vs []string) *Context {

	// Get N
	tmpN, _ := strconv.ParseUint(vs[0], 0, 32)
	N := uint64(tmpN)

	// Get Qi
	ModulusString := strings.Split(strings.TrimSpace(vs[1]), " ")

	Modulus := make([]uint64, len(ModulusString))
	for i := range ModulusString {
		tmp, _ := strconv.ParseInt(ModulusString[i], 0, 64)
		Modulus[i] = uint64(tmp)
	}

	// Generate the context from N and Qi
	context := NewContext()
	context.SetParameters(N, Modulus)
	context.GenNTTParams()

	return context
}

// creates a ring from a slice of string of coefficients
func constructPolynomialsFromString(vs []string, context *Context) *Poly {

	// Extracts the coefficients in CRT and in NTT
	coeffs := make([][]uint64, len(context.Modulus))
	for i := range context.Modulus {
		coeffsString := strings.Split(strings.TrimSpace(vs[i+2]), " ")
		coeffs[i] = make([]uint64, context.N)

		for j := uint64(0); j < context.N; j++ {
			tmp, _ := strconv.ParseInt(coeffsString[j], 0, 64)
			coeffs[i][j] = uint64(tmp)
		}
	}

	//Create the polynomials and assigns the CRT and NTT coefficients
	pol := context.NewPoly()
	pol.SetCoefficients(coeffs)

	return pol
}

func test_Vectors_NTT(t *testing.T) {

	for x := uint64(0); x < 7; x++ {

		vs := getParamsFromString(fmt.Sprintf(folder + files_60[x]))
		vs_ntt := getParamsFromString(fmt.Sprintf(folder + filesNTT_60[x]))

		context := constructContextFromString(vs)

		t.Run(fmt.Sprintf("Test_Vectors_NTT/N=%d/limbs=%d", context.N, len(context.Modulus)), func(t *testing.T) {

			Polx := constructPolynomialsFromString(vs, context)

			PolNTT := constructPolynomialsFromString(vs_ntt, context)

			CRTCoeffs := Polx.GetCoefficients()

			context.NTT(Polx, Polx)

			for i := range context.Modulus {
				for j := range Polx.Coeffs {
					if Polx.Coeffs[i][j] != PolNTT.Coeffs[i][j] {
						t.Errorf("error : NTT coeffs file n°%v: want %v, have %v", x, PolNTT.Coeffs[i][j], Polx.Coeffs[i][j])
						continue
					}
				}
			}

			context.InvNTT(Polx, Polx)

			for i := range context.Modulus {
				for j := range Polx.Coeffs {
					if Polx.Coeffs[i][j] != CRTCoeffs[i][j] {
						t.Errorf("error : InvNTT coeffs file n°%v: want %v, have %v", x, CRTCoeffs[i][j], Polx.Coeffs[i][j])
						continue
					}

				}
			}
		})
	}
}

func test_GenerateNTTPrimes(N, Qi uint64, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/GenerateNTTPrimes", N), func(t *testing.T) {

		primes, err := GenerateNTTPrimes(N, Qi, 100, 60, true)

		if err != nil {
			t.Errorf("error : generateNTTPrimes")
		}

		for _, q := range primes {
			if bits.Len64(q) != 60 || q&((N<<1)-1) != 1 || IsPrime(q) != true {
				t.Errorf("error : GenerateNTTPrimes for q = %v", q)
				break
			}
		}
	})
}

func test_ImportExportPolyString(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/ImportExportPolyString", context.N, len(context.Modulus)), func(t *testing.T) {

		p0 := context.NewUniformPoly()
		p1 := context.NewPoly()

		context.SetCoefficientsString(context.PolyToString(p0), p1)

		if context.Equal(p0, p1) != true {
			t.Errorf("error : import/export ring from/to string")
		}
	})
}

func test_Marshaler(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MarshalContext", context.N, len(context.Modulus)), func(t *testing.T) {

		data, _ := context.MarshalBinary()

		contextTest := NewContext()
		contextTest.UnMarshalBinary(data)

		if contextTest.N != context.N {
			t.Errorf("ERROR encoding/decoding N")
		}

		for i := range context.Modulus {
			if contextTest.Modulus[i] != context.Modulus[i] {
				t.Errorf("ERROR encoding/decoding Modulus")
			}
		}
	})

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MarshalPoly", context.N, len(context.Modulus)), func(t *testing.T) {

		p := context.NewUniformPoly()
		pTest := context.NewPoly()

		data, _ := p.MarshalBinary()

		pTest, _ = pTest.UnMarshalBinary(data)

		for i := range context.Modulus {
			for j := uint64(0); j < context.N; j++ {
				if p.Coeffs[i][j] != pTest.Coeffs[i][j] {
					t.Errorf("PolyBytes Import Error : want %v, have %v", p.Coeffs[i][j], pTest.Coeffs[i][j])
				}
			}
		}
	})
}

func test_GaussianPoly(sigma float64, context *Context, t *testing.T) {

	bound := int(sigma * 6)
	KYS := context.NewKYSampler(sigma, bound)
	TS := context.NewTernarySampler()

	pol := context.NewPoly()

	t.Run(fmt.Sprintf("N=%d/limbs=%d/NewGaussPoly", context.N, len(context.Modulus)), func(t *testing.T) {
		KYS.Sample(pol)
	})

	countOne := 0
	countZer := 0
	countMOn := 0
	t.Run(fmt.Sprintf("N=%d/limbs=%d/NewTernaryPoly", context.N, len(context.Modulus)), func(t *testing.T) {
		if err := TS.Sample(1.0/3, pol); err != nil {
			log.Fatal(err)
		}

		//fmt.Println(pol.Coeffs[0])

		for i := range pol.Coeffs[0] {
			if pol.Coeffs[0][i] == context.Modulus[0]-1 {
				countMOn += 1
			}

			if pol.Coeffs[0][i] == 0 {
				countZer += 1
			}

			if pol.Coeffs[0][i] == 1 {
				countOne += 1
			}
		}

		//fmt.Println("-1 :", countMOn)
		//fmt.Println(" 0 :", countZer)
		//fmt.Println(" 1 :", countOne)
	})
}

func test_BRed(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/BRed", context.N, len(context.Modulus)), func(t *testing.T) {
		for j, q := range context.Modulus {

			bigQ := NewUint(q)

			for i := 0; i < 65536; i++ {
				x := rand.Uint64() % q
				y := rand.Uint64() % q

				result := NewUint(x)
				result.Mul(result, NewUint(y))
				result.Mod(result, bigQ)

				test := BRed(x, y, q, context.bredParams[j])
				want := result.Uint64()

				if test != want {
					t.Errorf("error : 128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, want)
					break
				}
			}
		}
	})
}

func test_MRed(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MRed", context.N, len(context.Modulus)), func(t *testing.T) {
		for j := range context.Modulus {

			q := context.Modulus[j]

			bigQ := NewUint(q)

			for i := 0; i < 65536; i++ {

				x := rand.Uint64() % q
				y := rand.Uint64() % q

				result := NewUint(x)
				result.Mul(result, NewUint(y))
				result.Mod(result, bigQ)

				test := MRed(x, MForm(y, q, context.bredParams[j]), q, context.mredParams[j])
				want := result.Uint64()

				if test != want {
					t.Errorf("error : 128bit montgomery multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, want)
					break
				}

			}
		}
	})
}

func test_Shift(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/Shift", context.N, len(context.Modulus)), func(t *testing.T) {
		pWant := context.NewUniformPoly()
		pTest := context.NewPoly()

		context.Shift(pTest, 1, pWant)
		for i := range context.Modulus {
			for j := uint64(0); j < context.N; j++ {
				if pTest.Coeffs[i][(j+1)%context.N] != pWant.Coeffs[i][j] {
					t.Errorf("error : ShiftR")
					break
				}
			}
		}
	})
}

func test_GaloisShift(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/GaloisShift", context.N, len(context.Modulus)), func(t *testing.T) {

		pWant := context.NewUniformPoly()
		pTest := pWant.CopyNew()

		context.BitReverse(pTest, pTest)
		context.InvNTT(pTest, pTest)
		context.Rotate(pTest, 1, pTest)
		context.NTT(pTest, pTest)
		context.BitReverse(pTest, pTest)
		context.Reduce(pTest, pTest)

		context.Shift(pWant, 1, pWant)

		for i := range context.Modulus {
			for j := uint64(0); j < context.N; j++ {
				if pTest.Coeffs[i][j] != pWant.Coeffs[i][j] {
					t.Errorf("error : GaloisShiftR")
					break
				}
			}
		}
	})
}

func test_MForm(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MForm", context.N, len(context.Modulus)), func(t *testing.T) {

		polWant := context.NewUniformPoly()
		polTest := context.NewPoly()

		context.MForm(polWant, polTest)
		context.InvMForm(polTest, polTest)

		if context.Equal(polWant, polTest) != true {
			t.Errorf("error : 128bit MForm/InvMForm")
		}
	})

}

func test_MulPoly(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MulPoly", context.N, len(context.Modulus)), func(t *testing.T) {

		p1 := context.NewUniformPoly()
		p2 := context.NewUniformPoly()
		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.Reduce(p1, p1)
		context.Reduce(p2, p2)

		context.MulPolyNaive(p1, p2, p3Want)

		context.MulPoly(p1, p2, p3Test)

		for i := uint64(0); i < context.N; i++ {
			if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
				t.Errorf("ERROR MUL COEFF %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
				break
			}
		}
	})
}

func test_MulPoly_Montgomery(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MulPoly_Montgomery", context.N, len(context.Modulus)), func(t *testing.T) {
		p1 := context.NewUniformPoly()
		p2 := context.NewUniformPoly()
		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.MForm(p1, p1)
		context.MForm(p2, p2)

		context.MulPolyNaiveMontgomery(p1, p2, p3Want)
		context.MulPolyMontgomery(p1, p2, p3Test)

		context.InvMForm(p3Test, p3Test)
		context.InvMForm(p3Want, p3Want)

		for i := uint64(0); i < context.N; i++ {
			if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
				t.Errorf("ERROR MUL COEFF %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
				break
			}
		}
	})
}

func test_ExtendBasis(contextQ, contextP, contextQP *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d+%d/ExtendBasis", contextQ.N, len(contextQ.Modulus), len(contextP.Modulus)), func(t *testing.T) {

		basisextender := NewBasisExtender(contextQ, contextP)

		coeffs := make([]*Int, contextQ.N)
		for i := uint64(0); i < contextQ.N; i++ {
			coeffs[i] = RandInt(contextQ.ModulusBigint)
		}

		PolTest := contextQ.NewPoly()
		PolWant := contextQP.NewPoly()

		contextQ.SetCoefficientsBigint(coeffs, PolTest)
		contextQP.SetCoefficientsBigint(coeffs, PolWant)

		basisextender.ExtendBasis(PolTest, PolTest)

		for i := range contextQP.Modulus {
			for j := uint64(0); j < contextQ.N; j++ {
				if PolTest.Coeffs[i][j] != PolWant.Coeffs[i][j] {
					t.Errorf("error extendBasis, want %v - has %v", PolTest.Coeffs[i][j], PolWant.Coeffs[i][j])
					break
				}
			}
		}
	})
}

func test_SimpleScaling(T uint64, contextT, contextQ *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/T=%v/SimpleScaling", contextQ.N, len(contextQ.Modulus), T), func(t *testing.T) {

		rescaler := NewSimpleScaler(T, contextQ)

		coeffs := make([]*Int, contextQ.N)
		for i := uint64(0); i < contextQ.N; i++ {
			coeffs[i] = RandInt(contextQ.ModulusBigint)
		}

		coeffsWant := make([]*Int, contextQ.N)
		for i := range coeffs {
			coeffsWant[i] = Copy(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], NewUint(T))
			coeffsWant[i].DivRound(coeffsWant[i], contextQ.ModulusBigint)
			coeffsWant[i].Mod(coeffsWant[i], NewUint(T))
		}

		PolTest := contextQ.NewPoly()
		PolWant := contextT.NewPoly()

		contextQ.SetCoefficientsBigint(coeffs, PolTest)

		rescaler.Scale(PolTest, PolTest)

		contextT.SetCoefficientsBigint(coeffsWant, PolWant)

		for i := uint64(0); i < contextT.N; i++ {
			if PolWant.Coeffs[0][i] != PolTest.Coeffs[0][i] {
				t.Errorf("error : simple scaling, want %v have %v", PolWant.Coeffs[0][i], PolTest.Coeffs[0][i])
				break
			}
		}
	})
}

func test_ComplexScaling(T uint64, contextQ, contextP, contextQP *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d+%d/T=%d/ComplexScaling", contextQ.N, len(contextQ.Modulus), len(contextP.Modulus), T), func(t *testing.T) {

		complexRescaler := NewComplexScaler(T, contextQ, contextP)

		coeffs := make([]*Int, contextQ.N)
		for i := uint64(0); i < contextQ.N; i++ {
			coeffs[i] = RandInt(contextQP.ModulusBigint)
			coeffs[i].Div(coeffs[i], NewUint(10))
		}

		//coeffs[0].SetString("323702478295050366752968655518337494486119222600846780736994466085672189264880022173262")
		//coeffs[1].SetString("4888963624565002661780158142973932724425930851432042682757127743238542390594193351997723")
		//coeffs[2].SetString("4282984724308796735476641110474688971631286863366811526364125630744789177542625267835653")
		//coeffs[3].SetString("3443920543591757808834812628691980702836925006139799151579776655221448318527620147901056")
		//coeffs[4].SetString("2435992701314338051709289453765710565896473610780678347020688344988696126710239045561833")
		//coeffs[5].SetString("3577369130953681317513231564259583867075569841269600848965005473743746265741712227949357")
		//coeffs[6].SetString("3650829186364460588505284366902239740496861722443649340705308513232282111686882618728671")
		//coeffs[7].SetString("5413725621145553596020326503684835624178908809180347713780542414099674315138466376924388")
		//coeffs[8].SetString("2921703598815986183949445434906117845024974748774925667243895094611669622018594882674160")

		coeffsWant := make([]*Int, contextQ.N)
		for i := range coeffs {
			coeffsWant[i] = Copy(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], NewUint(T))
			coeffsWant[i].DivRound(coeffsWant[i], contextQ.ModulusBigint)
		}

		PolTest := contextQP.NewPoly()
		PolWant := contextQ.NewPoly()

		contextQP.SetCoefficientsBigint(coeffs, PolTest)

		complexRescaler.Scale(PolTest, PolTest)
		contextQ.SetCoefficientsBigint(coeffsWant, PolWant)

		for i := uint64(0); i < contextQ.N; i++ {
			for j := 0; j < len(contextQ.Modulus); j++ {
				if PolWant.Coeffs[j][i] != PolTest.Coeffs[j][i] {
					t.Errorf("error : complex scaling coeff %v Qi%v = %s, want %v have %v", i, j, coeffs[i].String(), PolWant.Coeffs[j][i], PolTest.Coeffs[j][i])
					break
				}
			}
		}
	})
}

func test_MultByMonomial(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MultByMonomial", context.N, len(context.Modulus)), func(t *testing.T) {

		p1 := context.NewUniformPoly()

		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.MultByMonomial(p1, 1, p3Test)
		context.MultByMonomial(p3Test, 8, p3Test)

		context.MultByMonomial(p1, 9, p3Want)

		for i := uint64(0); i < context.N; i++ {
			if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
				t.Errorf("ERROR MUL BY MONOMIAL %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
				break
			}
		}
	})
}
