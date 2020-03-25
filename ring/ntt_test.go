package ring

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
	"testing"
)

var folder = "test_data/"

// Name of the test vectors files
var files60 = []string{
	"test_pol_60____8_2",
	"test_pol_60___16_2",
	"test_pol_60___32_2",
	"test_pol_60___64_2",
	"test_pol_60__128_2",
	"test_pol_60__256_2",
	"test_pol_60__512_2",
}

// Name of the test vectors files
var filesNTT60 = []string{
	"test_pol_NTT_60____8_2",
	"test_pol_NTT_60___16_2",
	"test_pol_NTT_60___32_2",
	"test_pol_NTT_60___64_2",
	"test_pol_NTT_60__128_2",
	"test_pol_NTT_60__256_2",
	"test_pol_NTT_60__512_2",
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

func Test_NTT(t *testing.T) {

	for x := uint64(0); x < 7; x++ {

		vs := getParamsFromString(fmt.Sprintf(folder + files60[x]))
		vsNtt := getParamsFromString(fmt.Sprintf(folder + filesNTT60[x]))

		context := constructContextFromString(vs)

		t.Run(fmt.Sprintf("N=%d/limbs=%d", context.N, len(context.Modulus)), func(t *testing.T) {

			Polx := constructPolynomialsFromString(vs, context)

			PolNTT := constructPolynomialsFromString(vsNtt, context)

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
