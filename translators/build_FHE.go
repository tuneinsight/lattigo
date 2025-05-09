package main

import (
	"bufio"
	"flag"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func readFile(path string) (expected string, operations []string, err error) {
	file, err := os.Open(path)
	if err != nil {
		return "", nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := scanner.Text()

		switch line {
		case "# Expected":
			section = "expected"
			continue
		case "# Operations":
			section = "operations"
			continue
		case "":
			continue
		}

		switch section {
		case "expected":
			expected = line
		case "operations":
			operations = append(operations, line)
		}
	}

	return expected, operations, scanner.Err()
}

func (lattigo *LattigoFHE) parseOperation(line string) (int, *Term, string) {
	if line == "" || strings.HasPrefix(line, "#") {
		return -1, nil, ""
	}

	parts := strings.Split(line, " ")
	lineNum, _ := strconv.Atoi(strings.TrimSuffix(parts[0], ":"))
	op := parts[1]
	cs := parseIntArray(parts[2])
	isSecret, _ := strconv.ParseBool(parts[3])
	metadata := parts[4]

	term := &Term{
		op:       op,
		children: cs,
		secret:   isSecret,
	}
	if _, ok := lattigo.terms[lineNum]; !ok {
		lattigo.terms[lineNum] = term
	}
	return lineNum, term, metadata
}

type Term struct {
	op       string
	children []int
	secret   bool
}

type LattigoFHE struct {
	params *ckks.Parameters
	terms  map[int]*Term            // stores term info
	env    map[int]*rlwe.Ciphertext // stores ciphertexts
	ptEnv  map[int][]float64        // stores plaintexts
	n      int
	eval   *ckks.Evaluator
	enc    *rlwe.Encryptor
	ecd    *ckks.Encoder
	dec    *rlwe.Decryptor
	instructionsPath string
}

func NewLattigoFHE(n int, instructionsPath string) *LattigoFHE {
	return &LattigoFHE{
		terms: make(map[int]*Term),
		env:   make(map[int]*rlwe.Ciphertext),
		ptEnv: make(map[int][]float64),
		n:     n,
		instructionsPath: instructionsPath,
	}
}

func findUniqueRots(operations []string) []int {
	var rots []int
	for _, operation := range operations {
		if strings.Contains(operation, "ROT") {
			parts := strings.Split(operation, " ")
			if len(parts) > 4 {
				rot, _ := strconv.Atoi(parts[4])
				rots = append(rots, rot)
			}
		}
	}
	return rots
}

func (lattigo *LattigoFHE) createContext(depth int, rots []int) {
	logQ := make([]int, depth+1)
	for i := 0; i < depth+1; i++ {
		logQ[i] = 50
	}
	params, _ := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            int(math.Log2(float64(lattigo.n * 2))),
		LogQ:            logQ, // slice of 50s, length depth+1
		LogP:            []int{61, 61, 61},
		LogDefaultScale: 50,
	})
	lattigo.params = &params

	// Generate keys
	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	lattigo.enc = rlwe.NewEncryptor(params, pk)
	lattigo.ecd = ckks.NewEncoder(params)
	lattigo.dec = rlwe.NewDecryptor(params, sk)

	eval := ckks.NewEvaluator(params, evk)

	// Generate rotation keys
	galEls := make([]uint64, len(rots))
	for i, rot := range rots {
		galEls[i] = params.GaloisElement(rot)
	}
	lattigo.eval = eval.WithKey(rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(galEls, sk)...))
}

func parseFloatArray(s string) []float64 {
	s = strings.Trim(s, "[]")
	if s == "" {
		return nil
	}
	parts := strings.Split(s, ",")
	result := make([]float64, len(parts))
	for i, p := range parts {
		result[i], _ = strconv.ParseFloat(strings.TrimSpace(p), 64)
	}
	return result
}

func parseIntArray(s string) []int {
	s = strings.Trim(s, "[]")
	if s == "" {
		return nil
	}
	parts := strings.Split(s, ",")
	result := make([]int, len(parts))
	for i, p := range parts {
		result[i], _ = strconv.Atoi(strings.TrimSpace(p))
	}
	return result
}

func (lattigo *LattigoFHE) encode(values []float64) *rlwe.Ciphertext {
	pack := ckks.NewPlaintext(*lattigo.params, lattigo.params.MaxLevel())
	lattigo.ecd.Encode(values, pack)
	ct, _ := lattigo.enc.EncryptNew(pack)
	return ct
}

func (lattigo *LattigoFHE) evalAdd(ct1, ct2 *rlwe.Ciphertext) *rlwe.Ciphertext {
	ct, _ := lattigo.eval.AddNew(ct1, ct2)
	return ct
}

func (lattigo *LattigoFHE) evalMul(ct1, ct2 *rlwe.Ciphertext) *rlwe.Ciphertext {
	ct, _ := lattigo.eval.MulRelinNew(ct1, ct2)
	lattigo.eval.Rescale(ct, ct)
	return ct
}

func (lattigo *LattigoFHE) evalRot(ct1 *rlwe.Ciphertext, k int) *rlwe.Ciphertext {
	ct, _ := lattigo.eval.RotateNew(ct1, k)
	return ct
}

func (lattigo *LattigoFHE) evalOp(term *Term, metadata string) *rlwe.Ciphertext {
	switch term.op {
	case "FHEOp.PACK":
		values := parseFloatArray(metadata)
		return lattigo.encode(values)
	case "FHEOp.MASK":
		values := parseFloatArray(metadata)
		return lattigo.encode(values)
	case "FHEOp.ADD":
		return lattigo.evalAdd(lattigo.env[term.children[0]], lattigo.env[term.children[1]])
	case "FHEOp.MUL":
		return lattigo.evalMul(lattigo.env[term.children[0]], lattigo.env[term.children[1]])
	case "FHEOp.ROT":
		rot, _ := strconv.Atoi(metadata)
		return lattigo.evalRot(lattigo.env[term.children[0]], rot)
	default:
		return nil
	}
}

func (lattigo *LattigoFHE) preprocess(operations []string) {
	for _, line := range operations {
		lineNum, term, metadata := lattigo.parseOperation(line)
		if lineNum == -1 {
			continue
		}
		switch term.op {
		case "FHEOp.PACK":
			pt := parseFloatArray(metadata)
			if !term.secret {
				lattigo.ptEnv[lineNum] = pt
			}
			lattigo.env[lineNum] = lattigo.encode(pt)
		case "FHEOp.MASK":
			pt := parseFloatArray(metadata)
			lattigo.ptEnv[lineNum] = pt
			lattigo.env[lineNum] = lattigo.encode(pt)
		case "FHEOp.ADD":
			if a, oka := lattigo.ptEnv[term.children[0]]; oka {
				if b, okb := lattigo.ptEnv[term.children[1]]; okb {
					pt := make([]float64, lattigo.n)
					for i := 0; i < lattigo.n; i++ {
						pt[i] = a[i] + b[i]
					}
					lattigo.ptEnv[lineNum] = pt
					lattigo.env[lineNum] = lattigo.encode(pt)
				}
			}
		case "FHEOp.MUL":
			if a, oka := lattigo.ptEnv[term.children[0]]; oka {
				if b, okb := lattigo.ptEnv[term.children[1]]; okb {
					pt := make([]float64, lattigo.n)
					for i := 0; i < lattigo.n; i++ {
						pt[i] = a[i] * b[i]
					}
					lattigo.ptEnv[lineNum] = pt
					lattigo.env[lineNum] = lattigo.encode(pt)
				}
			}
		case "FHEOp.ROT":
			if a, oka := lattigo.ptEnv[term.children[0]]; oka {
				rot, _ := strconv.Atoi(metadata)
				pt := make([]float64, lattigo.n)
				for i := 0; i < lattigo.n; i++ {
					index := ((i+rot)%lattigo.n + lattigo.n) % lattigo.n
					pt[i] = a[index]
				}
				lattigo.ptEnv[lineNum] = pt
				lattigo.env[lineNum] = lattigo.encode(pt)
			}
		}
	}
}

func (lattigo *LattigoFHE) decryptToPlaintext(ct *rlwe.Ciphertext) []float64 {
	pt := lattigo.dec.DecryptNew(ct)
	decoded := make([]float64, lattigo.n)
	lattigo.ecd.Decode(pt, decoded)
	return decoded
}

type PrecisionStats struct {
	AvgPrecision float64
	StdDeviation float64
	MinPrecision float64
	MaxPrecision float64
}

func (lattigo *LattigoFHE) doPrecisionStats(lineNum int, term *Term, metadata string) (PrecisionStats, []float64) {
	want := make([]float64, lattigo.n)
	if want, ok := lattigo.ptEnv[lineNum]; !ok {
		switch term.op {
		case "FHEOp.ADD":
			a := lattigo.ptEnv[term.children[0]]
			b := lattigo.ptEnv[term.children[1]]
			for i := 0; i < min(len(a), len(b)); i++ {
				want[i] = a[i] + b[i]
			}
		case "FHEOp.MUL":
			a := lattigo.ptEnv[term.children[0]]
			b := lattigo.ptEnv[term.children[1]]
			for i := 0; i < min(len(a), len(b)); i++ {
				want[i] = a[i] * b[i]
			}
		case "FHEOp.ROT":
			rot, _ := strconv.Atoi(metadata)
			a := lattigo.ptEnv[term.children[0]]
			for i := 0; i < lattigo.n; i++ {
				index := ((i+rot)%lattigo.n + lattigo.n) % lattigo.n
				want[i] = a[index]
			}
		}
	}

	stats := ckks.GetPrecisionStats(*lattigo.params, lattigo.ecd, lattigo.dec, want, lattigo.env[lineNum], 0, false)

	precStats := PrecisionStats{
		AvgPrecision: stats.AVGLog2Prec.Real,
		StdDeviation: stats.STDLog2Prec.Real,
		MinPrecision: stats.MINLog2Prec.Real,
		MaxPrecision: stats.MAXLog2Prec.Real,
	}
	return precStats, want
}

func (lattigo *LattigoFHE) runInstructions(operations []string, statsPerLine bool) ([]*rlwe.Ciphertext, []float64, []PrecisionStats, time.Duration, error) {
	results := make([]*rlwe.Ciphertext, len(operations))
	allStats := make([]PrecisionStats, len(operations))
	startTime := time.Now()
	want := make([]float64, lattigo.n)

	for _, line := range operations {
		lineNum, term, metadata := lattigo.parseOperation(line)
		if lineNum == -1 {
			continue
		}
		if _, ok := lattigo.env[lineNum]; !ok {
			lattigo.env[lineNum] = lattigo.evalOp(term, metadata)
		}
		results[lineNum] = lattigo.env[lineNum]
		if statsPerLine {
			allStats[lineNum], want = lattigo.doPrecisionStats(lineNum, term, metadata)
		}
	}
	runtime := time.Since(startTime)

	return results, want, allStats, runtime, nil
}

func accurate(expected, got []float64) bool {
	const epsilon = 1e-6
	for i := 0; i < len(expected); i++ {
		if math.Abs(expected[i]-got[i]) > epsilon {
			return false
		}
	}
	return true
}

func (lattigo *LattigoFHE) Run(statsPerLine bool) error {
	expected_str, operations, err := readFile(lattigo.instructionsPath)
	if err != nil {
		return fmt.Errorf("error reading file: %v", err)
	}
	expected := parseFloatArray(expected_str)

	rots := findUniqueRots(operations)
	lattigo.createContext(8, rots)
	lattigo.preprocess(operations)

	results, _, stats, runtime, err := lattigo.runInstructions(operations, statsPerLine)
	if err != nil {
		return fmt.Errorf("error running instructions: %v", err)
	}

	lastResult := results[len(results)-1]
	pt_results := lattigo.decryptToPlaintext(lastResult)

	rounded := make([]float64, len(pt_results))
	for i, v := range pt_results {
		rounded[i] = math.Round(v)
	}
	fmt.Printf("\nOverall Statistics:\n")
	if accurate(expected, rounded) {
		fmt.Println("Passed! ")
	} else {
		fmt.Println("Failed... ")
		for i := 0; i < len(expected); i++ {
			fmt.Printf("Difference: %v\n", expected[i]-rounded[i])
		}
	}

	// Print stats
	if statsPerLine {
		var totalAvg, totalStd float64
		count := 0
		for _, stat := range stats {
			if stat.AvgPrecision != 0 {
				totalAvg += stat.AvgPrecision
				totalStd += stat.StdDeviation
				count++
			}
		}
		fmt.Printf("Average Precision: %.2f bits\n", totalAvg/float64(count))
		fmt.Printf("Average Std Deviation: %.2f bits\n", totalStd/float64(count))
	} else {
		finalStats := ckks.GetPrecisionStats(*lattigo.params, lattigo.ecd, lattigo.dec, expected, lastResult, 0, false)
		fmt.Printf("Final Result Precision: %.2f bits\n", finalStats.AVGLog2Prec.Real)
		fmt.Printf("Final Result Std Deviation: %.2f bits\n", finalStats.STDLog2Prec.Real)
	}
	fmt.Printf("Runtime: %v\n", runtime)

	return nil
}

func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage of %s:\n", os.Args[0])
		fmt.Fprintf(os.Stderr, "  Runs Lattigo FHE operations from Rotom instructions file\n\n")
		flag.PrintDefaults()
	}

	var n int
	var statsPerLine bool
	var instructionsPath string
	flag.IntVar(&n, "n", 4096, "The size of the input")
	flag.BoolVar(&statsPerLine, "spl", false, "Whether to get stats per line")
	flag.StringVar(&instructionsPath, "i", "/home/ubuntu/ajxi/fhe_compiler/instructions/fhe_terms.txt", "The path to the instructions file")
	flag.Parse()

	fhe := NewLattigoFHE(n, instructionsPath)
	if err := fhe.Run(statsPerLine); err != nil {
		fmt.Println(err)
	}
}
