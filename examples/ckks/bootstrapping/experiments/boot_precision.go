package main

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"runtime"
	"text/template"

	"github.com/ldsec/lattigo/v2/ckks"
)

var err error

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

var minLogSlots = int(4)

var paramSet = flag.Int("paramSet", 1, "index in BootStrappParams")
var nboot = flag.Int("nboot", 1, "number of bootstrapping (on the same ct for successive and on different ct for slotdist)")
var logslot = flag.Uint64("logslot", 15, "number of slots per ciphertext (max number for slotcount)")
var makePlot = flag.Bool("makeplot", false, "output a .tex plot")

func main() {

	flag.Parse()

	if flag.NArg() != 1 {
		fmt.Println("Usage: ./boot_precision [-flags] [successive|slotdist|slotcount]")
		flag.PrintDefaults()
		os.Exit(1)
	}
	exp := flag.Args()[0]

	if _, err := os.Stat("tpl"); *makePlot && os.IsNotExist(err) {
		log.Println("\"tpl\" folder not found with -makeplot, run the program from lattigo/experiments/boot_precision/")
		os.Exit(1)
	}

	btpParams := ckks.DefaultBootstrapParams[*paramSet].Copy()
	params, err := btpParams.Params()
	if err != nil {
		panic(err)
	}

	bReal, bImag := new(bytes.Buffer), new(bytes.Buffer)
	var stats []ckks.PrecisionStats

	switch exp {
	case "successive":
		params.SetLogSlots(*logslot)
		log.Printf("%% program args: paramSet=%d, nboot=%d, hw=%d, logslot=%d\n", *paramSet, *nboot, btpParams.H, *logslot)
		encoder, encryptor, evaluator, decryptor, bootstrapper := instanciateExperiment(params, btpParams)
		log.Println("Generating a plaintext of", params.Slots(), "random values...")
		values := make([]complex128, params.Slots())
		for i := range values {
			values[i] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
		}

		plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
		encoder.Encode(plaintext, values, params.LogSlots())
		ciphertext := encryptor.EncryptNew(plaintext)

		stats = make([]ckks.PrecisionStats, *nboot, *nboot)
		for i := range stats {
			log.Printf("Boot %d", i)
			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(params, encoder, decryptor, values, ciphertext, 0)
			if ciphertext.Scale() != params.Scale() {
				evaluator.SetScale(ciphertext, params.Scale())
			}
			runtime.GC()
		}

		formatSuccessive(stats, bReal, bImag)
		break

	case "slotdist":
		params.SetLogSlots(*logslot)
		log.Printf("%% program args: paramSet=%d, hw=%d, logslot=%d\n", *paramSet, btpParams.H, *logslot)
		encoder, encryptor, _, decryptor, bootstrapper := instanciateExperiment(params, btpParams)
		stats = make([]ckks.PrecisionStats, 1, 1) // Experiment seems stable enough
		for i := range stats {
			values := make([]complex128, params.Slots())
			for i := range values {
				values[i] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
			}

			plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
			encoder.Encode(plaintext, values, params.LogSlots())
			ciphertext := encryptor.EncryptNew(plaintext)

			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(params, encoder, decryptor, values, ciphertext, 0)

			runtime.GC()
		}
		formatSlotDist(stats, *logslot, bReal, bImag)
		break

	case "slotcount":
		stats = make([]ckks.PrecisionStats, int(params.LogN())-minLogSlots, int(params.LogN())-minLogSlots)
		log.Printf("%% program args: paramSet=%d, hw=%d\n", *paramSet, btpParams.H)
		for i, logSloti := 0, uint64(minLogSlots); logSloti <= params.LogN()-1; i, logSloti = i+1, logSloti+1 {
			log.Println("running experiment for logslot =", logSloti)
			params.SetLogSlots(logSloti)

			encoder, encryptor, _, decryptor, bootstrapper := instanciateExperiment(params, btpParams)

			values := make([]complex128, params.Slots())
			for j := range values {
				values[j] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
			}

			plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
			encoder.Encode(plaintext, values, params.LogSlots())
			ciphertext := encryptor.EncryptNew(plaintext)

			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(params, encoder, decryptor, values, ciphertext, 0)

			plaintext = nil
			ciphertext = nil
			bootstrapper = nil
			encoder = nil
			encryptor = nil
			decryptor = nil

			runtime.GC()
		}
		formatSlotCount(stats, bReal, bImag)
		break
	default:
		fmt.Println("Invalid experiment")
		os.Exit(1)
	}
	output(os.Stdout, exp, stats, *makePlot, bReal, bImag)
}

func instanciateExperiment(params *ckks.Parameters, btpParams *ckks.BootstrappingParameters) (encoder ckks.Encoder, encryptor ckks.Encryptor, evaluator ckks.Evaluator, decryptor ckks.Decryptor, bootstrapper *ckks.Bootstrapper) {

	keyGen := ckks.NewKeyGenerator(params)
	sk, pk := keyGen.GenKeyPairSparse(btpParams.H)

	encoder = ckks.NewEncoder(params)
	encryptor = ckks.NewEncryptorFromPk(params, pk)
	decryptor = ckks.NewDecryptor(params, sk)

	evaluator = ckks.NewEvaluator(params, ckks.EvaluationKey{nil, nil})

	log.Println("Generating the keys...")
	btpKey := keyGen.GenBootstrappingKey(params.LogSlots(), btpParams, sk)
	log.Println("Done")
	bootstrapper, err = ckks.NewBootstrapper(params, btpParams, *btpKey)
	if err != nil {
		panic(err)
	}

	return
}

func formatParams(params, succ int, hw, logSlot uint64) string {
	return fmt.Sprintf("%% program args: paramSet=%d, nboot=%d, hw=%d, logslot=%d\n", params, succ, hw, logSlot)
}

func formatSlotCount(stats []ckks.PrecisionStats, wReal, wImag io.Writer) {
	for logSlot, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(wReal, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", logSlot+minLogSlots, real(prec.MeanPrecision), real(prec.MaxPrecision-prec.MeanPrecision), real(prec.MeanPrecision-prec.MinPrecision))
	}
	for logSlot, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(wImag, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", logSlot+minLogSlots, imag(prec.MeanPrecision), imag(prec.MaxPrecision-prec.MeanPrecision), imag(prec.MeanPrecision-prec.MinPrecision))
	}
}

func formatSlotDist(stats []ckks.PrecisionStats, logSlot uint64, wReal, wImag io.Writer) {
	slotCount := 1 << logSlot

	for _, point := range stats[0].RealDist {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(wReal, "(%.2f, %.4f)\n", point.Prec, float64(point.Count)/float64(slotCount))
	}
	for _, point := range stats[0].ImagDist {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(wImag, "(%.2f, %.4f)\n", point.Prec, float64(point.Count)/float64(slotCount))
	}
}

func formatSuccessive(stats []ckks.PrecisionStats, wReal, wImag io.Writer) {
	for i, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(wReal, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", i, real(prec.MedianPrecision), real(prec.MaxPrecision-prec.MedianPrecision), real(prec.MedianPrecision-prec.MinPrecision))
	}
	for i, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(wImag, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", i, imag(prec.MedianPrecision), imag(prec.MaxPrecision-prec.MedianPrecision), imag(prec.MedianPrecision-prec.MinPrecision))
	}
}

func output(out io.Writer, exp string, stats []ckks.PrecisionStats, makePlot bool, rReal, rImag *bytes.Buffer) {
	if makePlot {
		t := template.Must(template.ParseFiles("tpl" + string(os.PathSeparator) + exp + ".tex.tpl"))
		err := t.Execute(out, struct {
			DataReal string
			DataImag string
		}{
			rReal.String(),
			rImag.String(),
		})
		if err != nil {
			panic(err)
		}
		return
	}
	fmt.Fprintln(out, "% Real\n", rReal.String(), "% Imag\n", rImag.String())
}
