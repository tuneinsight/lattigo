package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/ldsec/lattigo/ckks"
	"io"
	"log"
	"math/rand"
	"os"
	"text/template"
)

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

var paramSet = flag.Int("paramSet", 1, "index in BootStrappParams")
var nboot = flag.Int("nboot", 1, "number of bootstrapping (on the same ct for successive and on different ct for slotdist)")
var logslot = flag.Uint64("logslot", 8, "number of slots per ciphertext (max number for slotcount)")
var hw = flag.Uint64("hw", 192, "secret key hamming weight")
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

	params := ckks.BootstrappParams[*paramSet].Copy()
	params.LogSlots = *logslot
	bReal, bImag := new(bytes.Buffer), new(bytes.Buffer)
	fmt.Println(formatParams(*paramSet, *nboot, *hw, *logslot))
	var stats []ckks.PrecisionStats
	switch exp {
	case "successive":
		encoder, encryptor, evaluator, decryptor, bootstrapper := instanciateExperiment(params)
		log.Println("Generating a plaintext of", params.Slots, "random values...")
		values := make([]complex128, params.Slots)
		for i := range values {
			values[i] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
		}

		plaintext := ckks.NewPlaintext(&params.Parameters, params.MaxLevel, params.Scale)
		encoder.Encode(plaintext, values, params.Slots)
		ciphertext := encryptor.EncryptNew(plaintext)

		stats = make([]ckks.PrecisionStats, *nboot, *nboot)
		for i := range stats {
			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(&params.Parameters, encoder, decryptor, values, ciphertext)
			fmt.Println(stats[i])
			evaluator.SetScale(ciphertext, params.Scale)
			fmt.Println(ciphertext.Level(), ciphertext.Scale())
		}

		formatSuccessive(stats, bReal, bImag)
		break

	case "slotdist":

		encoder, encryptor, _, decryptor, bootstrapper := instanciateExperiment(params)
		stats = make([]ckks.PrecisionStats, 1, 1) // Experiment seems stable enough
		for i := range stats {
			values := make([]complex128, params.Slots)
			for i := range values {
				values[i] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
			}

			plaintext := ckks.NewPlaintext(&params.Parameters, params.MaxLevel, params.Scale)
			encoder.Encode(plaintext, values, params.Slots)
			ciphertext := encryptor.EncryptNew(plaintext)

			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(&params.Parameters, encoder, decryptor, values, ciphertext)
		}
		formatSlotDist(stats, *logslot, bReal, bImag)
		break

	case "slotcount":
		stats = make([]ckks.PrecisionStats, *logslot-2, *logslot-2)
		for i, logSloti := 0, uint64(4); logSloti <= *logslot; i, logSloti = i+1, logSloti+1 {
			log.Println("running experiment for logslot =", logSloti)
			params := ckks.BootstrappParams[*paramSet].Copy()
			params.LogSlots = logSloti
			encoder, encryptor, _, decryptor, bootstrapper := instanciateExperiment(params)

			values := make([]complex128, params.Slots)
			for j := range values {
				values[j] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
			}

			plaintext := ckks.NewPlaintext(&params.Parameters, params.MaxLevel, params.Scale)
			encoder.Encode(plaintext, values, params.Slots)
			ciphertext := encryptor.EncryptNew(plaintext)

			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(&params.Parameters, encoder, decryptor, values, ciphertext)
		}
		formatSlotCount(stats, bReal, bImag)
		break
	default:
		fmt.Println("Invalid experiment")
		os.Exit(1)
	}
	output(os.Stdout, exp, stats, *makePlot, bReal, bImag)
}

func instanciateExperiment(params *ckks.BootParams) (encoder ckks.Encoder, encryptor ckks.Encryptor, evaluator ckks.Evaluator, decryptor ckks.Decryptor, bootstrapper *ckks.BootContext) {

	if err := params.Gen(); err != nil {
		log.Fatal(err)
	}
	keyGen := ckks.NewKeyGenerator(&params.Parameters)
	sk, pk := keyGen.GenKeyPairSparse(*hw)

	encoder = ckks.NewEncoder(&params.Parameters)
	encryptor = ckks.NewEncryptorFromPk(&params.Parameters, pk)
	decryptor = ckks.NewDecryptor(&params.Parameters, sk)

	evaluator = ckks.NewEvaluator(&params.Parameters)

	bootstrapper = ckks.NewBootContext(params)
	log.Println("Generating the keys...")
	bootstrapper.GenBootKeys(sk)
	return
}

func formatParams(params, succ int, hw, logSlot uint64) string {
	return fmt.Sprintf("%% program args: paramSet=%d, nboot=%d, hw=%d, logslot=%d\n", params, succ, hw, logSlot)
}

func formatSlotCount(stats []ckks.PrecisionStats, wReal, wImag io.Writer) {
	for logSlot, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(wReal, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", logSlot+3, real(prec.MedianPrecision), real(prec.MaxPrecision-prec.MedianPrecision), real(prec.MedianPrecision-prec.MinPrecision))
	}
	for logSlot, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(wImag, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", logSlot+3, real(prec.MedianPrecision), real(prec.MaxPrecision-prec.MedianPrecision), real(prec.MedianPrecision-prec.MinPrecision))
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
