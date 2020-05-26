package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/ldsec/lattigo/ckks"
	"log"
	"math"
	"os"
)

var paramSet = flag.Int("paramSet", 1, "index in BootStrappParams")
var nboot = flag.Int("nboot", 1, "number of bootstrapping (on the same ct for successive and on different ct for slotdist)")
var logslot = flag.Uint64("logslot", 10, "number of slots per ciphertext (max number for slotcount)")
var hw = flag.Uint64("hw", 128, "secret key hamming weight")

func main() {

	flag.Parse()

	if flag.NArg() != 1 {
		fmt.Println("Usage: ./boot_precision [-flags] [successive|slotdist|slotcount]")
		flag.PrintDefaults()
		os.Exit(1)
	}
	exp := flag.Args()[0]

	params := ckks.BootstrappParams[*paramSet].Copy()
	params.LogSlots = *logslot

	fmt.Println(formatParams(*paramSet, *nboot, *hw, *logslot))
	switch exp {
	case "successive":
		encoder, encryptor, decryptor, bootstrapper := instanciateExperiment(params)
		log.Println("Generating a plaintext of", params.Slots, "random values...")
		values := make([]complex128, params.Slots)
		for i := range values {
			values[i] = complex(ckks.RandomFloat(-1, 1), ckks.RandomFloat(-1, 1))
		}

		plaintext := ckks.NewPlaintext(&params.Parameters, params.MaxLevel, params.Scale)
		encoder.Encode(plaintext, values, params.Slots)
		ciphertext := encryptor.EncryptNew(plaintext)

		stats := make([]ckks.PrecisionStats, *nboot, *nboot)
		for i := range stats {
			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(&params.Parameters, encoder, decryptor, values, ciphertext)
		}

		fmt.Println(formatSuccessive(stats))
		break


	case "slotdist":

		encoder, encryptor, decryptor, bootstrapper := instanciateExperiment(params)
		stats := make([]ckks.PrecisionStats, *nboot, *nboot)
		for i := range stats {
			values := make([]complex128, params.Slots)
			for i := range values {
				values[i] = complex(ckks.RandomFloat(-1, 1), ckks.RandomFloat(-1, 1))
			}

			plaintext := ckks.NewPlaintext(&params.Parameters, params.MaxLevel, params.Scale)
			encoder.Encode(plaintext, values, params.Slots)
			ciphertext := encryptor.EncryptNew(plaintext)

			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(&params.Parameters, encoder, decryptor, values, ciphertext)
		}

		fmt.Println(formatSlotDist(stats, *logslot))
		break

	case "slotcount":
		stats := make([]ckks.PrecisionStats, *logslot-2, *logslot-2)
		for i, logSloti := 0, uint64(3); logSloti <= *logslot; i, logSloti = i+1, logSloti+1 {
			log.Println("running experiment for logslot =", logSloti)
			params := ckks.BootstrappParams[*paramSet].Copy()
			params.LogSlots = logSloti
			encoder, encryptor, decryptor, bootstrapper := instanciateExperiment(params)

			values := make([]complex128, params.Slots)
			for j := range values {
				values[j] = complex(ckks.RandomFloat(-1, 1), ckks.RandomFloat(-1, 1))
			}

			plaintext := ckks.NewPlaintext(&params.Parameters, params.MaxLevel, params.Scale)
			encoder.Encode(plaintext, values, params.Slots)
			ciphertext := encryptor.EncryptNew(plaintext)

			ciphertext = bootstrapper.Bootstrapp(ciphertext)
			stats[i] = ckks.GetPrecisionStats(&params.Parameters, encoder, decryptor, values, ciphertext)
		}
		fmt.Println(formatSlotCount(stats))
		break
	}
}

func instanciateExperiment(params *ckks.BootParams) (encoder ckks.Encoder, encryptor ckks.Encryptor, decryptor ckks.Decryptor, bootstrapper *ckks.BootContext) {

	if err := params.Gen(); err != nil {
		log.Fatal(err)
	}
	keyGen := ckks.NewKeyGenerator(&params.Parameters)
	sk, pk := keyGen.GenKeyPairSparse(*hw)

	encoder = ckks.NewEncoder(&params.Parameters)
	encryptor = ckks.NewEncryptorFromPk(&params.Parameters, pk)
	decryptor = ckks.NewDecryptor(&params.Parameters, sk)

	bootstrapper = ckks.NewBootContext(params)
	log.Println("Generating the keys...")
	bootstrapper.GenBootKeys(sk)
	return
}

func formatParams(params, succ int, hw, logSlot uint64) string {
	return fmt.Sprintf("%% paramSet=%d, nboot=%d, hw=%d, logslot=%d\n", params, succ, hw, logSlot)
}

func formatSlotCount(stats []ckks.PrecisionStats) string {
	w := new(bytes.Buffer)
	fmt.Fprintln(w, "% Real")
	for logSlot, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(w, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", logSlot+3, math.Log2(1/real(prec.Median)), math.Log2(1/real(prec.Max)), math.Log2(1/real(prec.Min)))
	}
	fmt.Fprintln(w, "% Imag")
	for logSlot, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(w, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", logSlot+3, math.Log2(1/imag(prec.Median)), math.Log2(1/imag(prec.Max)), math.Log2(1/imag(prec.Min)))
	}
	return w.String()
}

func formatSlotDist(stats []ckks.PrecisionStats, logSlot uint64) string {

	aggrReal := make(map[float64]uint64)
	aggrImag := make(map[float64]uint64)

	for _, stat := range stats {
		for precBin, count := range stat.RealDist {
			aggrReal[precBin] += count
		}
		for precBin, count := range stat.ImagDist {
			aggrImag[precBin] += count
		}
	}

	w := new(bytes.Buffer)
	slotCount := 1 << logSlot
	fmt.Fprintln(w, "% Real")
	var tot uint64
	for precBin, count := range aggrReal {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(w, "(%.2f, %.6f)", precBin, float64(count)/float64(slotCount))
		tot += count
	}

	fmt.Fprintln(w, "\n% Imag")
	for precBin, count := range aggrImag {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(w, "(%.2f, %.6f)",  precBin, float64(count)/float64(slotCount))
	}
	fmt.Fprintln(w)
	return w.String()
}

func formatSuccessive(stats []ckks.PrecisionStats) string {
	w := new(bytes.Buffer)
	fmt.Fprintln(w, "% Real")
	for i, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(w, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", i, math.Log2(1/real(prec.Median)), math.Log2(1/real(prec.Max)), math.Log2(1/real(prec.Min)))
	}
	fmt.Fprintln(w, "% Imag")
	for i, prec := range stats {
		// (1,  19.77) += (0, 13.1) -= (0, 4.87)
		fmt.Fprintf(w, "(%d, %.2f) += (0, %.2f) -= (0, %.2f)\n", i, math.Log2(1/imag(prec.Median)), math.Log2(1/imag(prec.Max)), math.Log2(1/imag(prec.Min)))
	}
	return w.String()
}
