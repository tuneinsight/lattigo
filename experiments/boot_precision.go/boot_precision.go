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

var succ = flag.Int("succ", 3, "number of successive bootstrapping")
var slot = flag.Uint64("logSlot", 15, "number of slots per ciphertext")
var params = flag.Int("params", 1, "index in BootStrappParams")
var hw = flag.Uint64("hw", 128, "secret key hamming weight")

func main() {

	flag.Parse()

	if flag.NArg() != 1 {
		fmt.Println("Usage: ./boot_precision [-flags] [successive]")
		os.Exit(1)
	}
	exp := flag.Args()[0]

	n_boot := *succ
	params := ckks.BootstrappParams[*params]
	params.LogSlots = *slot
	params.Gen()

	//ckkscontext := newContext(&params.Parameters)

	keyGen := ckks.NewKeyGenerator(&params.Parameters)
	sk, pk := keyGen.GenKeyPairSparse(*hw)

	encoder := ckks.NewEncoder(&params.Parameters)
	encryptorPk := ckks.NewEncryptorFromPk(&params.Parameters, pk)
	decryptor := ckks.NewDecryptor(&params.Parameters, sk)

	bootstrapper := ckks.NewBootContext(params)
	log.Println("Generating the keys...")
	bootstrapper.GenBootKeys(sk)

	values := make([]complex128, params.Slots)
	for i := range values {
		values[i] = complex(ckks.RandomFloat(-1, 1), ckks.RandomFloat(-1, 1))
	}

	plaintext := ckks.NewPlaintext(&params.Parameters, params.MaxLevel, params.Scale)
	encoder.Encode(plaintext, values, params.Slots)
	ciphertext := encryptorPk.EncryptNew(plaintext)

	stats := make([]ckks.PrecisionStats, n_boot, n_boot)
	for i := range stats {
		ciphertext = bootstrapper.Bootstrapp(ciphertext)
		stats[i] = ckks.GetPrecisionStats(&params.Parameters, encoder, decryptor, values, ciphertext)
	}

	switch exp {
	case "successive":
		fmt.Println(formatSuccessive(stats))
	}
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
