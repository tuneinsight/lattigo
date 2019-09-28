package main

import (
	"fmt"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"log"
	"math"
	"math/bits"
)

var N uint64
var T uint64
var Qi []uint64
var Pi []uint64
var Sigma float64
var bfvContext *bfv.BfvContext

func ObliviousRiding() {

	NbDrivers := uint64(2048)

	params := bfv.DefaultParams[0]

	params.T = 67084289

	bfvContext, err := bfv.NewBfvContextWithParam(&params)
	if err != nil {
		log.Fatal(err)
	}

	kgen := bfvContext.NewKeyGenerator()

	Sk, Pk := kgen.NewKeyPair()

	Decryptor, err := bfvContext.NewDecryptor(Sk)
	if err != nil {
		log.Fatal(err)
	}

	EncryptorPk, err := bfvContext.NewEncryptorFromPk(Pk)
	if err != nil {
		log.Fatal(err)
	}

	EncryptorSk, err := bfvContext.NewEncryptorFromSk(Sk)
	if err != nil {
		log.Fatal(err)
	}

	evalKey := kgen.NewRelinKey(Sk, 1, 60)

	evaluator := bfvContext.NewEvaluator()

	fmt.Println("===========================================")
	fmt.Println("Homomorphic computations on batched integers")
	fmt.Println("===========================================")
	fmt.Println()
	fmt.Printf("Parameters : N=%d, T=%d, logQ = %d (%d limbs), sigma = %f \n", bfvContext.N(), bfvContext.T(), bfvContext.LogQ(), len(params.Qi), bfvContext.Sigma())
	fmt.Println()

	maxvalue := uint64(math.Sqrt(float64(params.T))) // max values = floor(sqrt(plaintext modulus))
	mask := uint64(1<<uint64(bits.Len64(maxvalue))) - 1

	fmt.Printf("Generating %d Drivers and 1 Rider randomly positioned on a grid of %d x %d units \n", NbDrivers, maxvalue, maxvalue)
	fmt.Println()

	Drivers := make([][]uint64, params.N>>1)
	Rider := make([]uint64, params.N)

	riderposition := []uint64{ring.RandUniform(maxvalue, mask), ring.RandUniform(maxvalue, mask)}

	for i := uint64(0); i < params.N>>1; i++ {
		Drivers[i] = make([]uint64, params.N)
		Drivers[i][(i << 1)] = ring.RandUniform(maxvalue, mask)
		Drivers[i][(i<<1)+1] = ring.RandUniform(maxvalue, mask)

		Rider[(i << 1)] = riderposition[0]
		Rider[(i<<1)+1] = riderposition[1]
	}

	batchEncoder, err := bfvContext.NewBatchEncoder()
	if err != nil {
		log.Fatal(err)
	}

	DriversPlaintexts := make([]*bfv.Plaintext, NbDrivers)

	for i := uint64(0); i < NbDrivers; i++ {
		DriversPlaintexts[i] = bfvContext.NewPlaintext()
		batchEncoder.EncodeUint(Drivers[i], DriversPlaintexts[i])
	}

	RiderPlaintext := bfvContext.NewPlaintext()
	batchEncoder.EncodeUint(Rider, RiderPlaintext)

	fmt.Printf("Encrypting %d Drivers (x, y) and 1 Rider (%d, %d) \n", NbDrivers, riderposition[0], riderposition[1])
	fmt.Println()

	DriversCiphertexts := make([]*bfv.Ciphertext, NbDrivers)

	for i := uint64(0); i < NbDrivers; i++ {
		if DriversCiphertexts[i], err = EncryptorPk.EncryptNew(DriversPlaintexts[i]); err != nil {
			log.Fatal(err)
		}
	}

	RiderCiphertext, err := EncryptorSk.EncryptNew(RiderPlaintext)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Compute encrypted Distance = ((CtD1 + CtD2 + CtD3 + CtD4...) - CtR)^2 ...")
	fmt.Println()

	if err := evaluator.Neg(RiderCiphertext, RiderCiphertext); err != nil {
		log.Fatal(err)
	}

	for i := uint64(0); i < NbDrivers; i++ {
		if err := evaluator.Add(RiderCiphertext, DriversCiphertexts[i], RiderCiphertext); err != nil {
			log.Fatal(err)
		}
	}

	res, _ := evaluator.MulNew(RiderCiphertext, RiderCiphertext)
	RiderCiphertext = res.Ciphertext()

	RiderCiphertext, err = evaluator.RelinearizeNew(RiderCiphertext, evalKey)
	if err != nil {
		log.Fatal(err)
	}

	result := batchEncoder.DecodeUint(Decryptor.DecryptNew(RiderCiphertext))

	errors := 0
	closest := []uint64{0, params.T, params.T, params.T}

	for i := uint64(0); i < NbDrivers; i++ {

		r1, r1exp := result[i<<1]+result[(i<<1)+1], Distance(Drivers[i][i<<1], Drivers[i][(i<<1)+1], Rider[0], Rider[1])

		if r1 == r1exp {
			if closest[3] > r1 {
				closest[0] = i
				closest[1] = Drivers[i][i<<1]
				closest[2] = Drivers[i][(i<<1)+1]
				closest[3] = r1exp
			}
		}

		if r1 != r1exp {
			errors += 1
		}

		if i < 4 || i > NbDrivers-5 {
			fmt.Printf("Distance with Driver %d : %8d = (%4d - %4d)^2 + (%4d - %4d)^2: %t \n", i, r1, Drivers[i][i<<1], Rider[0], Drivers[i][(i<<1)+1], Rider[1], r1 == r1exp)
		}

		if i == NbDrivers>>1 {
			fmt.Println("...")
		}
	}

	fmt.Println()

	fmt.Printf("Finished with %.2f%% errors \n", 100*float64(errors)/float64(NbDrivers))

	fmt.Println()

	fmt.Printf("Closest Driver to Rider is n°%d (%d, %d) with a distance of %d units \n", closest[0], closest[1], closest[2], uint64(math.Sqrt(float64(closest[3]))))
}

func Distance(a, b, c, d uint64) uint64 {
	if a > c {
		a, c = c, a
	}
	if b > d {
		b, d = d, b
	}
	x, y := a-c, b-d
	return x*x + y*y
}

func main() {
	ObliviousRiding()
}
