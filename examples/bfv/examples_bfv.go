package main

import (
	"fmt"
	"log"
	"math"
	"math/bits"

	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

var N uint64
var T uint64
var Qi []uint64
var Pi []uint64
var Sigma float64
var bfvContext *bfv.BfvContext

func ObliviousRiding() {

	// This example will simulate a situation where an anonymous rider wants to find the closest available rider within a given area.
	// The application is inspired by the paper https://oride.epfl.ch/
	//
	// 		A. Pham, I. Dacosta, G. Endignoux, J. Troncoso-Pastoriza, K. Huguenin, and J.-P. Hubaux. ORide: A Privacy-Preserving yet Accountable Ride-Hailing Service.
	//		In Proceedings of the 26th USENIX Security Symposium, Vancouver, BC, Canada, August 2017.
	//
	// Each area is represented as a rectangular grid where each driver reports his coordinates in real time to a server assigned to the area.
	//
	// First, the rider generates an ephemeral key pair (sk, pk), which he uses to encrypt his coordinates. He then sends the tuple (pk, enc(coordinates)) to the
	// server handling the area he is in.
	//
	// Once the public key and the encrypted rider coordinates of the rider have been received by the server, the rider's public key is
	// transferred to all the drivers within the area, with a randomized different index for each of them, that indicates in which coefficient each driver must
	// encode his coordinates.
	//
	// Each driver encodes his coordinates in the designated coefficient and uses the received public key to encrypt his encoded coordinates.
	// He then sends back the encrypted coordinates to the server.
	//
	// Once the encrypted coordinates of the drivers have been received, the server homomorphically computes the squared distance: (x0 - x1)^2 + (y0 - y1)^2 between
	// the rider and each of the drivers, and sends back the encrypted result to the rider.
	//
	// The rider decrypts the result and chooses the closest driver.

	// Number of drivers in the area
	NbDrivers := uint64(2048)

	// BFV parameters (128 bit security)
	params := bfv.DefaultParams[0]

	// Plaintext modulus
	params.T = 0x3ee0001

	bfvContext, err := bfv.NewBfvContextWithParam(&params)
	if err != nil {
		log.Fatal(err)
	}

	encoder := bfvContext.NewEncoder()

	// Rider's keygen
	kgen := bfvContext.NewKeyGenerator()

	Sk, Pk := kgen.NewKeyPair()

	Decryptor := bfvContext.NewDecryptor(Sk)

	EncryptorPk := bfvContext.NewEncryptorFromPk(Pk)

	EncryptorSk := bfvContext.NewEncryptorFromSk(Sk)

	evaluator := bfvContext.NewEvaluator()

	fmt.Println("============================================")
	fmt.Println("Homomorphic computations on batched integers")
	fmt.Println("============================================")
	fmt.Println()
	fmt.Printf("Parameters : N=%d, T=%d, logQ = %d (%d limbs), sigma = %f \n", bfvContext.N(), bfvContext.T(), bfvContext.LogQ(), len(params.Qi), bfvContext.Sigma())
	fmt.Println()

	maxvalue := uint64(math.Sqrt(float64(params.T))) // max values = floor(sqrt(plaintext modulus))
	mask := uint64(1<<uint64(bits.Len64(maxvalue))) - 1

	fmt.Printf("Generating %d Drivers and 1 Rider randomly positioned on a grid of %d x %d units \n", NbDrivers, maxvalue, maxvalue)
	fmt.Println()

	Drivers := make([][]uint64, params.N>>1)

	// Rider coordinates [x, y, x, y, ....., x, y]
	riderposition := []uint64{ring.RandUniform(maxvalue, mask), ring.RandUniform(maxvalue, mask)}

	Rider := make([]uint64, params.N)
	for i := uint64(0); i < params.N>>1; i++ {
		Rider[(i << 1)] = riderposition[0]
		Rider[(i<<1)+1] = riderposition[1]
	}

	RiderPlaintext := bfvContext.NewPlaintext()
	encoder.EncodeUint(Rider, RiderPlaintext)

	// Drivers coordinates [0, 0, ..., x, y, ..., 0, 0]
	for i := uint64(0); i < params.N>>1; i++ {
		Drivers[i] = make([]uint64, params.N)
		Drivers[i][(i << 1)] = ring.RandUniform(maxvalue, mask)
		Drivers[i][(i<<1)+1] = ring.RandUniform(maxvalue, mask)
	}

	DriversPlaintexts := make([]*bfv.Plaintext, NbDrivers)
	for i := uint64(0); i < NbDrivers; i++ {
		DriversPlaintexts[i] = bfvContext.NewPlaintext()
		encoder.EncodeUint(Drivers[i], DriversPlaintexts[i])
	}

	fmt.Printf("Encrypting %d Drivers (x, y) and 1 Rider (%d, %d) \n", NbDrivers, riderposition[0], riderposition[1])
	fmt.Println()

	RiderCiphertext := EncryptorSk.EncryptNew(RiderPlaintext)

	DriversCiphertexts := make([]*bfv.Ciphertext, NbDrivers)
	for i := uint64(0); i < NbDrivers; i++ {
		DriversCiphertexts[i] = EncryptorPk.EncryptNew(DriversPlaintexts[i])
	}

	fmt.Println("Computing encrypted Distance = ((CtD1 + CtD2 + CtD3 + CtD4...) - CtR)^2 ...")
	fmt.Println()

	evaluator.Neg(RiderCiphertext, RiderCiphertext)
	for i := uint64(0); i < NbDrivers; i++ {
		evaluator.Add(RiderCiphertext, DriversCiphertexts[i], RiderCiphertext)
	}

	result := encoder.DecodeUint(Decryptor.DecryptNew(evaluator.MulNew(RiderCiphertext, RiderCiphertext)))

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
