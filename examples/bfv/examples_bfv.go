package main

import (
	"fmt"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/bits"
)

func obliviousRiding() {

	// This example will simulate a situation where an anonymous rider
	// wants to find the closest available rider within a given area.
	// The application is inspired by the paper https://oride.epfl.ch/
	//
	// 		A. Pham, I. Dacosta, G. Endignoux, J. Troncoso-Pastoriza,
	//		K. Huguenin, and J.-P. Hubaux. ORide: A Privacy-Preserving
	//		yet Accountable Ride-Hailing Service. In Proceedings of the
	//		26th USENIX Security Symposium, Vancouver, BC, Canada, August 2017.
	//
	// Each area is represented as a rectangular grid where each driver
	// anyonymously signs in (i.e. the server only knows the driver is located
	// in the area).
	//
	// First, the rider generates an ephemeral key pair (riderSk, riderPk), which he
	// uses to encrypt his coordinates. He then sends the tuple (riderPk, enc(coordinates))
	// to the server handling the area he is in.
	//
	// Once the public key and the encrypted rider coordinates of the rider
	// have been received by the server, the rider's public key is transferred
	// to all the drivers within the area, with a randomized different index
	// for each of them, that indicates in which coefficient each driver must
	// encode his coordinates.
	//
	// Each driver encodes his coordinates in the designated coefficient and
	// uses the received public key to encrypt his encoded coordinates.
	// He then sends back the encrypted coordinates to the server.
	//
	// Once the encrypted coordinates of the drivers have been received, the server
	// homomorphically computes the squared distance: (x0 - x1)^2 + (y0 - y1)^2 between
	// the rider and each of the drivers, and sends back the encrypted result to the rider.
	//
	// The rider decrypts the result and chooses the closest driver.

	// Number of drivers in the area
	nbDrivers := uint64(2048) //max is N

	// BFV parameters (128 bit security)
	params := bfv.DefaultParams[0]

	// Plaintext modulus
	params.T = 0x3ee0001

	encoder := bfv.NewEncoder(&params)

	bfvContext := bfv.NewContextWithParam(&params)

	// Rider's keygen
	kgen := bfvContext.NewKeyGenerator()

	riderSk, riderPk := kgen.NewKeyPair()

	decryptor := bfv.NewDecryptor(riderSk, &params)

	encryptorRiderPk := bfv.NewEncryptorFromPk(riderPk, &params)

	encryptorRiderSk := bfv.NewEncryptorFromSk(riderSk, &params)

	evaluator := bfv.NewEvaluator(&params)

	fmt.Println("============================================")
	fmt.Println("Homomorphic computations on batched integers")
	fmt.Println("============================================")
	fmt.Println()
	fmt.Printf("Parameters : N=%d, T=%d, logQ = %d (%d limbs), sigma = %f \n",
		bfvContext.N(), bfvContext.T(), bfvContext.LogQ(), len(params.Qi), bfvContext.Sigma())
	fmt.Println()

	maxvalue := uint64(math.Sqrt(float64(params.T)))    // max values = floor(sqrt(plaintext modulus))
	mask := uint64(1<<uint64(bits.Len64(maxvalue))) - 1 // binary mask uperbound for the uniform sampling

	fmt.Printf("Generating %d driversData and 1 Rider randomly positioned on a grid of %d x %d units \n",
		nbDrivers, maxvalue, maxvalue)
	fmt.Println()

	// Rider coordinates [x, y, x, y, ....., x, y]
	riderPosX, riderPosY := ring.RandUniform(maxvalue, mask), ring.RandUniform(maxvalue, mask)

	Rider := make([]uint64, params.N)
	for i := uint64(0); i < nbDrivers; i++ {
		Rider[(i << 1)] = riderPosX
		Rider[(i<<1)+1] = riderPosY
	}

	ringCtx := bfv.NewRingContext(&params)

	riderPlaintext := bfv.NewPlaintext(ringCtx)
	encoder.EncodeUint(Rider, riderPlaintext)

	// driversData coordinates [0, 0, ..., x, y, ..., 0, 0]
	driversData := make([][]uint64, nbDrivers)

	driversPlaintexts := make([]*bfv.Plaintext, nbDrivers)
	for i := uint64(0); i < nbDrivers; i++ {
		driversData[i] = make([]uint64, params.N)
		driversData[i][(i << 1)] = ring.RandUniform(maxvalue, mask)
		driversData[i][(i<<1)+1] = ring.RandUniform(maxvalue, mask)
		driversPlaintexts[i] = bfv.NewPlaintext(ringCtx)
		encoder.EncodeUint(driversData[i], driversPlaintexts[i])
	}

	fmt.Printf("Encrypting %d driversData (x, y) and 1 Rider (%d, %d) \n",
		nbDrivers, riderPosX, riderPosY)
	fmt.Println()

	RiderCiphertext := encryptorRiderSk.EncryptNew(riderPlaintext)

	DriversCiphertexts := make([]*bfv.Ciphertext, nbDrivers)
	for i := uint64(0); i < nbDrivers; i++ {
		DriversCiphertexts[i] = encryptorRiderPk.EncryptNew(driversPlaintexts[i])
	}

	fmt.Println("Computing encrypted distance = ((CtD1 + CtD2 + CtD3 + CtD4...) - CtR)^2 ...")
	fmt.Println()

	evaluator.Neg(RiderCiphertext, RiderCiphertext)
	for i := uint64(0); i < nbDrivers; i++ {
		evaluator.Add(RiderCiphertext, DriversCiphertexts[i], RiderCiphertext)
	}

	result := encoder.DecodeUint(decryptor.DecryptNew(evaluator.MulNew(RiderCiphertext, RiderCiphertext)))

	minIndex, minPosX, minPosY, minDist := uint64(0), params.T, params.T, params.T

	errors := 0

	for i := uint64(0); i < nbDrivers; i++ {

		driverPosX, driverPosY := driversData[i][i<<1], driversData[i][(i<<1)+1]

		computedDist := result[i<<1] + result[(i<<1)+1]
		expectedDist := distance(driverPosX, driverPosY, riderPosX, riderPosY)

		if computedDist == expectedDist {
			if computedDist < minDist {
				minIndex = i
				minPosX, minPosY = driverPosX, driverPosY
				minDist = computedDist
			}
		} else {
			errors++
		}

		if i < 4 || i > nbDrivers-5 {
			fmt.Printf("Distance with Driver %d : %8d = (%4d - %4d)^2 + (%4d - %4d)^2 --> correct: %t\n",
				i, computedDist, driverPosX, riderPosX, driverPosY, riderPosY, computedDist == expectedDist)
		}

		if i == nbDrivers>>1 {
			fmt.Println("...")
		}
	}

	fmt.Printf("\nFinished with %.2f%% errors\n\n", 100*float64(errors)/float64(nbDrivers))
	fmt.Printf("Closest Driver to Rider is n°%d (%d, %d) with a distance of %d units\n",
		minIndex, minPosX, minPosY, uint64(math.Sqrt(float64(minDist))))
}

func distance(a, b, c, d uint64) uint64 {
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
	obliviousRiding()
}
