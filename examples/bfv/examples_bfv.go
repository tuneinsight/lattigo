package main

import (
	"fmt"
	"log"
	"math"
	"math/bits"

	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

func obliviousRiding() {

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
	nbDrivers := 2048 // N/2?

	// BFV parameters (128 bit security)
	params := bfv.DefaultParams[0]

	// Plaintext modulus
	params.T = 67084289 // Where does that come from? :)

	bfvContext, err := bfv.NewBfvContextWithParam(&params)
	if err != nil {
		log.Fatal(err)
	}

	batchEncoder, err := bfvContext.NewBatchEncoder()
	if err != nil {
		log.Fatal(err)
	}

	// Rider's keygen
	kgen := bfvContext.NewKeyGenerator()

	riderSk, riderPk := kgen.NewKeyPair()

	decryptor, err := bfvContext.NewDecryptor(riderSk)
	if err != nil {
		log.Fatal(err)
	}

	encryptorPk, err := bfvContext.NewEncryptorFromPk(riderPk)
	if err != nil {
		log.Fatal(err)
	}

	encryptorSk, err := bfvContext.NewEncryptorFromSk(riderSk)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("============================================")
	fmt.Println("Homomorphic computations on batched integers")
	fmt.Println("============================================")
	fmt.Println()
	fmt.Printf("Parameters : N=%d, T=%d, logQ = %d (%d limbs), sigma = %f\n\n",
		bfvContext.N(), bfvContext.T(), bfvContext.LogQ(), len(params.Qi), bfvContext.Sigma())

	maxValue := uint64(math.Sqrt(float64(params.T))) // max values = floor(sqrt(plaintext modulus))
	mask := uint64(1<<uint64(bits.Len64(maxValue))) - 1

	fmt.Printf("Generating %d Drivers and 1 Rider randomly positioned on a grid of %d x %d units\n\n",
		nbDrivers, maxValue, maxValue)

	riderPosX, riderPosY := ring.RandUniform(maxValue, mask), ring.RandUniform(maxValue, mask)

	// Encode the rider position (x, y) as an N-value array [x, y, x, y, ...,
	// x, y], as (N-1) consecutive copies
	riderData := make([]uint64, 2*nbDrivers)
	for i := 0; i < nbDrivers; i++ {
		riderData[i<<1] = riderPosX
		riderData[(i<<1)+1] = riderPosY
	}

	riderPlaintext := bfvContext.NewPlaintext()
	batchEncoder.EncodeUint(riderData, riderPlaintext)

	// Encode each driver i's position (x_i, y_i) as an N-value array [0, 0,
	// ..., x_i, y_i, ..., 0, 0] with the values at the i-th coordinates (i.e.
	// positions 2i and 2i+1)
	driversData := make([][]uint64, nbDrivers)
	driversPlaintexts := make([]*bfv.Plaintext, nbDrivers)
	for i := 0; i < nbDrivers; i++ {
		driverData := make([]uint64, 2*nbDrivers)
		driverData[i<<1] = ring.RandUniform(maxValue, mask)
		driverData[(i<<1)+1] = ring.RandUniform(maxValue, mask)

		driversPlaintexts[i] = bfvContext.NewPlaintext()
		batchEncoder.EncodeUint(driverData, driversPlaintexts[i])

		driversData[i] = driverData
	}

	fmt.Printf("Encrypting %d Drivers (x, y) and 1 Rider (%d, %d)\n\n",
		nbDrivers, riderPosX, riderPosY)

	// Encryption performed by the rider
	riderCiphertext, err := encryptorSk.EncryptNew(riderPlaintext)
	if err != nil {
		log.Fatal(err)
	}

	// Encryption performed by each driver individually
	driversCiphertexts := make([]*bfv.Ciphertext, nbDrivers)
	for i := 0; i < nbDrivers; i++ {
		driversCiphertexts[i], err = encryptorPk.EncryptNew(driversPlaintexts[i])
		if err != nil {
			log.Fatal(err)
		}
	}

	fmt.Println("Computing encrypted Distance = ((CtD1 + CtD2 + CtD3 + CtD4...) - CtR)^2 ...")
	fmt.Println()

	// Computation performed by the server
	evaluator := bfvContext.NewEvaluator()

	// tmpCiphertext = -riderCiphertext
	// --> [-x, -y, -x, -y, ..., -x, -y]
	tmpCiphertext, err := evaluator.NegNew(riderCiphertext)
	if err != nil {
		log.Fatal(err)
	}

	// tmpCiphertext += driversCiphertexts[i] for all drivers
	// --> [x_0-x, y_0-x, x_1-x, y_1-x, ..., x_n-x, y_n-y]
	for i := 0; i < nbDrivers; i++ {
		err := evaluator.Add(tmpCiphertext, driversCiphertexts[i], tmpCiphertext)
		if err != nil {
			log.Fatal(err)
		}
	}

	// resultCiphertext = tmpCiphertext * tmpCiphertext (element-wise multiplication?)
	// --> [(x_0-x)^2, (y_0-y)^2, (x_1-x)^2, (y_1-y)^2, ..., (x_n-x)^2, (y_n-y)^2]
	resultCiphertext, err := evaluator.MulNew(tmpCiphertext, tmpCiphertext)
	if err != nil {
		log.Fatal(err)
	}

	// Decryption performed by the rider
	resultPlaintext := decryptor.DecryptNew(resultCiphertext)
	result := batchEncoder.DecodeUint(resultPlaintext)

	// Compute the shortest distance, and check that the computations performed on the ciphertext are correct
	minIndex, minPosX, minPosY, minDist := 0, params.T, params.T, params.T
	errors := 0

	for i := 0; i < nbDrivers; i++ {
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
