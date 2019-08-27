package main

import (
	"fmt"
	"github.com/lca1/lattigo/bfv"
	"github.com/lca1/lattigo/ring"
	"log"
	"math/bits"
)

var N uint64
var T uint64
var Qi []uint64
var Pi []uint64
var Sigma float64
var bfvContext *bfv.BfvContext

func Plaintext_Batching() {

	N = 8
	T = 786433
	Qi = []uint64{1152921504053723137, 1152921504050839553}
	Pi = []uint64{576460752308273153, 576460752315482113, 576460752319021057}
	Sigma = float64(3.19)

	bfvContext, err := bfv.NewBfvContextWithParam(N, T, Qi, Pi, Sigma)
	if err != nil {
		log.Fatal(err)
	}

	kgen := bfvContext.NewKeyGenerator()

	Sk := kgen.NewSecretKey()

	Decryptor, err := bfvContext.NewDecryptor(Sk)
	if err != nil {
		log.Fatal(err)
	}

	pks := make([]*bfv.PublicKey, 5, 5)
	encryptors := make([]*bfv.Encryptor, 5, 5)
	for i := range pks {
		pks[i], err = kgen.NewPublicKey(Sk)
		if err != nil {
			log.Fatal(err)
		}
		encryptors[i], err = bfvContext.NewEncryptor(pks[i])
		if err != nil {
			log.Fatal(err)
		}
	}

	evalKey, err := kgen.NewRelinKey(Sk, 2, 60)
	if err != nil {
		log.Fatal(err)
	}

	evaluator, err := bfvContext.NewEvaluator()
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("===========================================")
	fmt.Println("Homomorphic computations on batched integers")
	fmt.Println("===========================================")
	fmt.Println()
	fmt.Printf("Parameters : N=%d, T=%d, Qi = %dx60, sigma = %f \n", N, T, len(Qi), Sigma)
	fmt.Println()

	maxvalue := uint64(880)
	mask := uint64(1<<uint64(bits.Len64(maxvalue))) - 1

	//			       (i, x, y)
	Driver1 := []int64{0, int64(ring.RandUniform(maxvalue, mask)), int64(ring.RandUniform(maxvalue, mask))}
	Driver2 := []int64{1, int64(ring.RandUniform(maxvalue, mask)), int64(ring.RandUniform(maxvalue, mask))}
	Driver3 := []int64{2, int64(ring.RandUniform(maxvalue, mask)), int64(ring.RandUniform(maxvalue, mask))}
	Driver4 := []int64{3, int64(ring.RandUniform(maxvalue, mask)), int64(ring.RandUniform(maxvalue, mask))}
	Rider := []int64{int64(ring.RandUniform(maxvalue, mask)), int64(ring.RandUniform(maxvalue, mask))}

	coeffsD1 := []int64{0, 0, 0, 0, 0, 0, 0, 0}
	coeffsD2 := []int64{0, 0, 0, 0, 0, 0, 0, 0}
	coeffsD3 := []int64{0, 0, 0, 0, 0, 0, 0, 0}
	coeffsD4 := []int64{0, 0, 0, 0, 0, 0, 0, 0}
	coeffsR := []int64{0, 0, 0, 0, 0, 0, 0, 0}

	coeffsD1[Driver1[0]<<1], coeffsD1[1+(Driver1[0])<<1] = Driver1[1], Driver1[2]
	coeffsD2[Driver2[0]<<1], coeffsD2[1+(Driver2[0])<<1] = Driver2[1], Driver2[2]
	coeffsD3[Driver3[0]<<1], coeffsD3[1+(Driver3[0])<<1] = Driver3[1], Driver3[2]
	coeffsD4[Driver4[0]<<1], coeffsD4[1+(Driver4[0])<<1] = Driver4[1], Driver4[2]

	for i := uint64(0); i < N; i += 2 {
		coeffsR[i] = Rider[0]
		coeffsR[i+1] = Rider[1]
	}

	mD1 := bfvContext.NewPlaintext()
	mD2 := bfvContext.NewPlaintext()
	mD3 := bfvContext.NewPlaintext()
	mD4 := bfvContext.NewPlaintext()
	mR := bfvContext.NewPlaintext()

	batchEncoder := bfvContext.NewBatchEncoder()

	batchEncoder.EncodeInt(coeffsD1, mD1)
	batchEncoder.EncodeInt(coeffsD2, mD2)
	batchEncoder.EncodeInt(coeffsD3, mD3)
	batchEncoder.EncodeInt(coeffsD4, mD4)
	batchEncoder.EncodeInt(coeffsR, mR)

	fmt.Printf("Driver 1 : [")
	for i := uint64(0); i < N; i++ {
		fmt.Printf("%8d", coeffsD1[i])
	}
	fmt.Printf(" ] -> Encrypt with Pk1 -> CtD1\n")

	fmt.Printf("Driver 2 : [")
	for i := uint64(0); i < N; i++ {
		fmt.Printf("%8d", coeffsD2[i])
	}

	fmt.Printf(" ] -> Encrypt with Pk2 -> CtD2\n")

	fmt.Printf("Driver 3 : [")
	for i := uint64(0); i < N; i++ {
		fmt.Printf("%8d", coeffsD3[i])
	}

	fmt.Printf(" ] -> Encrypt with Pk3 -> CtD3\n")

	fmt.Printf("Driver 4 : [")
	for i := uint64(0); i < N; i++ {
		fmt.Printf("%8d", coeffsD4[i])
	}

	fmt.Printf(" ] -> Encrypt with Pk4 -> CtD4\n")

	fmt.Printf("Rider    : [")
	for i := uint64(0); i < N; i++ {
		fmt.Printf("%8d", coeffsR[i])
	}
	fmt.Printf(" ] -> Encrypt with Pk5 -> CtR\n")

	fmt.Println()

	CtD1, err := encryptors[0].EncryptNew(mD1)
	if err != nil {
		log.Fatal(err)
	}

	CtD2, err := encryptors[1].EncryptNew(mD2)
	if err != nil {
		log.Fatal(err)
	}

	CtD3, err := encryptors[2].EncryptNew(mD3)
	if err != nil {
		log.Fatal(err)
	}

	CtD4, err := encryptors[3].EncryptNew(mD4)
	if err != nil {
		log.Fatal(err)
	}

	CtR, err := encryptors[4].EncryptNew(mR)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Compute encrypted Distance = ((CtD1 + CtD2 + CtD3 + CtD4) - CtR)^2 ...")

	if err := evaluator.Add(CtD1, CtD2, CtD1); err != nil {
		log.Fatal(err)
	}
	if err := evaluator.Add(CtD1, CtD3, CtD1); err != nil {
		log.Fatal(err)
	}
	if err := evaluator.Add(CtD1, CtD4, CtD1); err != nil {
		log.Fatal(err)
	}

	evaluator.Sub(CtD1, CtR, CtD1)

	CtD1 = evaluator.MulNew(CtD1, CtD1).(*bfv.Ciphertext)

	CtD1, err = evaluator.RelinearizeNew(CtD1, evalKey)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Done!")
	fmt.Println()

	mR, err = Decryptor.DecryptNew(CtD1)
	if err != nil {
		log.Fatal(err)
	}

	result, err := batchEncoder.DecodeUint(mR)
	if err != nil {
		log.Fatal(err)
	}

	r1, r1exp := result[0]+result[1], uint64((coeffsD1[0]-coeffsR[0])*(coeffsD1[0]-coeffsR[0])+(coeffsD1[1]-coeffsR[1])*(coeffsD1[1]-coeffsR[1]))
	r2, r2exp := result[2]+result[3], uint64((coeffsD2[2]-coeffsR[2])*(coeffsD2[2]-coeffsR[2])+(coeffsD2[3]-coeffsR[3])*(coeffsD2[3]-coeffsR[3]))
	r3, r3exp := result[4]+result[5], uint64((coeffsD3[4]-coeffsR[4])*(coeffsD3[4]-coeffsR[4])+(coeffsD3[5]-coeffsR[5])*(coeffsD3[5]-coeffsR[5]))
	r4, r4exp := result[6]+result[7], uint64((coeffsD4[6]-coeffsR[6])*(coeffsD4[6]-coeffsR[6])+(coeffsD4[7]-coeffsR[7])*(coeffsD4[7]-coeffsR[7]))

	fmt.Printf("Distance with Driver %d : %6d = (%3d - %3d)^2 + (%3d - %3d)^2: %t \n", 1, r1, coeffsD1[0], coeffsR[0], coeffsD1[1], coeffsR[1], r1 == r1exp)
	fmt.Printf("Distance with Driver %d : %6d = (%3d - %3d)^2 + (%3d - %3d)^2: %t \n", 2, r2, coeffsD2[2], coeffsR[2], coeffsD2[3], coeffsR[3], r2 == r2exp)
	fmt.Printf("Distance with Driver %d : %6d = (%3d - %3d)^2 + (%3d - %3d)^2: %t \n", 3, r3, coeffsD3[4], coeffsR[4], coeffsD3[5], coeffsR[5], r3 == r3exp)
	fmt.Printf("Distance with Driver %d : %6d = (%3d - %3d)^2 + (%3d - %3d)^2: %t \n", 4, r4, coeffsD4[6], coeffsR[6], coeffsD4[7], coeffsR[7], r4 == r4exp)
	fmt.Println()
}

func Homomorphic_Inner_product() {

	fmt.Println("===================================================================")
	fmt.Println("Homomorphic computations on batched integers (sum of inner product)")
	fmt.Println("===================================================================")
	fmt.Println()

	N = 8192
	T = 281474976317441
	Qi = []uint64{1152921504066306049, 1152921504057917441, 1152921504053723137, 1152921504050839553}
	Pi = []uint64{576460752568975361, 576460752573431809, 576460752580902913, 576460752585490433, 576460752586407937}
	Sigma = float64(3.19)

	fmt.Printf("Parameters : N=%d, T=%d, Qi = %dx60, sigma = %f \n", N, T, len(Qi), Sigma)
	fmt.Println()

	bfvContext = bfv.NewBfvContext()
	if err := bfvContext.SetParameters(N, T, Qi, Pi, Sigma); err != nil {
		log.Fatal(err)
	}

	kgen := bfvContext.NewKeyGenerator()

	encoder := bfvContext.NewBatchEncoder()

	Sk := kgen.NewSecretKey()

	Pk, err := kgen.NewPublicKey(Sk)
	if err != nil {
		log.Fatal(err)
	}

	rlk, err := kgen.NewRelinKey(Sk, 2, 30)
	if err != nil {
		log.Fatal(err)
	}

	rotationKey, err := kgen.NewRotationKeysPow2(Sk, 30, true)
	if err != nil {
		log.Fatal(err)
	}

	Decryptor, err := bfvContext.NewDecryptor(Sk)
	if err != nil {
		log.Fatal(err)
	}

	Encryptor, err := bfvContext.NewEncryptor(Pk)
	if err != nil {
		log.Fatal(err)
	}

	Evaluator, err := bfvContext.NewEvaluator()
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Generating to random arrays (a and b) of 784 integers in the range [-500000, 500000]")
	fmt.Println()
	coeffsD1 := make([]int64, N)
	coeffsD2 := make([]int64, N)
	var tmp0, tmp1, sign1, sign2 int64

	for i := int64(0); i < 784; i++ {
		tmp0 = int64(ring.RandUniform(2, 3))
		tmp1 = int64(ring.RandUniform(2, 3))
		sign1 = (tmp0 * -1) ^ (tmp0 ^ 1)
		sign2 = (tmp1 * -1) ^ (tmp1 ^ 1)
		coeffsD1[i] = int64(ring.RandUniform(500000, 0x7ffff)) * sign1
		coeffsD2[i] = int64(ring.RandUniform(500000, 0x7ffff)) * sign2
	}

	fmt.Printf("a : [")
	for i := uint64(0); i < 8; i++ {
		fmt.Printf("%8d", coeffsD1[i])
	}
	fmt.Printf("  ... ]\n")

	fmt.Printf("b : [")
	for i := uint64(0); i < 8; i++ {
		fmt.Printf("%8d", coeffsD2[i])
	}
	fmt.Printf("  ... ]\n")
	fmt.Println()

	// Computes the sum of the inner product
	sum := int64(0)
	innerProduct := make([]int64, N)
	for i := uint64(0); i < 784; i++ {
		innerProduct[i] = coeffsD1[i] * coeffsD2[i]
		sum += innerProduct[i]
	}

	// Creates two new plaintexts and encode the coeffs array on those plaintexts
	mD1 := bfvContext.NewPlaintext()
	if err := encoder.EncodeInt(coeffsD1, mD1); err != nil {
		log.Fatal(err)
	}

	mD2 := bfvContext.NewPlaintext()
	if err := encoder.EncodeInt(coeffsD2, mD2); err != nil {
		log.Fatal(err)
	}

	// Encrypts the first plaintext
	CtD1, err := Encryptor.EncryptNew(mD1)
	if err != nil {
		log.Fatal(err)
	}

	// Encrypts the second plaintext
	CtD2, err := Encryptor.EncryptNew(mD2)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Computation of the inner product followed by an inner sum...")
	// Inner product
	CtD1 = Evaluator.MulNew(CtD1, CtD2).(*bfv.Ciphertext)

	CtD1, err = Evaluator.RelinearizeNew(CtD1, rlk)
	if err != nil {
		log.Fatal(err)
	}

	// Inner Sum
	if err := Evaluator.InnerSum(CtD1, rotationKey, CtD1); err != nil {
		log.Fatal(err)
	}
	fmt.Println("Done!")
	fmt.Println()

	// Decrypts
	mR, err := Decryptor.DecryptNew(CtD1)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("Decrypt and decode : [")

	coeffs, err := encoder.DecodeInt(mR)
	if err != nil {
		log.Fatal(err)
	}

	for i := uint64(0); i < 4; i++ {
		fmt.Printf("%16d", coeffs[i])
	}
	fmt.Printf("  ... ]\n")
	fmt.Println("Inner sum of inner product :", coeffs[0], "==", sum, coeffs[0] == sum)

}

func main() {
	Plaintext_Batching()
	//Homomorphic_Inner_product()
}
