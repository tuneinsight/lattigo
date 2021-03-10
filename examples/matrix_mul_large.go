package main


import (
	"fmt"
	"time"
	"math/rand"
	"math/bits"
	"unsafe"
	"math"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/ckks"
)


type uint128 struct{
	hi uint64
	lo uint64
}


func MulCoeffsAndAdd128(a, b []uint64, c []uint128){

	var hi, lo, carry uint64

	for j := 0; j < len(a); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&a[j]))
		y := (*[8]uint64)(unsafe.Pointer(&b[j]))
		z := (*[8]uint128)(unsafe.Pointer(&c[j]))

		hi, lo = bits.Mul64(x[0], y[0])
		z[0].lo, carry = bits.Add64(z[0].lo, lo, 0)
		z[0].hi += hi + carry

		hi, lo = bits.Mul64(x[1], y[1])
		z[1].lo, carry = bits.Add64(z[1].lo, lo, 0)
		z[1].hi += hi + carry

		hi, lo = bits.Mul64(x[2], y[2])
		z[2].lo, carry = bits.Add64(z[2].lo, lo, 0)
		z[2].hi += hi + carry

		hi, lo = bits.Mul64(x[3], y[3])
		z[3].lo, carry = bits.Add64(z[3].lo, lo, 0)
		z[3].hi += hi + carry

		hi, lo = bits.Mul64(x[4], y[4])
		z[4].lo, carry = bits.Add64(z[4].lo, lo, 0)
		z[4].hi += hi + carry

		hi, lo = bits.Mul64(x[5], y[5])
		z[5].lo, carry = bits.Add64(z[5].lo, lo, 0)
		z[5].hi += hi + carry

		hi, lo = bits.Mul64(x[6], y[6])
		z[6].lo, carry = bits.Add64(z[6].lo, lo, 0)
		z[6].hi += hi + carry

		hi, lo = bits.Mul64(x[7], y[7])
		z[7].lo, carry = bits.Add64(z[7].lo, lo, 0)
		z[7].hi += hi + carry
	}
}


func ReduceAndAddUint128(in []uint128, out []uint64, qInv, q uint64){

	var hhi uint64

	for j := 0; j < len(in); j = j + 8 {

		x := (*[8]uint128)(unsafe.Pointer(&in[j]))
		y := (*[8]uint64)(unsafe.Pointer(&out[j]))

		hhi, _ = bits.Mul64(x[0].lo * qInv, q)
		y[0] += x[0].hi - hhi + q

		hhi, _ = bits.Mul64(x[1].lo*qInv, q)
		y[1] += x[1].hi - hhi + q 

		hhi, _ = bits.Mul64(x[2].lo*qInv, q)
		y[2] += x[2].hi - hhi + q 

		hhi, _ = bits.Mul64(x[3].lo*qInv, q)
		y[3] += x[3].hi - hhi + q 

		hhi, _ = bits.Mul64(x[4].lo*qInv, q)
		y[4] += x[4].hi - hhi + q 

		hhi, _ = bits.Mul64(x[5].lo*qInv, q)
		y[5] += x[5].hi - hhi + q 

		hhi, _ = bits.Mul64(x[6].lo*qInv, q)
		y[6] += x[6].hi - hhi + q 

		hhi, _ = bits.Mul64(x[7].lo*qInv, q)
		y[7] += x[7].hi - hhi + q 
	}
}


func main(){

	rand.Seed(1)

	LogN := uint64(12)
	LogSlots := LogN-1

	LogModuli := ckks.LogModuli{
		LogQi: []uint64{50, 35, 35},
		LogPi: []uint64{55},
	}

	Scale := float64(1 << 35)

	params, err := ckks.NewParametersFromLogModuli(LogN, &LogModuli)
	if err != nil {
		panic(err)
	}
	params.SetScale(Scale)
	params.SetLogSlots(LogSlots)

	// Set if the matrix is in plaintext or ciphertext (asymptotically twice the complexity)
	matrixInPlaintext := true

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptorFromSk(params, sk)
	decryptor := ckks.NewDecryptor(params, sk)
	rlk := kgen.GenRelinearizationKey(sk)

	rotations := []int{}

	for i := 1; i < int(params.Slots()); i<<=1{
		rotations = append(rotations, i)
	}

	rotkey := kgen.GenRotationKeysForRotations(rotations, false, sk)

	evaluator := ckks.NewEvaluator(params, ckks.EvaluationKey{rlk, rotkey})
	ringQ, _ := ring.NewRing(params.N(), params.Qi())


	values := make([]complex128, params.Slots())
	for j := range values{
		values[j] = complex((rand.Float64()-0.5)*2, 0)
	}

	A := encoder.EncodeNTTNew(values, params.LogSlots())
	ctA := encryptor.EncryptNew(A)

	for j := range values{
		values[j] = complex((rand.Float64()-0.5)*2, 0)
	}

	// Matrix in plaintext
	B := encoder.EncodeNTTNew(values, params.LogSlots())
	C := encryptor.EncryptNew(B)


	/*
	// Plaintext matrix must be in the Montgomery domain
	ringQ.MFormLvl(B.Level(), B.Value()[0], B.Value()[0])

	// Ciphertext matrix must be in the Montgomery domain
	ringQ.MFormLvl(C.Level(), C.Value()[0], C.Value()[0])
	ringQ.MFormLvl(C.Level(), C.Value()[1], C.Value()[1])
	*/


	// Level at which the plaintext/ciphertext are multiplied to the ciphetext
	ctLevelEval := params.MaxLevel() - 1


	// Multiplication of (s x n) x (n x m)

	//  ________________________     ____________________________________________________________________ ... _______
	// |						|   |																				 |
	// |		 s x n			| x |																				 |
	// |________________________|	|																				 |
	//								|																				 |
	//								|									n x m										 |
	//								|																				 |
	//								|																				 |
	//								|																				 |
	//								|																				 |
	//								|_____________________________________________________________________ ... ______|

	var s uint64 = 1
	var n uint64 = 1000
	var m uint64 = 1000000

	m_over_slots := uint64(math.Ceil(float64(m)/float64(params.Slots())))

	// Accumulators
	acc0 := make([][][]uint128, m_over_slots)
	acc1 := make([][][]uint128, m_over_slots)
	acc2 := make([][][]uint128, m_over_slots)
	for i := range acc0{
		acc0[i] = make([][]uint128, ctLevelEval+1) 
		acc1[i] = make([][]uint128, ctLevelEval+1) 
		acc2[i] = make([][]uint128, ctLevelEval+1) 

		for j := range acc0[i]{
			acc0[i][j] = make([]uint128, params.N()) // c0 * pt (plaintext matrix) or c0 * c0' (ciphertext matrix)
			acc1[i][j] = make([]uint128, params.N()) // c1 * pt (plaintext matrix) or c0 * c1' + c1 * c0' (ciphertext matrix)
			acc2[i][j] = make([]uint128, params.N()) // c1 * c1' (only if ciphertext x ciphertext)
		}
	}

	fmt.Printf("logN : %d\n", params.LogN())
	fmt.Printf("s : %d\n", s)
	fmt.Printf("n : %d\n", n)
	fmt.Printf("m : %d\n", m)
	fmt.Printf("nb values (s*n): %d\n", s*n)
	fmt.Printf("nb rotations (logN * s * n): %d\n", params.LogN()*n*s)
	fmt.Printf("nb mulpt&add ((m/N) * n): %d\n", m_over_slots*s*n)


	mask := make([]complex128, params.Slots()) 
	maskPt := ckks.NewPlaintext(params, ctLevelEval+1, float64(params.Qi()[ctA.Level()]))

	start := time.Now()

	tmp2 := ckks.NewCiphertext(params, 1, 2, params.Scale())
	// Computes

	ringQ.MFormLvl(ctA.Level(), ctA.Value()[0], ctA.Value()[0])
	ringQ.MFormLvl(ctA.Level(), ctA.Value()[1], ctA.Value()[1])

	for row := uint64(0); row < s; row++{
		for i := uint64(0); i < n; i++ {

			// Masking of one value [a, b, c, d] * [0, 1, 0, 0] = [0, b, 0, 0]
			tmp := ckks.NewCiphertext(params, 1, ctA.Level(), params.Scale())

			// TODO : precompute all the masks (altho, its fast, less than 1% of the total computation)
			mask[i] = complex(1, 0) // its a 1 in the mask at the position n
			encoder.EncodeNTT(maskPt, mask, params.LogSlots()) // encodes the mask

			// Multiplies with the mask
			evaluator.MulRelin(ctA, maskPt, tmp) 

			// Replicates the value in a ciphertext [0, b, 0, 0] -> [b, b, b, b]
			for i := uint64(0); i < params.LogSlots(); i++ {
				evaluator.Rotate(tmp, 1<<i, tmp2)
				evaluator.Add(tmp, tmp2, tmp)
			}

			// Rescales after to reduce the noisegrowth
			evaluator.Rescale(tmp, params.Scale(), tmp) 

			//ringQ.MFormLvl(tmp.Level(), tmp.Value()[0], tmp.Value()[0])
			//ringQ.MFormLvl(tmp.Level(), tmp.Value()[1], tmp.Value()[1])

			for j := uint64(0); j < m_over_slots; j++ {
				// Point-wise product with the extracted value [b, b, b, b] and the row encoded matrix m/N times
				// and addition to the accumulator
				if matrixInPlaintext {
					for x := uint64(0); x < ctLevelEval+1; x++{
						MulCoeffsAndAdd128(tmp.Value()[0].Coeffs[x], B.Value()[0].Coeffs[x], acc0[j][x]) // c0 * pt
						MulCoeffsAndAdd128(tmp.Value()[1].Coeffs[x], B.Value()[0].Coeffs[x], acc1[j][x]) // c1 * pt
					}
				}else{
					for x := uint64(0); x < ctLevelEval+1; x++{
						MulCoeffsAndAdd128(tmp.Value()[0].Coeffs[x], C.Value()[0].Coeffs[x], acc0[j][x]) // c0 * c0'
						MulCoeffsAndAdd128(tmp.Value()[0].Coeffs[x], C.Value()[1].Coeffs[x], acc1[j][x]) // c0 * c1'
						MulCoeffsAndAdd128(tmp.Value()[1].Coeffs[x], C.Value()[0].Coeffs[x], acc1[j][x]) // c1 * c0'
						MulCoeffsAndAdd128(tmp.Value()[1].Coeffs[x], C.Value()[1].Coeffs[x], acc2[j][x]) // c1 * c1'
					}
				}
			}
		}

		
		for j := uint64(0); j < m_over_slots; j++ {

			var ctRes *ckks.Ciphertext

			if matrixInPlaintext {

				ctRes = ckks.NewCiphertext(params, 1, 1, A.Scale() * B.Scale())

				for x := uint64(0); x < ctLevelEval+1; x++{
					mredParams := ringQ.MredParams[x]
					qi := ringQ.Modulus[x]
					ReduceAndAddUint128(acc0[j][x], ctRes.Value()[0].Coeffs[x], mredParams, qi)
					ReduceAndAddUint128(acc1[j][x], ctRes.Value()[1].Coeffs[x], mredParams, qi)
				}
			}else{
				ctRes = ckks.NewCiphertext(params, 2, 1, A.Scale() * B.Scale())

				for x := uint64(0); x < ctLevelEval+1; x++{
					mredParams := ringQ.MredParams[x]
					qi := ringQ.Modulus[x]
					ReduceAndAddUint128(acc0[j][x], ctRes.Value()[0].Coeffs[x], mredParams, qi)
					ReduceAndAddUint128(acc1[j][x], ctRes.Value()[1].Coeffs[x], mredParams, qi)
					ReduceAndAddUint128(acc2[j][x], ctRes.Value()[2].Coeffs[x], mredParams, qi)
				}

				evaluator.Relinearize(ctRes, ctRes)
			}
			
			valuesTest := encoder.Decode(decryptor.DecryptNew(ctRes), params.LogSlots())

			fmt.Println(valuesTest[:4])

			var error float64
			for i := range valuesTest{
				error += math.Abs(imag(valuesTest[i])) // should be close to zero
			}

			error /= float64(params.Slots())

			fmt.Printf("Values %6d - %6d : %f\n", j*params.Slots(), (j+1)*params.Slots(), error)

		}
	}

	fmt.Printf("Done in %s \n", time.Since(start))
}