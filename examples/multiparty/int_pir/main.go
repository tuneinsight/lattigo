package main

import (
	"log"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func check(err error) {
	if err != nil {
		panic(err)
	}
}

func runTimed(f func()) time.Duration {
	start := time.Now()
	f()
	return time.Since(start)
}

func runTimedParty(f func(), N int) time.Duration {
	start := time.Now()
	f()
	return time.Duration(time.Since(start).Nanoseconds() / int64(N))
}

type party struct {
	sk         *rlwe.SecretKey
	rlkEphemSk *rlwe.SecretKey

	ckgShare    multiparty.PublicKeyGenShare
	rkgShareOne multiparty.RelinearizationKeyGenShare
	rkgShareTwo multiparty.RelinearizationKeyGenShare
	gkgShare    multiparty.GaloisKeyGenShare
	cksShare    multiparty.KeySwitchShare

	input []uint64
}

type maskTask struct {
	query           *rlwe.Ciphertext
	mask            *rlwe.Plaintext
	row             *rlwe.Ciphertext
	res             *rlwe.Ciphertext
	elapsedmaskTask time.Duration
}

var elapsedCKGCloud time.Duration
var elapsedCKGParty time.Duration
var elapsedRKGCloud time.Duration
var elapsedRKGParty time.Duration
var elapsedGKGCloud time.Duration
var elapsedGKGParty time.Duration
var elapsedCKSCloud time.Duration
var elapsedPCKSParty time.Duration
var elapsedRequestParty time.Duration
var elapsedRequestCloud time.Duration
var elapsedRequestCloudCPU time.Duration

func main() {

	// This example simulates a SMC instance of a private information retrieval (PIR) problem.
	// The problem statement is as follows: a cloud stores data of several parties
	// encrypted under a shared public-key. An external party wants to retrieve
	// the plaintext content of one of the ciphertexts while ensuring the following
	// security property: no information other than the fact that a request was made must
	// be disclosed to the cloud, to the owners of the shared public-key or to anyone else.
	//
	// For more details see
	//    Multiparty Homomorphic Encryption: From Theory to Practice (<https://eprint.iacr.org/2020/304>)

	l := log.New(os.Stderr, "", 0)

	// $go run main.go arg1 arg2
	// arg1: number of parties
	// arg2: number of Go routines
	// MinDelta number of parties for n=8192: 512 parties (this is a memory intensive process)

	N := 3 // Default number of parties
	var err error
	if len(os.Args[1:]) >= 1 {
		N, err = strconv.Atoi(os.Args[1])
		check(err)
	}

	NGoRoutine := 1 // Default number of Go routines
	if len(os.Args[1:]) >= 2 {
		NGoRoutine, err = strconv.Atoi(os.Args[2])
		check(err)
	}

	// Index of the ciphertext to retrieve.
	queryIndex := 2

	// Creating encryption parameters
	// LogN = 13 & LogQP = 218
	params, err := bgv.NewParametersFromLiteral(bgv.ParametersLiteral{
		LogN:             13,
		LogQ:             []int{54, 54, 54},
		LogP:             []int{55},
		PlaintextModulus: 65537,
	})
	if err != nil {
		panic(err)
	}

	// Common reference polynomial generator that uses the PRNG
	crs, err := sampling.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}

	// Instantiation of each of the protocols needed for the PIR example

	// Create each party, and allocate the memory for all the shares that the protocols will need
	P := genparties(params, N)

	// 1) Collective public key generation
	pk := ckgphase(params, crs, P)

	// 2) Collective RelinearizationKey generation
	RelinearizationKey := rkgphase(params, crs, P)

	// 3) Collective GaloisKeys generation
	galKeys := gkgphase(params, crs, P)

	// Instantiates EvaluationKeySet
	evk := rlwe.NewMemEvaluationKeySet(RelinearizationKey, galKeys...)

	l.Printf("\tSetup done (cloud: %s, party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedGKGCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedGKGParty)

	// Pre-loading memory
	encoder := bgv.NewEncoder(params)
	l.Println("> Memory alloc Phase")
	encInputs := make([]*rlwe.Ciphertext, N)
	plainMask := make([]*rlwe.Plaintext, N)

	// Ciphertexts to be retrieved
	for i := range encInputs {
		encInputs[i] = bgv.NewCiphertext(params, 1, params.MaxLevel())
	}

	// Plaintext masks: plainmask[i] = encode([0, ..., 0, 1_i, 0, ..., 0])
	// (zero with a 1 at the i-th position).
	for i := range plainMask {
		maskCoeffs := make([]uint64, params.N())
		maskCoeffs[i] = 1
		plainMask[i] = bgv.NewPlaintext(params, params.MaxLevel())
		if err := encoder.Encode(maskCoeffs, plainMask[i]); err != nil {
			panic(err)
		}
	}

	// Ciphertexts encrypted under collective public key and stored in the cloud
	l.Println("> Encrypt Phase")
	encryptor := rlwe.NewEncryptor(params, pk)
	pt := bgv.NewPlaintext(params, params.MaxLevel())
	elapsedEncryptParty := runTimedParty(func() {
		for i, pi := range P {
			if err := encoder.Encode(pi.input, pt); err != nil {
				panic(err)
			}
			if err := encryptor.Encrypt(pt, encInputs[i]); err != nil {
				panic(err)
			}
		}
	}, N)

	elapsedEncryptCloud := time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	// Request phase
	encQuery := genquery(params, queryIndex, encoder, encryptor)

	result := requestphase(params, queryIndex, NGoRoutine, encQuery, encInputs, plainMask, evk)

	// Collective (partial) decryption (key switch)
	encOut := cksphase(params, P, result)

	l.Println("> ResulPlaintextModulus:")

	// Decryption by the external party
	decryptor := rlwe.NewDecryptor(params, P[0].sk)
	ptres := bgv.NewPlaintext(params, params.MaxLevel())
	elapsedDecParty := runTimed(func() {
		decryptor.Decrypt(encOut, ptres)
	})

	res := make([]uint64, params.MaxSlots())
	if err := encoder.Decode(ptres, res); err != nil {
		panic(err)
	}

	l.Printf("\t%v...%v\n", res[:8], res[params.N()-8:])
	l.Printf("> Finished (total cloud: %s, total party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedGKGCloud+elapsedEncryptCloud+elapsedRequestCloudCPU+elapsedCKSCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedGKGParty+elapsedEncryptParty+elapsedRequestParty+elapsedPCKSParty+elapsedDecParty)
}

func cksphase(params bgv.Parameters, P []*party, result *rlwe.Ciphertext) *rlwe.Ciphertext {
	l := log.New(os.Stderr, "", 0)

	l.Println("> KeySwitch Phase")

	cks, err := multiparty.NewKeySwitchProtocol(params, ring.DiscreteGaussian{Sigma: 1 << 30, Bound: 6 * (1 << 30)}) // Collective public-key re-encryption
	if err != nil {
		panic(err)
	}

	for _, pi := range P {
		pi.cksShare = cks.AllocateShare(params.MaxLevel())
	}

	zero := rlwe.NewSecretKey(params)
	cksCombined := cks.AllocateShare(params.MaxLevel())
	elapsedPCKSParty = runTimedParty(func() {
		for _, pi := range P[1:] {
			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			cks.GenShare(pi.sk, zero, result, &pi.cksShare)
		}
	}, len(P)-1)

	encOut := bgv.NewCiphertext(params, 1, params.MaxLevel())
	elapsedCKSCloud = runTimed(func() {
		for _, pi := range P {
			if err := cks.AggregateShares(pi.cksShare, cksCombined, &cksCombined); err != nil {
				panic(err)
			}
		}
		cks.KeySwitch(result, cksCombined, encOut)
	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKSCloud, elapsedPCKSParty)

	return encOut
}

func genparties(params bgv.Parameters, N int) []*party {

	P := make([]*party, N)

	kgen := rlwe.NewKeyGenerator(params)

	for i := range P {
		pi := &party{}
		pi.sk = kgen.GenSecretKeyNew()

		pi.input = make([]uint64, params.N())
		for j := range pi.input {
			pi.input[j] = uint64(i)
		}

		P[i] = pi
	}

	return P
}

func ckgphase(params bgv.Parameters, crs sampling.PRNG, P []*party) *rlwe.PublicKey {

	l := log.New(os.Stderr, "", 0)

	l.Println("> PublicKeyGen Phase")

	ckg := multiparty.NewPublicKeyGenProtocol(params) // Public key generation

	ckgCombined := ckg.AllocateShare()
	for _, pi := range P {
		pi.ckgShare = ckg.AllocateShare()
	}

	crp := ckg.SampleCRP(crs)

	elapsedCKGParty = runTimedParty(func() {
		for _, pi := range P {
			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			ckg.GenShare(pi.sk, crp, &pi.ckgShare)
		}
	}, len(P))

	pk := rlwe.NewPublicKey(params)

	elapsedCKGCloud = runTimed(func() {
		for _, pi := range P {
			ckg.AggregateShares(pi.ckgShare, ckgCombined, &ckgCombined)
		}
		ckg.GenPublicKey(ckgCombined, crp, pk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKGCloud, elapsedCKGParty)

	return pk
}

func rkgphase(params bgv.Parameters, crs sampling.PRNG, P []*party) *rlwe.RelinearizationKey {
	l := log.New(os.Stderr, "", 0)

	l.Println("> RelinearizationKeyGen Phase")

	rkg := multiparty.NewRelinearizationKeyGenProtocol(params) // Relineariation key generation

	_, rkgCombined1, rkgCombined2 := rkg.AllocateShare()

	for _, pi := range P {
		pi.rlkEphemSk, pi.rkgShareOne, pi.rkgShareTwo = rkg.AllocateShare()
	}

	crp := rkg.SampleCRP(crs)

	elapsedRKGParty = runTimedParty(func() {
		for _, pi := range P {
			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			rkg.GenShareRoundOne(pi.sk, crp, pi.rlkEphemSk, &pi.rkgShareOne)
		}
	}, len(P))

	elapsedRKGCloud = runTimed(func() {
		for _, pi := range P {
			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			rkg.AggregateShares(pi.rkgShareOne, rkgCombined1, &rkgCombined1)
		}
	})

	elapsedRKGParty += runTimedParty(func() {
		for _, pi := range P {
			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			rkg.GenShareRoundTwo(pi.rlkEphemSk, pi.sk, rkgCombined1, &pi.rkgShareTwo)
		}
	}, len(P))

	rlk := rlwe.NewRelinearizationKey(params)
	elapsedRKGCloud += runTimed(func() {
		for _, pi := range P {
			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			rkg.AggregateShares(pi.rkgShareTwo, rkgCombined2, &rkgCombined2)
		}
		rkg.GenRelinearizationKey(rkgCombined1, rkgCombined2, rlk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedRKGCloud, elapsedRKGParty)

	return rlk
}

func gkgphase(params bgv.Parameters, crs sampling.PRNG, P []*party) (galKeys []*rlwe.GaloisKey) {

	l := log.New(os.Stderr, "", 0)

	l.Println("> RTG Phase")

	gkg := multiparty.NewGaloisKeyGenProtocol(params) // Rotation keys generation

	for _, pi := range P {
		pi.gkgShare = gkg.AllocateShare()
	}

	galEls := append(params.GaloisElementsForInnerSum(1, params.N()>>1), params.GaloisElementForRowRotation())
	galKeys = make([]*rlwe.GaloisKey, len(galEls))

	gkgShareCombined := gkg.AllocateShare()

	for i, galEl := range galEls {

		gkgShareCombined.GaloisElement = galEl

		crp := gkg.SampleCRP(crs)

		elapsedGKGParty += runTimedParty(func() {
			for _, pi := range P {
				/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
				if err := gkg.GenShare(pi.sk, galEl, crp, &pi.gkgShare); err != nil {
					panic(err)
				}
			}

		}, len(P))

		elapsedGKGCloud += runTimed(func() {

			if err := gkg.AggregateShares(P[0].gkgShare, P[1].gkgShare, &gkgShareCombined); err != nil {
				panic(err)
			}

			for _, pi := range P[2:] {
				if err := gkg.AggregateShares(pi.gkgShare, gkgShareCombined, &gkgShareCombined); err != nil {
					panic(err)
				}
			}

			galKeys[i] = rlwe.NewGaloisKey(params)

			if err := gkg.GenGaloisKey(gkgShareCombined, crp, galKeys[i]); err != nil {
				panic(err)
			}
		})
	}
	l.Printf("\tdone (cloud: %s, party %s)\n", elapsedGKGCloud, elapsedGKGParty)

	return
}

func genquery(params bgv.Parameters, queryIndex int, encoder *bgv.Encoder, encryptor *rlwe.Encryptor) *rlwe.Ciphertext {
	// Query ciphertext
	queryCoeffs := make([]uint64, params.N())
	queryCoeffs[queryIndex] = 1
	query := bgv.NewPlaintext(params, params.MaxLevel())
	var encQuery *rlwe.Ciphertext
	elapsedRequestParty += runTimed(func() {
		var err error
		if err = encoder.Encode(queryCoeffs, query); err != nil {
			panic(err)
		}
		if encQuery, err = encryptor.EncryptNew(query); err != nil {
			panic(err)
		}
	})

	return encQuery
}

func requestphase(params bgv.Parameters, queryIndex, NGoRoutine int, encQuery *rlwe.Ciphertext, encInputs []*rlwe.Ciphertext, plainMask []*rlwe.Plaintext, evk rlwe.EvaluationKeySet) *rlwe.Ciphertext {

	l := log.New(os.Stderr, "", 0)

	l.Println("> Request Phase")

	// Buffer for the intermediate computation done by the cloud
	encPartial := make([]*rlwe.Ciphertext, len(encInputs))
	for i := range encPartial {
		encPartial[i] = bgv.NewCiphertext(params, 2, params.MaxLevel())
	}

	evaluator := bgv.NewEvaluator(params, evk)

	// Split the task among the Go routines
	tasks := make(chan *maskTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {
			evaluator := evaluator.ShallowCopy() // creates a shallow evaluator copy for this goroutine
			tmp := bgv.NewCiphertext(params, 1, params.MaxLevel())
			for task := range tasks {
				task.elapsedmaskTask = runTimed(func() {
					// 1) Multiplication BFV-style of the query with the plaintext mask
					if err := evaluator.MulScaleInvariant(task.query, task.mask, tmp); err != nil {
						panic(err)
					}

					// 2) Inner sum (populate all the slots with the sum of all the slots)
					if err := evaluator.InnerSum(tmp, 1, params.N()>>1, tmp); err != nil {
						panic(err)
					}

					if tmpRot, err := evaluator.RotateRowsNew(tmp); err != nil {

					} else {
						if err := evaluator.Add(tmp, tmpRot, tmp); err != nil {
							panic(err)
						}
					}

					// 3) Multiplication of 2) with the i-th ciphertext stored in the cloud
					if err := evaluator.Mul(tmp, task.row, task.res); err != nil {
						panic(err)
					}
				})
			}
			//l.Println("\t evaluator", i, "down")
			workers.Done()
		}(i)
		//l.Println("\t evaluator", i, "started")
	}

	taskList := make([]*maskTask, 0)

	elapsedRequestCloud += runTimed(func() {
		for i := range encInputs {
			task := &maskTask{
				query: encQuery,
				mask:  plainMask[i],
				row:   encInputs[i],
				res:   encPartial[i],
			}
			taskList = append(taskList, task)
			tasks <- task
		}
		close(tasks)
		workers.Wait()
	})

	for _, t := range taskList {
		elapsedRequestCloudCPU += t.elapsedmaskTask
	}

	resultDeg2 := bgv.NewCiphertext(params, 2, params.MaxLevel())
	result := bgv.NewCiphertext(params, 1, params.MaxLevel())

	// Summation of all the partial result among the different Go routines
	finalAddDuration := runTimed(func() {
		for i := 0; i < len(encInputs); i++ {
			if err := evaluator.Add(resultDeg2, encPartial[i], resultDeg2); err != nil {
				panic(err)
			}
		}
		if err := evaluator.Relinearize(resultDeg2, result); err != nil {
			panic(err)
		}
	})

	elapsedRequestCloud += finalAddDuration
	elapsedRequestCloudCPU += finalAddDuration

	l.Printf("\tdone (cloud: %s/%s, party: %s)\n",
		elapsedRequestCloud, elapsedRequestCloudCPU, elapsedRequestParty)

	return result
}
