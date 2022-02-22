package main

import (
	"log"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/dckks"
	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
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

	ckgShare    *drlwe.CKGShare
	rkgShareOne *drlwe.RKGShare
	rkgShareTwo *drlwe.RKGShare
	rtgShare    *drlwe.RTGShare
	cksShare    *drlwe.CKSShare

	input []complex128
}

type maskTask struct {
	query           *ckks.Ciphertext
	mask            *ckks.Plaintext
	row             *ckks.Ciphertext
	res             *ckks.Ciphertext
	elapsedmaskTask time.Duration
}

var elapsedCKGCloud time.Duration
var elapsedCKGParty time.Duration
var elapsedRKGCloud time.Duration
var elapsedRKGParty time.Duration
var elapsedRTGCloud time.Duration
var elapsedRTGParty time.Duration
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

	// Creating encryption parameters from a default params with N=8192
	paramsDef := ckks.PN13QP218
	params, err := ckks.NewParametersFromLiteral(paramsDef)
	if err != nil {
		panic(err)
	}

	// Common reference polynomial generator that uses the PRNG
	crs, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}

	// Instantiation of each of the protocols needed for the PIR example

	// Create each party, and allocate the memory for all the shares that the protocols will need
	P := genparties(params, N)

	// 1) Collective public key generation
	pk := ckgphase(params, crs, P)

	// 2) Collective relinearization key generation
	rlk := rkgphase(params, crs, P)

	// 3) Collective rotation keys generation
	rtk := rtkphase(params, crs, P)

	l.Printf("\tSetup done (cloud: %s, party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedRTGCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedRTGParty)

	// Pre-loading memory
	encoder := ckks.NewEncoder(params)
	l.Println("> Memory alloc Phase")
	encInputs := make([]*ckks.Ciphertext, params.Slots())
	plainMask := make([]*ckks.Plaintext, params.Slots())

	// Ciphertexts to be retrieved
	for i := range encInputs {
		encInputs[i] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())
	}

	// Plaintext masks: plainmask[i] = encode([0, ..., 0, 1_i, 0, ..., 0])
	// (zero with a 1 at the i-th position).
	for i := range plainMask {
		maskCoeffs := make([]float64, params.Slots())
		maskCoeffs[i] = 1
		plainMask[i] = ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
		encoder.EncodeCoeffs(maskCoeffs, plainMask[i])
	}

	// Ciphertexts encrypted under CKG and stored in the cloud
	l.Println("> Encrypt Phase")
	encryptor := ckks.NewEncryptor(params, pk)
	pt := ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
	elapsedEncryptParty := runTimedParty(func() {
		for i, pi := range P {
			encoder.Encode(pi.input, pt, params.LogSlots())
			encryptor.Encrypt(pt, encInputs[i])
		}
	}, N)

	elapsedEncryptCloud := time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	// Request phase
	encQuery := genquery(params, queryIndex, encoder, encryptor)

	result := requestphase(params, queryIndex, NGoRoutine, encQuery, encInputs, plainMask, rlk, rtk)

	// Collective (partial) decryption (key switch)
	encOut := cksphase(params, P, result)

	l.Println("> Result:")

	// Decryption by the external party
	decryptor := ckks.NewDecryptor(params, P[0].sk)
	ptres := ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
	elapsedDecParty := runTimed(func() {
		decryptor.Decrypt(encOut, ptres)
	})

	res := encoder.Decode(ptres, params.LogSlots())

	l.Printf("\t%v\n", res[:16])
	l.Printf("> Finished (total cloud: %s, total party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedRTGCloud+elapsedEncryptCloud+elapsedRequestCloudCPU+elapsedCKSCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedRTGParty+elapsedEncryptParty+elapsedRequestParty+elapsedPCKSParty+elapsedDecParty)
}

func cksphase(params ckks.Parameters, P []*party, result *ckks.Ciphertext) *ckks.Ciphertext {
	l := log.New(os.Stderr, "", 0)

	l.Println("> CKS Phase")

	cks := dckks.NewCKSProtocol(params, 3.19) // Collective public-key re-encryption

	for _, pi := range P {
		pi.cksShare = cks.AllocateShare(params.MaxLevel())
	}

	zero := ckks.NewSecretKey(params)
	cksCombined := cks.AllocateShare(params.MaxLevel())
	elapsedPCKSParty = runTimedParty(func() {
		for _, pi := range P[1:] {
			cks.GenShare(pi.sk, zero, result.Value[1], pi.cksShare)
		}
	}, len(P)-1)

	encOut := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())
	elapsedCKSCloud = runTimed(func() {
		for _, pi := range P {
			cks.AggregateShare(pi.cksShare, cksCombined, cksCombined)
		}
		cks.KeySwitch(result, cksCombined, encOut)
	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKSCloud, elapsedPCKSParty)

	return encOut
}

func genparties(params ckks.Parameters, N int) []*party {

	P := make([]*party, N)

	kgen := ckks.NewKeyGenerator(params)

	for i := range P {
		pi := &party{}
		pi.sk = kgen.GenSecretKey()

		pi.input = make([]complex128, params.Slots())
		for j := range pi.input {
			pi.input[j] = complex(float64(i), 0)
		}

		P[i] = pi
	}

	return P
}

func ckgphase(params ckks.Parameters, crs utils.PRNG, P []*party) *rlwe.PublicKey {

	l := log.New(os.Stderr, "", 0)

	l.Println("> CKG Phase")

	ckg := dckks.NewCKGProtocol(params) // Public key generation

	ckgCombined := ckg.AllocateShare()
	for _, pi := range P {
		pi.ckgShare = ckg.AllocateShare()
	}

	crp := ckg.SampleCRP(crs)

	elapsedCKGParty = runTimedParty(func() {
		for _, pi := range P {
			ckg.GenShare(pi.sk, crp, pi.ckgShare)
		}
	}, len(P))

	pk := ckks.NewPublicKey(params)

	elapsedCKGCloud = runTimed(func() {
		for _, pi := range P {
			ckg.AggregateShare(pi.ckgShare, ckgCombined, ckgCombined)
		}
		ckg.GenPublicKey(ckgCombined, crp, pk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKGCloud, elapsedCKGParty)

	return pk
}

func rkgphase(params ckks.Parameters, crs utils.PRNG, P []*party) *rlwe.RelinearizationKey {
	l := log.New(os.Stderr, "", 0)

	l.Println("> RKG Phase")

	rkg := dckks.NewRKGProtocol(params) // Relineariation key generation

	_, rkgCombined1, rkgCombined2 := rkg.AllocateShare()

	for _, pi := range P {
		pi.rlkEphemSk, pi.rkgShareOne, pi.rkgShareTwo = rkg.AllocateShare()
	}

	crp := rkg.SampleCRP(crs)

	elapsedRKGParty = runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundOne(pi.sk, crp, pi.rlkEphemSk, pi.rkgShareOne)
		}
	}, len(P))

	elapsedRKGCloud = runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShare(pi.rkgShareOne, rkgCombined1, rkgCombined1)
		}
	})

	elapsedRKGParty += runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundTwo(pi.rlkEphemSk, pi.sk, rkgCombined1, pi.rkgShareTwo)
		}
	}, len(P))

	rlk := ckks.NewRelinearizationKey(params)
	elapsedRKGCloud += runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShare(pi.rkgShareTwo, rkgCombined2, rkgCombined2)
		}
		rkg.GenRelinearizationKey(rkgCombined1, rkgCombined2, rlk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedRKGCloud, elapsedRKGParty)

	return rlk
}

func rtkphase(params ckks.Parameters, crs utils.PRNG, P []*party) *rlwe.RotationKeySet {

	l := log.New(os.Stderr, "", 0)

	l.Println("> RTG Phase")

	rtg := dckks.NewRotKGProtocol(params) // Rotation keys generation

	for _, pi := range P {
		pi.rtgShare = rtg.AllocateShare()
	}

	galEls := params.GaloisElementsForRowInnerSum()
	rotKeySet := ckks.NewRotationKeySet(params, galEls)

	for _, galEl := range galEls {

		rtgShareCombined := rtg.AllocateShare()

		crp := rtg.SampleCRP(crs)

		elapsedRTGParty += runTimedParty(func() {
			for _, pi := range P {
				rtg.GenShare(pi.sk, galEl, crp, pi.rtgShare)
			}
		}, len(P))

		elapsedRTGCloud += runTimed(func() {
			for _, pi := range P {
				rtg.AggregateShare(pi.rtgShare, rtgShareCombined, rtgShareCombined)
			}
			rtg.GenRotationKey(rtgShareCombined, crp, rotKeySet.Keys[galEl])
		})
	}
	l.Printf("\tdone (cloud: %s, party %s)\n", elapsedRTGCloud, elapsedRTGParty)

	return rotKeySet
}

func genquery(params ckks.Parameters, queryIndex int, encoder ckks.Encoder, encryptor ckks.Encryptor) *ckks.Ciphertext {
	// Query ciphertext
	queryCoeffs := make([]float64, params.Slots())
	queryCoeffs[queryIndex] = 1
	query := ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
	var encQuery *ckks.Ciphertext
	elapsedRequestParty += runTimed(func() {
		encoder.EncodeCoeffs(queryCoeffs, query)
		encQuery = encryptor.EncryptNew(query)
	})

	return encQuery
}

func requestphase(params ckks.Parameters, queryIndex, NGoRoutine int, encQuery *ckks.Ciphertext, encInputs []*ckks.Ciphertext, plainMask []*ckks.Plaintext, rlk *rlwe.RelinearizationKey, rtk *rlwe.RotationKeySet) *ckks.Ciphertext {

	l := log.New(os.Stderr, "", 0)

	l.Println("> Request Phase")

	// Buffer for the intermediate computation done by the cloud
	encPartial := make([]*ckks.Ciphertext, len(encInputs))
	for i := range encPartial {
		encPartial[i] = ckks.NewCiphertext(params, 2, params.MaxLevel(), params.DefaultScale())
	}

	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rtk})

	// Split the task among the Go routines
	tasks := make(chan *maskTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {
			evaluator := evaluator.ShallowCopy() // creates a shallow evaluator copy for this goroutine
			tmp := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())
			for task := range tasks {
				task.elapsedmaskTask = runTimed(func() {
					// 1) Multiplication of the query with the plaintext mask
					evaluator.Mul(task.query, task.mask, tmp)
					// 2) Inner sum (populate all the slots with the sum of all the slots)
					evaluator.InnerSum(tmp, 1, 1, tmp)
					// 3) Multiplication of 2) with the i-th ciphertext stored in the cloud
					evaluator.Mul(tmp, task.row, task.res)
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

	resultDeg2 := ckks.NewCiphertext(params, 2, params.MaxLevel(), params.DefaultScale())
	result := ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())

	// Summation of all the partial result among the different Go routines
	finalAddDuration := runTimed(func() {
		for i := 0; i < len(encInputs); i++ {
			evaluator.Add(resultDeg2, encPartial[i], resultDeg2)
		}
		evaluator.Relinearize(resultDeg2, result)
	})

	elapsedRequestCloud += finalAddDuration
	elapsedRequestCloudCPU += finalAddDuration

	l.Printf("\tdone (cloud: %s/%s, party: %s)\n",
		elapsedRequestCloud, elapsedRequestCloudCPU, elapsedRequestParty)

	return result
}
