package main

import (
	"log"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/dbfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
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

	ckgShare    drlwe.PublicKeyGenShare
	rkgShareOne drlwe.RelinKeyGenShare
	rkgShareTwo drlwe.RelinKeyGenShare
	gkgShare    drlwe.GaloisKeyGenShare
	cksShare    drlwe.KeySwitchShare

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
	params, err := bfv.NewParametersFromLiteral(bfv.ParametersLiteral{
		LogN: 13,
		LogQ: []int{54, 54, 54},
		LogP: []int{55},
		T:    65537,
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
	relinKey := rkgphase(params, crs, P)

	// 3) Collective GaloisKeys generation
	galKeys := gkgphase(params, crs, P)

	// Instantiates EvaluationKeySet
	evk := rlwe.NewMemEvaluationKeySet(relinKey, galKeys...)

	l.Printf("\tSetup done (cloud: %s, party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedGKGCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedGKGParty)

	// Pre-loading memory
	encoder := bfv.NewEncoder(params)
	l.Println("> Memory alloc Phase")
	encInputs := make([]*rlwe.Ciphertext, N)
	plainMask := make([]*rlwe.Plaintext, N)

	// Ciphertexts to be retrieved
	for i := range encInputs {
		encInputs[i] = bfv.NewCiphertext(params, 1, params.MaxLevel())
	}

	// Plaintext masks: plainmask[i] = encode([0, ..., 0, 1_i, 0, ..., 0])
	// (zero with a 1 at the i-th position).
	for i := range plainMask {
		maskCoeffs := make([]uint64, params.N())
		maskCoeffs[i] = 1
		plainMask[i] = bfv.NewPlaintext(params, params.MaxLevel())
		if err := encoder.Encode(maskCoeffs, plainMask[i]); err != nil {
			panic(err)
		}
	}

	// Ciphertexts encrypted under collective public key and stored in the cloud
	l.Println("> Encrypt Phase")
	encryptor := bfv.NewEncryptor(params, pk)
	pt := bfv.NewPlaintext(params, params.MaxLevel())
	elapsedEncryptParty := runTimedParty(func() {
		for i, pi := range P {
			if err := encoder.Encode(pi.input, pt); err != nil {
				panic(err)
			}
			encryptor.Encrypt(pt, encInputs[i])
		}
	}, N)

	elapsedEncryptCloud := time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	// Request phase
	encQuery := genquery(params, queryIndex, encoder, encryptor)

	result := requestphase(params, queryIndex, NGoRoutine, encQuery, encInputs, plainMask, evk)

	// Collective (partial) decryption (key switch)
	encOut := cksphase(params, P, result)

	l.Println("> Result:")

	// Decryption by the external party
	decryptor := bfv.NewDecryptor(params, P[0].sk)
	ptres := bfv.NewPlaintext(params, params.MaxLevel())
	elapsedDecParty := runTimed(func() {
		decryptor.Decrypt(encOut, ptres)
	})

	res := make([]uint64, params.PlaintextSlots())
	if err := encoder.Decode(ptres, res); err != nil {
		panic(err)
	}

	l.Printf("\t%v...%v\n", res[:8], res[params.N()-8:])
	l.Printf("> Finished (total cloud: %s, total party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedGKGCloud+elapsedEncryptCloud+elapsedRequestCloudCPU+elapsedCKSCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedGKGParty+elapsedEncryptParty+elapsedRequestParty+elapsedPCKSParty+elapsedDecParty)
}

func cksphase(params bfv.Parameters, P []*party, result *rlwe.Ciphertext) *rlwe.Ciphertext {
	l := log.New(os.Stderr, "", 0)

	l.Println("> KeySwitch Phase")

	cks := dbfv.NewKeySwitchProtocol(params, ring.DiscreteGaussian{Sigma: 1 << 30, Bound: 6 * (1 << 30)}) // Collective public-key re-encryption

	for _, pi := range P {
		pi.cksShare = cks.AllocateShare(params.MaxLevel())
	}

	zero := rlwe.NewSecretKey(params.Parameters)
	cksCombined := cks.AllocateShare(params.MaxLevel())
	elapsedPCKSParty = runTimedParty(func() {
		for _, pi := range P[1:] {
			cks.GenShare(pi.sk, zero, result, &pi.cksShare)
		}
	}, len(P)-1)

	encOut := bfv.NewCiphertext(params, 1, params.MaxLevel())
	elapsedCKSCloud = runTimed(func() {
		for _, pi := range P {
			cks.AggregateShares(&pi.cksShare, &cksCombined, &cksCombined)
		}
		cks.KeySwitch(result, cksCombined, encOut)
	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKSCloud, elapsedPCKSParty)

	return encOut
}

func genparties(params bfv.Parameters, N int) []*party {

	P := make([]*party, N)

	kgen := bfv.NewKeyGenerator(params)

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

func ckgphase(params bfv.Parameters, crs sampling.PRNG, P []*party) *rlwe.PublicKey {

	l := log.New(os.Stderr, "", 0)

	l.Println("> PublicKeyGen Phase")

	ckg := dbfv.NewPublicKeyGenProtocol(params) // Public key generation

	ckgCombined := ckg.AllocateShare()
	for _, pi := range P {
		pi.ckgShare = ckg.AllocateShare()
	}

	crp := ckg.SampleCRP(crs)

	elapsedCKGParty = runTimedParty(func() {
		for _, pi := range P {
			ckg.GenShare(pi.sk, crp, &pi.ckgShare)
		}
	}, len(P))

	pk := rlwe.NewPublicKey(params.Parameters)

	elapsedCKGCloud = runTimed(func() {
		for _, pi := range P {
			ckg.AggregateShares(&pi.ckgShare, &ckgCombined, &ckgCombined)
		}
		ckg.GenPublicKey(ckgCombined, crp, pk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKGCloud, elapsedCKGParty)

	return pk
}

func rkgphase(params bfv.Parameters, crs sampling.PRNG, P []*party) *rlwe.RelinearizationKey {
	l := log.New(os.Stderr, "", 0)

	l.Println("> RelinKeyGen Phase")

	rkg := dbfv.NewRelinKeyGenProtocol(params) // Relineariation key generation

	_, rkgCombined1, rkgCombined2 := rkg.AllocateShare()

	for _, pi := range P {
		pi.rlkEphemSk, pi.rkgShareOne, pi.rkgShareTwo = rkg.AllocateShare()
	}

	crp := rkg.SampleCRP(crs)

	elapsedRKGParty = runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundOne(pi.sk, crp, pi.rlkEphemSk, &pi.rkgShareOne)
		}
	}, len(P))

	elapsedRKGCloud = runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShares(&pi.rkgShareOne, &rkgCombined1, &rkgCombined1)
		}
	})

	elapsedRKGParty += runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundTwo(pi.rlkEphemSk, pi.sk, rkgCombined1, &pi.rkgShareTwo)
		}
	}, len(P))

	rlk := rlwe.NewRelinearizationKey(params.Parameters)
	elapsedRKGCloud += runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShares(&pi.rkgShareTwo, &rkgCombined2, &rkgCombined2)
		}
		rkg.GenRelinearizationKey(rkgCombined1, rkgCombined2, rlk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedRKGCloud, elapsedRKGParty)

	return rlk
}

func gkgphase(params bfv.Parameters, crs sampling.PRNG, P []*party) (galKeys []*rlwe.GaloisKey) {

	l := log.New(os.Stderr, "", 0)

	l.Println("> RTG Phase")

	gkg := dbfv.NewGaloisKeyGenProtocol(params) // Rotation keys generation

	for _, pi := range P {
		pi.gkgShare = gkg.AllocateShare()
	}

	galEls := append(params.GaloisElementsForInnerSum(1, params.N()>>1), params.GaloisElementInverse())
	galKeys = make([]*rlwe.GaloisKey, len(galEls))

	gkgShareCombined := gkg.AllocateShare()

	for i, galEl := range galEls {

		gkgShareCombined.GaloisElement = galEl

		crp := gkg.SampleCRP(crs)

		elapsedGKGParty += runTimedParty(func() {
			for _, pi := range P {
				gkg.GenShare(pi.sk, galEl, crp, &pi.gkgShare)
			}

		}, len(P))

		elapsedGKGCloud += runTimed(func() {

			gkg.AggregateShares(&P[0].gkgShare, &P[1].gkgShare, &gkgShareCombined)

			for _, pi := range P[2:] {
				gkg.AggregateShares(&pi.gkgShare, &gkgShareCombined, &gkgShareCombined)
			}

			galKeys[i] = rlwe.NewGaloisKey(params.Parameters)

			gkg.GenGaloisKey(gkgShareCombined, crp, galKeys[i])
		})
	}
	l.Printf("\tdone (cloud: %s, party %s)\n", elapsedGKGCloud, elapsedGKGParty)

	return
}

func genquery(params bfv.Parameters, queryIndex int, encoder *bfv.Encoder, encryptor rlwe.EncryptorInterface) *rlwe.Ciphertext {
	// Query ciphertext
	queryCoeffs := make([]uint64, params.N())
	queryCoeffs[queryIndex] = 1
	query := bfv.NewPlaintext(params, params.MaxLevel())
	var encQuery *rlwe.Ciphertext
	elapsedRequestParty += runTimed(func() {
		if err := encoder.Encode(queryCoeffs, query); err != nil {
			panic(err)
		}
		encQuery = encryptor.EncryptNew(query)
	})

	return encQuery
}

func requestphase(params bfv.Parameters, queryIndex, NGoRoutine int, encQuery *rlwe.Ciphertext, encInputs []*rlwe.Ciphertext, plainMask []*rlwe.Plaintext, evk rlwe.EvaluationKeySet) *rlwe.Ciphertext {

	l := log.New(os.Stderr, "", 0)

	l.Println("> Request Phase")

	// Buffer for the intermediate computation done by the cloud
	encPartial := make([]*rlwe.Ciphertext, len(encInputs))
	for i := range encPartial {
		encPartial[i] = bfv.NewCiphertext(params, 2, params.MaxLevel())
	}

	evaluator := bfv.NewEvaluator(params, evk)

	// Split the task among the Go routines
	tasks := make(chan *maskTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {
			evaluator := evaluator.ShallowCopy() // creates a shallow evaluator copy for this goroutine
			tmp := bfv.NewCiphertext(params, 1, params.MaxLevel())
			for task := range tasks {
				task.elapsedmaskTask = runTimed(func() {
					// 1) Multiplication of the query with the plaintext mask
					evaluator.Mul(task.query, task.mask, tmp)

					// 2) Inner sum (populate all the slots with the sum of all the slots)
					evaluator.InnerSum(tmp, 1, params.N()>>1, tmp)
					evaluator.Add(tmp, evaluator.RotateRowsNew(tmp), tmp)

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

	resultDeg2 := bfv.NewCiphertext(params, 2, params.MaxLevel())
	result := bfv.NewCiphertext(params, 1, params.MaxLevel())

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
