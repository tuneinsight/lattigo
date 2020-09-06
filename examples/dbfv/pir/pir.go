package main

import (
	"log"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/dbfv"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
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

func main() {

	// This examples simulates a SMC instance of a private information retrieval.
	// The situation is as follow : a cloud stores data of several parties
	// encrypted under a shared public-key. An external party wants to retrieve
	// the plaintext content of one of the ciphertexts while ensuring the following
	// security property : no information other than that a request was made must
	// leak to the cloud, to the owners of the shared public-key or to anyone else.
	//
	// For more details cf : paper

	l := log.New(os.Stderr, "", 0)

	// $go run pir.go arg1 arg2
	// arg1 : number of parties
	// arg2 : number of Go routines
	// MinDelta number of parties for n=8192 : 512 parties (!memory intensive!)

	N := 3 // default number of parties
	var err error
	if len(os.Args[1:]) >= 1 {
		N, err = strconv.Atoi(os.Args[1])
		check(err)
	}

	NGoRoutine := 1 // default number of Go routines
	if len(os.Args[1:]) >= 2 {
		NGoRoutine, err = strconv.Atoi(os.Args[2])
		check(err)
	}

	// Index of the ciphertext to retrieve.
	queryIndex := 2

	type party struct {
		sk         *bfv.SecretKey
		rlkEphemSk *ring.Poly

		ckgShare    dbfv.CKGShare
		rkgShareOne dbfv.RKGShare
		rkgShareTwo dbfv.RKGShare
		rtgShare    dbfv.RTGShare
		cksShare    dbfv.CKSShare

		input []uint64
	}

	params := bfv.DefaultParams[bfv.PN13QP218].WithT(65537) // default params with N=8192

	// Common reference polynomial generator keyed with
	// "lattigo"
	crsGen := dbfv.NewCRPGenerator(params, []byte{'l', 'a', 't', 't', 'i', 'g', 'o'})

	// Generation of the common reference polynomials
	crs := crsGen.ReadNew()                     // for the public-key
	crp := make([]*ring.Poly, params.Beta())    // for the relinearization keys
	crpRot := make([]*ring.Poly, params.Beta()) // for the rotation keys
	for i := uint64(0); i < params.Beta(); i++ {
		crp[i] = crsGen.ReadNew()
	}
	for i := uint64(0); i < params.Beta(); i++ {
		crpRot[i] = crsGen.ReadNew()
	}

	// Collective secret key = sum(individual secret keys)
	colSk := bfv.NewSecretKey(params)

	// Instantiation of each of the protocols needed for the pir example
	ckg := dbfv.NewCKGProtocol(params)       // public key generation
	rkg := dbfv.NewEkgProtocol(params)       // relineariation key generation
	rtg := dbfv.NewRotKGProtocol(params)     // rotation keys generation
	cks := dbfv.NewCKSProtocol(params, 3.19) // collective public-key re-encryption

	kgen := bfv.NewKeyGenerator(params)

	contextKeys, _ := ring.NewContextWithParams(1<<params.LogN(), append(params.Qi(), params.Pi()...))
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, contextKeys, 0.5, true)

	// Creates each party, and allocates the memory for all the shares that the protocols will need
	P := make([]*party, N, N)
	for i := range P {
		pi := &party{}
		pi.sk = kgen.GenSecretKey()
		pi.rlkEphemSk = ternarySamplerMontgomery.ReadNew()
		contextKeys.NTT(pi.rlkEphemSk, pi.rlkEphemSk)
		pi.input = make([]uint64, 1<<params.LogN(), 1<<params.LogN())
		for j := range pi.input {
			pi.input[j] = uint64(i)
		}
		contextKeys.Add(colSk.Get(), pi.sk.Get(), colSk.Get()) //TODO: doc says "return"

		pi.ckgShare = ckg.AllocateShares()
		pi.rkgShareOne, pi.rkgShareTwo = rkg.AllocateShares()
		pi.rtgShare = rtg.AllocateShare()
		pi.cksShare = cks.AllocateShare()

		P[i] = pi
	}

	var elapsedCKGCloud time.Duration
	var elapsedCKGParty time.Duration
	var elapsedRKGCloud time.Duration
	var elapsedRKGParty time.Duration
	var elapsedRTGCloud time.Duration
	var elapsedRTGParty time.Duration

	// 1) Collective public key generation
	l.Println("> CKG Phase")
	pk := bfv.NewPublicKey(params)
	elapsedCKGParty = runTimedParty(func() {
		for _, pi := range P {
			ckg.GenShare(pi.sk.Get(), crs, pi.ckgShare)
		}
	}, N)
	ckgCombined := ckg.AllocateShares()
	elapsedCKGCloud = runTimed(func() {
		for _, pi := range P {
			ckg.AggregateShares(pi.ckgShare, ckgCombined, ckgCombined)
		}
		ckg.GenPublicKey(ckgCombined, crs, pk)
	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKGCloud, elapsedCKGParty)

	// 2) Collective relinearization key generation
	l.Println("> RKG Phase")
	elapsedRKGParty = runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundOne(pi.rlkEphemSk, pi.sk.Get(), crp, pi.rkgShareOne)
		}
	}, N)

	rkgCombined1, rkgCombined2 := rkg.AllocateShares()

	elapsedRKGCloud = runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShareRoundOne(pi.rkgShareOne, rkgCombined1, rkgCombined1)
		}
	})

	elapsedRKGParty += runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundTwo(rkgCombined1, pi.rlkEphemSk, pi.sk.Get(), crp, pi.rkgShareTwo)
		}
	}, N)

	rlk := bfv.NewRelinKey(params, 1)
	elapsedRKGCloud += runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShareRoundTwo(pi.rkgShareTwo, rkgCombined2, rkgCombined2)
		}
		rkg.GenRelinearizationKey(rkgCombined1, rkgCombined2, rlk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedRKGCloud, elapsedRKGParty)

	// 3) Collective rotation keys geneneration
	l.Println("> RTG Phase")
	rtk := bfv.NewRotationKeys()
	for _, rot := range []bfv.Rotation{bfv.RotationRight, bfv.RotationLeft, bfv.RotationRow} {
		for k := uint64(1); (rot == bfv.RotationRow && k == 1) || (rot != bfv.RotationRow && k < 1<<(params.LogN()-1)); k <<= 1 {

			rtgShareCombined := rtg.AllocateShare()
			rtgShareCombined.Type = rot
			rtgShareCombined.K = k

			elapsedRTGParty += runTimedParty(func() {
				for _, pi := range P {
					rtg.GenShare(rot, k, pi.sk.Get(), crpRot, &pi.rtgShare)
				}
			}, N)

			elapsedRTGCloud += runTimed(func() {
				for _, pi := range P {
					rtg.Aggregate(pi.rtgShare, rtgShareCombined, rtgShareCombined)
				}
				rtg.Finalize(rtgShareCombined, crpRot, rtk)
			})
		}
	}
	l.Printf("\tdone (cloud: %s, party %s)\n", elapsedRTGCloud, elapsedRTGParty)

	l.Printf("\tSetup done (cloud: %s, party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedRTGCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedRTGParty)

	// Pre-loading memory
	encoder := bfv.NewEncoder(params)
	l.Println("> Memory alloc Phase")
	encInputs := make([]*bfv.Ciphertext, N, N)
	plainMask := make([]*bfv.Plaintext, N, N)
	encPartial := make([]*bfv.Ciphertext, N, N)

	// Ciphertexts to be retrieved.
	for i := range encInputs {
		encInputs[i] = bfv.NewCiphertext(params, 1)
	}

	// Plaintext masks : plainmask[i] = encode([0, ..., 0, 1_i, 0, ..., 0])
	// (zero with a 1 at the ith position).
	for i := range plainMask {
		maskCoeffs := make([]uint64, 1<<params.LogN())
		maskCoeffs[i] = 1
		plainMask[i] = bfv.NewPlaintext(params)
		encoder.EncodeUint(maskCoeffs, plainMask[i])
	}

	// Buffer for the intermediate compuation done by the cloud.
	for i := range encPartial {
		encPartial[i] = bfv.NewCiphertext(params, 2)
	}

	// Ciphertexts encrypted under CPK and stored in the cloud.
	l.Println("> Encrypt Phase")
	encryptor := bfv.NewEncryptorFromPk(params, pk)
	pt := bfv.NewPlaintext(params)
	elapsedEncryptParty := runTimedParty(func() {
		for i, pi := range P {
			encoder.EncodeUint(pi.input, pt)
			encryptor.Encrypt(pt, encInputs[i])
		}
	}, N)

	elapsedEncryptCloud := time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	l.Println("> Request Phase")
	var elapsedRequestParty time.Duration
	var elapsedRequestCloud time.Duration
	var elapsedRequestCloudCPU time.Duration

	// Query ciphertext
	queryCoeffs := make([]uint64, 1<<params.LogN())
	queryCoeffs[queryIndex] = 1
	query := bfv.NewPlaintext(params)
	var encQuery *bfv.Ciphertext
	elapsedRequestParty += runTimed(func() {
		encoder.EncodeUint(queryCoeffs, query)
		encQuery = encryptor.EncryptNew(query)
	})

	type MaskTask struct {
		query           *bfv.Ciphertext
		mask            *bfv.Plaintext
		row             *bfv.Ciphertext
		res             *bfv.Ciphertext
		elapsedMaskTask time.Duration
	}

	// Splits the task among the Go routines
	tasks := make(chan *MaskTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {
			evaluator := bfv.NewEvaluator(params)
			tmp := bfv.NewCiphertext(params, 1)
			for task := range tasks {
				task.elapsedMaskTask = runTimed(func() {
					// 1) Multiplication of the query with the plaintext mask.
					evaluator.Mul(task.query, task.mask, tmp)
					// 2) Inner sum (populates all the slots with the sum of all the slots).
					evaluator.InnerSum(tmp, rtk, tmp)
					// 3) Multiplication of 2) with the ith ciphertext stored in the cloud.
					evaluator.Mul(tmp, task.row, task.res)
				})
			}
			//l.Println("\t evaluator", i, "down")
			workers.Done()
		}(i)
		//l.Println("\t evaluator", i, "started")
	}

	taskList := make([]*MaskTask, 0)

	elapsedRequestCloud += runTimed(func() {
		for i := range encInputs {
			task := &MaskTask{
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
		elapsedRequestCloudCPU += t.elapsedMaskTask
	}

	evaluator := bfv.NewEvaluator(params)
	resultDeg2 := bfv.NewCiphertext(params, 2)
	result := bfv.NewCiphertext(params, 1)

	// Summation of all the partial result among the different Go routines
	finalAddDuration := runTimed(func() {
		for i := 0; i < N; i++ {
			evaluator.Add(resultDeg2, encPartial[i], resultDeg2)
		}
		evaluator.Relinearize(resultDeg2, rlk, result)
	})

	elapsedRequestCloud += finalAddDuration
	elapsedRequestCloudCPU += finalAddDuration

	l.Printf("\tdone (cloud: %s/%s, party: %s)\n",
		elapsedRequestCloud, elapsedRequestCloudCPU, elapsedRequestParty)

	// Collective (partial) decryption.
	l.Println("> CKS Phase")
	zero := params.NewPolyQ()
	cksCombined := cks.AllocateShare()
	elapsedPCKSParty := runTimedParty(func() {
		for _, pi := range P[1:] {
			cks.GenShare(pi.sk.Get(), zero, result, pi.cksShare)
		}
	}, N-1)

	encOut := bfv.NewCiphertext(params, 1)
	elapsedCKSCloud := runTimed(func() {
		for _, pi := range P {
			cks.AggregateShares(pi.cksShare, cksCombined, cksCombined)
		}
		cks.KeySwitch(cksCombined, result, encOut)
	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKSCloud, elapsedPCKSParty)

	l.Println("> Result:")

	// Decryption by the external party
	decryptor := bfv.NewDecryptor(params, P[0].sk)
	ptres := bfv.NewPlaintext(params)
	elapsedDecParty := runTimed(func() {
		decryptor.Decrypt(encOut, ptres)
	})

	res := encoder.DecodeUint(ptres)

	l.Printf("\t%v\n", res[:16])
	l.Printf("> Finished (total cloud: %s, total party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedRTGCloud+elapsedEncryptCloud+elapsedRequestCloudCPU+elapsedCKSCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedRTGParty+elapsedEncryptParty+elapsedRequestParty+elapsedPCKSParty+elapsedDecParty)

}
