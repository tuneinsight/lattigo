package main

import (
	"log"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/dbfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
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
	// For more details about the PSI example see
	//     Multiparty Homomorphic Encryption: From Theory to Practice (<https://eprint.iacr.org/2020/304>)

	l := log.New(os.Stderr, "", 0)

	// $go run main.go arg1 arg2
	// arg1: number of parties
	// arg2: number of Go routines

	// Largest for n=8192: 512 parties
	N := 8 // Default number of parties
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

	type party struct {
		sk         *bfv.SecretKey
		rlkEphemSk *ring.Poly

		ckgShare    dbfv.CKGShare
		rkgShareOne dbfv.RKGShare
		rkgShareTwo dbfv.RKGShare
		pcksShare   dbfv.PCKSShare

		input []uint64
	}

	// Use the defaultparams logN=14, logQP=438 with a plaintext modulus T=65537
	params := bfv.DefaultParams[bfv.PN14QP438].WithT(65537)

	// PRNG keyed with "lattigo"
	lattigoPRNG, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}

	// Ring for the common reference polynomials sampling
	ringQP, _ := ring.NewRing(1<<params.LogN(), append(params.Qi(), params.Pi()...))

	// Common reference polynomial generator that uses the PRNG
	crsGen := ring.NewUniformSampler(lattigoPRNG, ringQP)
	crs := crsGen.ReadNew()                  // for the public-key
	crp := make([]*ring.Poly, params.Beta()) // for the relinearization keys
	for i := uint64(0); i < params.Beta(); i++ {
		crp[i] = crsGen.ReadNew()
	}

	// Target private and public keys
	tsk, tpk := bfv.NewKeyGenerator(params).GenKeyPair()

	expRes := make([]uint64, 1<<params.LogN())
	for i := range expRes {
		expRes[i] = 1
	}

	ckg := dbfv.NewCKGProtocol(params)
	rkg := dbfv.NewEkgProtocol(params)
	pcks := dbfv.NewPCKSProtocol(params, 3.19)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, ringQP, 0.5, true)

	// Create each party, and allocate the memory for all the shares that the protocols will need
	P := make([]*party, N)
	for i := range P {
		pi := &party{}
		pi.sk = bfv.NewKeyGenerator(params).GenSecretKey()

		pi.rlkEphemSk = ternarySamplerMontgomery.ReadNew()
		ringQP.NTT(pi.rlkEphemSk, pi.rlkEphemSk)

		pi.input = make([]uint64, 1<<params.LogN())
		for i := range pi.input {
			if utils.RandFloat64(0, 1) > 0.3 || i == 4 {
				pi.input[i] = 1
			}
			expRes[i] *= pi.input[i]
		}

		pi.ckgShare = ckg.AllocateShares()
		pi.rkgShareOne, pi.rkgShareTwo = rkg.AllocateShares()
		pi.pcksShare = pcks.AllocateShares()

		P[i] = pi
	}

	var elapsedCKGCloud time.Duration
	var elapsedCKGParty time.Duration
	var elapsedRKGCloud time.Duration
	var elapsedRKGParty time.Duration

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

	l.Printf("\tdone (cloud: %s, party: %s)\n",
		elapsedRKGCloud, elapsedRKGParty)
	l.Printf("\tSetup done (cloud: %s, party: %s)\n",
		elapsedRKGCloud+elapsedCKGCloud, elapsedRKGParty+elapsedCKGParty)

	// Pre-loading memory
	l.Println("> Memory alloc Phase")
	encInputs := make([]*bfv.Ciphertext, N)
	for i := range encInputs {
		encInputs[i] = bfv.NewCiphertext(params, 1)
	}

	encLvls := make([][]*bfv.Ciphertext, 0)
	encLvls = append(encLvls, encInputs)
	for nLvl := N / 2; nLvl > 0; nLvl = nLvl >> 1 {
		encLvl := make([]*bfv.Ciphertext, nLvl)
		for i := range encLvl {
			encLvl[i] = bfv.NewCiphertext(params, 2)
		}
		encLvls = append(encLvls, encLvl)
	}
	encRes := encLvls[len(encLvls)-1][0]

	// Each party encrypts its input vector
	l.Println("> Encrypt Phase")
	encryptor := bfv.NewEncryptorFromPk(params, pk)
	encoder := bfv.NewEncoder(params)
	pt := bfv.NewPlaintext(params)
	elapsedEncryptParty := runTimedParty(func() {
		for i, pi := range P {
			encoder.EncodeUint(pi.input, pt)
			encryptor.Encrypt(pt, encInputs[i])
		}
	}, N)

	elapsedEncryptCloud := time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	type MultTask struct {
		wg              *sync.WaitGroup
		op1             *bfv.Ciphertext
		op2             *bfv.Ciphertext
		res             *bfv.Ciphertext
		elapsedMultTask time.Duration
	}

	// Split the task among the Go routines
	tasks := make(chan *MultTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)
	//l.Println("> Spawning", NGoRoutine, "evaluator goroutine")
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {
			evaluator := bfv.NewEvaluator(params)
			for task := range tasks {
				task.elapsedMultTask = runTimed(func() {
					// 1) Multiplication of two input vectors
					evaluator.Mul(task.op1, task.op2, task.res)
					// 2) Relinearization
					evaluator.Relinearize(task.res, rlk, task.res)
				})
				task.wg.Done()
			}
			//l.Println("\t evaluator", i, "down")
			workers.Done()
		}(i)
		//l.Println("\t evaluator", i, "started")
	}

	// Start the tasks
	taskList := make([]*MultTask, 0)
	l.Println("> Eval Phase")
	elapsedEvalCloud := runTimed(func() {
		for i, lvl := range encLvls[:len(encLvls)-1] {
			nextLvl := encLvls[i+1]
			l.Println("\tlevel", i, len(lvl), "->", len(nextLvl))
			wg := &sync.WaitGroup{}
			wg.Add(len(nextLvl))
			for j, nextLvlCt := range nextLvl {
				task := MultTask{wg, lvl[2*j], lvl[2*j+1], nextLvlCt, 0}
				taskList = append(taskList, &task)
				tasks <- &task
			}
			wg.Wait()
		}
	})
	elapsedEvalCloudCPU := time.Duration(0)
	for _, t := range taskList {
		elapsedEvalCloudCPU += t.elapsedMultTask
	}
	elapsedEvalParty := time.Duration(0)
	l.Printf("\tdone (cloud: %s (wall: %s), party: %s)\n",
		elapsedEvalCloudCPU, elapsedEvalCloud, elapsedEvalParty)

	//l.Println("> Shutting down workers")
	close(tasks)
	workers.Wait()

	// Collective key switching from the collective secret key to
	// the target public key
	l.Println("> PCKS Phase")
	elapsedPCKSParty := runTimedParty(func() {
		for _, pi := range P {
			pcks.GenShare(pi.sk.Get(), tpk, encRes, pi.pcksShare)
		}
	}, N)

	pcksCombined := pcks.AllocateShares()
	encOut := bfv.NewCiphertext(params, 1)
	elapsedPCKSCloud := runTimed(func() {
		for _, pi := range P {
			pcks.AggregateShares(pi.pcksShare, pcksCombined, pcksCombined)
		}
		pcks.KeySwitch(pcksCombined, encRes, encOut)

	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedPCKSCloud, elapsedPCKSParty)

	// Decrypt the result with the target secret key
	l.Println("> Result:")
	decryptor := bfv.NewDecryptor(params, tsk)
	ptres := bfv.NewPlaintext(params)
	elapsedDecParty := runTimed(func() {
		decryptor.Decrypt(encOut, ptres)
	})

	// Check the result
	res := encoder.DecodeUintNew(ptres)
	l.Printf("\t%v\n", res[:16])
	for i := range expRes {
		if expRes[i] != res[i] {
			//l.Printf("\t%v\n", expRes)
			l.Println("\tincorrect")
			return
		}
	}
	l.Println("\tcorrect")
	l.Printf("> Finished (total cloud: %s, total party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedEncryptCloud+elapsedEvalCloud+elapsedPCKSCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedEncryptParty+elapsedEvalParty+elapsedPCKSParty+elapsedDecParty)

}
