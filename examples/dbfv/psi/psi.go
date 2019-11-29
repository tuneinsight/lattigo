package main

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/dbfv"
	"github.com/ldsec/lattigo/ring"
	"log"
	"math/rand"
	"os"
	"strconv"
	"sync"
	"time"
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

	l := log.New(os.Stderr, "", 0)

	// Largest for n=8192 : 512 parties
	N := 8
	var err error
	if len(os.Args[1:]) >= 1 {
		N, err = strconv.Atoi(os.Args[1])
		check(err)
	}

	NGoRoutine := 1
	if len(os.Args[1:]) >= 2 {
		NGoRoutine, err = strconv.Atoi(os.Args[2])
		check(err)
	}

	type party struct {
		sk         *bfv.SecretKey
		rlkEphemSk *ring.Poly

		ckgShare      dbfv.CKGShare
		rkgShareOne   dbfv.RKGShareRoundOne
		rkgShareTwo   dbfv.RKGShareRoundTwo
		rkgShareThree dbfv.RKGShareRoundThree
		pcksShare     dbfv.PCKSShare

		input []uint64
	}

	params := bfv.DefaultParams[14]
	params.T = 65537
	bfvctx := bfv.NewContext(params)

	crsGen := ring.NewCRPGenerator([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'}, bfvctx.ContextKeys())
	crs := crsGen.ClockNew()
	crp := make([]*ring.Poly, bfvctx.Beta())
	for i := uint64(0); i < bfvctx.Beta(); i++ {
		crp[i] = crsGen.ClockNew()
	}

	tsk, tpk := bfv.NewKeyGenerator(params).NewKeyPair()
	colSk := &bfv.SecretKey{}
	colSk.Set(bfvctx.ContextKeys().NewPoly())

	expRes := make([]uint64, 1<<params.LogN, 1<<params.LogN)
	for i := range expRes {
		expRes[i] = 1
	}

	ckg := dbfv.NewCKGProtocol(params)
	rkg := dbfv.NewEkgProtocol(params)
	pcks := dbfv.NewPCKSProtocol(params, 3.19)

	P := make([]*party, N, N)
	for i := range P {
		pi := &party{}
		pi.sk = bfv.NewKeyGenerator(params).NewSecretKey()
		pi.rlkEphemSk = bfvctx.ContextKeys().SampleTernaryMontgomeryNTTNew(1.0 / 3)
		pi.input = make([]uint64, 1<<params.LogN, 1<<params.LogN)
		for i := range pi.input {
			if rand.Float32() > 0.3 || i == 4 {
				pi.input[i] = 1
			}
			expRes[i] *= pi.input[i]
		}
		bfvctx.ContextKeys().Add(colSk.Get(), pi.sk.Get(), colSk.Get()) //TODO: doc says "return"

		pi.ckgShare = ckg.AllocateShares()
		pi.rkgShareOne, pi.rkgShareTwo, pi.rkgShareThree = rkg.AllocateShares()
		pi.pcksShare = pcks.AllocateShares()

		P[i] = pi
	}

	var elapsedCKGCloud time.Duration
	var elapsedCKGParty time.Duration
	var elapsedRKGCloud time.Duration
	var elapsedRKGParty time.Duration

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

	l.Println("> RKG Phase")
	elapsedRKGParty = runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundOne(pi.rlkEphemSk, pi.sk.Get(), crp, pi.rkgShareOne)
		}
	}, N)

	rkgCombined1, rkgCombined2, rkgCombined3 := rkg.AllocateShares()

	elapsedRKGCloud = runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShareRoundOne(pi.rkgShareOne, rkgCombined1, rkgCombined1)
		}
	})

	elapsedRKGParty += runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundTwo(rkgCombined1, pi.sk.Get(), crp, pi.rkgShareTwo)
		}
	}, N)

	elapsedRKGCloud += runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShareRoundTwo(pi.rkgShareTwo, rkgCombined2, rkgCombined2)
		}
	})

	elapsedRKGParty += runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundThree(rkgCombined2, pi.rlkEphemSk, pi.sk.Get(), pi.rkgShareThree)
		}
	}, N)

	rlk := bfv.NewRelinKey(params, 1)
	elapsedRKGCloud += runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShareRoundThree(pi.rkgShareThree, rkgCombined3, rkgCombined3)
		}
		rkg.GenRelinearizationKey(rkgCombined2, rkgCombined3, rlk)
	})
	l.Printf("\tdone (cloud: %s, party: %s)\n",
		elapsedRKGCloud, elapsedRKGParty)
	l.Printf("\tSetup done (cloud: %s, party: %s)\n",
		elapsedRKGCloud+elapsedCKGCloud, elapsedRKGParty+elapsedCKGParty)

	// Pre-loading memory
	l.Println("> Memory alloc Phase")
	encInputs := make([]*bfv.Ciphertext, N, N)
	for i := range encInputs {
		encInputs[i] = bfv.NewCiphertext(params, 1)
	}

	encLvls := make([][]*bfv.Ciphertext, 0)
	encLvls = append(encLvls, encInputs)
	for nLvl := N / 2; nLvl > 0; nLvl = nLvl >> 1 {
		encLvl := make([]*bfv.Ciphertext, nLvl, nLvl)
		for i := range encLvl {
			encLvl[i] = bfv.NewCiphertext(params, 2)
		}
		encLvls = append(encLvls, encLvl)
	}
	encRes := encLvls[len(encLvls)-1][0]

	l.Println("> Encrypt Phase")
	encryptor := bfv.NewEncryptorFromPk(params, pk)
	encoder := bfv.NewEncoder(params)
	pt := bfv.NewPlaintextFromParams(params)
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

	tasks := make(chan *MultTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)
	//l.Println("> Spawning", NGoRoutine, "evaluator goroutine")
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {
			evaluator := bfv.NewEvaluator(params)
			for task := range tasks {
				task.elapsedMultTask = runTimed(func() {
					evaluator.Mul(task.op1, task.op2, task.res)
					evaluator.Relinearize(task.res, rlk, task.res)
				})
				task.wg.Done()
			}
			//l.Println("\t evaluator", i, "down")
			workers.Done()
		}(i)
		//l.Println("\t evaluator", i, "started")
	}

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

	l.Println("> Result:")
	decryptor := bfv.NewDecryptor(params, tsk)
	ptres := bfv.NewPlaintextFromParams(params)
	elapsedDecParty := runTimed(func() {
		decryptor.Decrypt(encOut, ptres) // TODO : manage error
	})

	//check(err)
	res := encoder.DecodeUint(ptres)
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
