package main

import (
	"log"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v3/bfv"
	"github.com/tuneinsight/lattigo/v3/dbfv"
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
	sk            *rlwe.SecretKey
	thresholdizer *drlwe.Thresholdizer
	combiner      drlwe.Combiner
	gen           *drlwe.ShamirPolynomial
	rlkEphemSk    *rlwe.SecretKey
	tsks          *drlwe.ShamirSecretShare
	tsk           *rlwe.SecretKey
	tpk           drlwe.ShamirPublicKey
	ckgShare      *drlwe.CKGShare
	rkgShareOne   *drlwe.RKGShare
	rkgShareTwo   *drlwe.RKGShare
	pcksShare     *drlwe.PCKSShare

	input []uint64
}
type multTask struct {
	wg              *sync.WaitGroup
	op1             *bfv.Ciphertext
	op2             *bfv.Ciphertext
	res             *bfv.Ciphertext
	elapsedmultTask time.Duration
}

var elapsedEncryptParty time.Duration
var elapsedEncryptCloud time.Duration
var elapsedCKGCloud time.Duration
var elapsedCKGParty time.Duration
var elapsedRKGCloud time.Duration
var elapsedRKGParty time.Duration
var elapsedPCKSCloud time.Duration
var elapsedPCKSParty time.Duration
var elapsedEvalCloudCPU time.Duration
var elapsedEvalCloud time.Duration
var elapsedEvalParty time.Duration
var elapsedThreshCombCloud time.Duration
var elapsedThreshCombParty time.Duration
var elapsedThreshShareCloud time.Duration
var elapsedThreshShareParty time.Duration
var elapsedThreshInitParty time.Duration

var l *log.Logger

func main() {
	// For more details about the PSI example see
	//     Multiparty Homomorphic Encryption: From Theory to Practice (<https://eprint.iacr.org/2020/304>)

	l = log.New(os.Stderr, "", 0)

	// $go run main.go arg1 arg2 arg3
	// arg1: number of parties
	// arg2: number of Go routines
	// arg3: threshold value

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

	t := N / 2 // Default threshold value
	if len(os.Args[1:]) >= 3 {
		t, err = strconv.Atoi(os.Args[3])
		check(err)
	}

	paramsDef := bfv.PN14QP438
	paramsDef.T = 65537
	params, err := bfv.NewParametersFromLiteral(paramsDef)
	if err != nil {
		panic(err)
	}

	// PRNG keyed with "lattigo"
	crs, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}

	encoder := bfv.NewEncoder(params)

	// Target private and public keys
	tsk, tpk := bfv.NewKeyGenerator(params).GenKeyPair()

	// Create each party, and allocate the memory for all the shares that the protocols will need
	P := genparties(params, N)

	genThresholdizers(params, P, t)

	thresholdGenShares(params, P)
	// Inputs & expected result
	expRes := genInputs(params, P)

	// 1) Collective public key generation
	pk := ckgphase(params, crs, P)
	// 2) Collective relinearization key generation
	rlk := rkgphase(params, crs, P)

	l.Printf("\tdone (cloud: %s, party: %s)\n",
		elapsedRKGCloud, elapsedRKGParty)
	l.Printf("Setup done (cloud: %s, party: %s)\n",
		elapsedRKGCloud+elapsedCKGCloud+elapsedThreshShareCloud, elapsedRKGParty+elapsedCKGParty+elapsedThreshShareParty+elapsedThreshInitParty)

	encInputs := encPhase(params, P, pk, encoder)

	encRes := evalPhase(params, NGoRoutine, encInputs, rlk)

	activeParties := thresholdCombine(params, P, uint64(t))
	//Active players have now a different secret key share

	//Only active players take part in the key switching protocol
	encOut := pcksPhase(params, tpk, encRes, activeParties)

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
		elapsedCKGCloud+elapsedRKGCloud+elapsedEncryptCloud+elapsedEvalCloud+elapsedPCKSCloud+elapsedThreshShareCloud+elapsedThreshCombCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedEncryptParty+elapsedEvalParty+elapsedPCKSParty+
			elapsedDecParty+elapsedThreshShareParty+elapsedThreshCombParty+elapsedThreshInitParty)

}

func thresholdGenShares(params bfv.Parameters, P []*party) {

	l.Println("> Threshold Shares Generation")

	// Each party generates the polynomial shares for its own Secret Key
	shares := make(map[*party]map[*party]*drlwe.ShamirSecretShare, len(P))
	elapsedThreshShareParty = runTimedParty(func() {
		for _, pi := range P {

			shares[pi] = make(map[*party]*drlwe.ShamirSecretShare)

			for _, pj := range P {
				shares[pi][pj] = pi.thresholdizer.AllocateThresholdSecretShare()
				pi.thresholdizer.GenShamirSecretShare(pj.tpk, pi.gen, shares[pi][pj])
			}
		}
	}, len(P))

	elapsedThreshShareCloud = runTimed(func() {
		for _, pi := range P {
			for _, pj := range P {
				pi.thresholdizer.AggregateShares(pi.tsks, shares[pj][pi], pi.tsks)
			}
		}
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedThreshShareCloud, elapsedThreshShareParty)
}

func thresholdCombine(params bfv.Parameters, P []*party, t uint64) (activeParties []*party) {

	l.Println("> ThresholdCombine")
	activeParties = make([]*party, int(t))
	activePoints := make([]drlwe.ShamirPublicKey, int(t))

	//Determining which players are active and their key
	elapsedThreshCombCloud = runTimed(func() {
		for i := range P[:t] {
			activeParties[i] = P[i]
			activePoints[i] = P[i].tpk
		}
	})

	//Combining
	elapsedThreshCombParty = runTimedParty(func() {
		for _, pi := range activeParties {
			pi.combiner.GenAdditiveShare(activePoints, pi.tpk, pi.tsks, pi.tsk)
		}
	}, len(P))

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedThreshCombCloud, elapsedThreshCombParty)

	return
}

func encPhase(params bfv.Parameters, P []*party, pk *rlwe.PublicKey, encoder bfv.Encoder) (encInputs []*bfv.Ciphertext) {

	encInputs = make([]*bfv.Ciphertext, len(P))
	for i := range encInputs {
		encInputs[i] = bfv.NewCiphertext(params, 1)
	}

	// Each party encrypts its input vector
	l.Println("> Encrypt Phase")
	encryptor := bfv.NewEncryptor(params, pk)

	pt := bfv.NewPlaintext(params)
	elapsedEncryptParty = runTimedParty(func() {
		for i, pi := range P {
			encoder.Encode(pi.input, pt)
			encryptor.Encrypt(pt, encInputs[i])
		}
	}, len(P))

	elapsedEncryptCloud = time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	return
}

func evalPhase(params bfv.Parameters, NGoRoutine int, encInputs []*bfv.Ciphertext, rlk *rlwe.RelinearizationKey) (encRes *bfv.Ciphertext) {

	encLvls := make([][]*bfv.Ciphertext, 0)
	encLvls = append(encLvls, encInputs)
	for nLvl := len(encInputs) / 2; nLvl > 0; nLvl = nLvl >> 1 {
		encLvl := make([]*bfv.Ciphertext, nLvl)
		for i := range encLvl {
			encLvl[i] = bfv.NewCiphertext(params, 2)
		}
		encLvls = append(encLvls, encLvl)
	}
	encRes = encLvls[len(encLvls)-1][0]

	evaluator := bfv.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: nil})
	// Split the task among the Go routines
	tasks := make(chan *multTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)
	//l.Println("> Spawning", NGoRoutine, "evaluator goroutine")
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {
			evaluator := evaluator.ShallowCopy() // creates a shallow evaluator copy for this goroutine
			for task := range tasks {
				task.elapsedmultTask = runTimed(func() {
					// 1) Multiplication of two input vectors
					evaluator.Mul(task.op1, task.op2, task.res)
					// 2) Relinearization
					evaluator.Relinearize(task.res, task.res)
				})
				task.wg.Done()
			}
			//l.Println("\t evaluator", i, "down")
			workers.Done()
		}(i)
		//l.Println("\t evaluator", i, "started")
	}

	// Start the tasks
	taskList := make([]*multTask, 0)
	l.Println("> Eval Phase")
	elapsedEvalCloud = runTimed(func() {
		for i, lvl := range encLvls[:len(encLvls)-1] {
			nextLvl := encLvls[i+1]
			l.Println("\tlevel", i, len(lvl), "->", len(nextLvl))
			wg := &sync.WaitGroup{}
			wg.Add(len(nextLvl))
			for j, nextLvlCt := range nextLvl {
				task := multTask{wg, lvl[2*j], lvl[2*j+1], nextLvlCt, 0}
				taskList = append(taskList, &task)
				tasks <- &task
			}
			wg.Wait()
		}
	})
	elapsedEvalCloudCPU = time.Duration(0)
	for _, t := range taskList {
		elapsedEvalCloudCPU += t.elapsedmultTask
	}
	elapsedEvalParty = time.Duration(0)
	l.Printf("\tdone (cloud: %s (wall: %s), party: %s)\n",
		elapsedEvalCloudCPU, elapsedEvalCloud, elapsedEvalParty)

	//l.Println("> Shutting down workers")
	close(tasks)
	workers.Wait()

	return
}

func genparties(params bfv.Parameters, N int) []*party {

	// Create each party, and allocate the memory for all the shares that the protocols will need
	P := make([]*party, N)
	for i := range P {
		pi := &party{}
		pi.sk = bfv.NewKeyGenerator(params).GenSecretKey()
		pi.tsk = bfv.NewSecretKey(params)
		P[i] = pi
	}

	return P
}

func genThresholdizers(params bfv.Parameters, P []*party, t int) {
	l.Println("> Thresholdizers Initialization")
	elapsedThreshInitParty = runTimedParty(func() {
		var err error
		for i, pi := range P {
			pi.thresholdizer = drlwe.NewThresholdizer(params.Parameters)
			pi.combiner = drlwe.NewCombiner(params.Parameters, t)
			pi.tpk = drlwe.ShamirPublicKey(i + 1)
			pi.gen, err = pi.thresholdizer.GenShamirPolynomial(t, pi.sk)
			if err != nil {
				panic(err)
			}
			pi.tsks = pi.thresholdizer.AllocateThresholdSecretShare()
		}
	}, len(P))
	l.Printf("\tdone (party : %s)\n", elapsedThreshInitParty)
}

func genInputs(params bfv.Parameters, P []*party) (expRes []uint64) {

	expRes = make([]uint64, params.N())
	for i := range expRes {
		expRes[i] = 1
	}

	for _, pi := range P {

		pi.input = make([]uint64, params.N())
		for i := range pi.input {
			if utils.RandFloat64(0, 1) > 0.3 || i == 4 {
				pi.input[i] = 1
			}
			expRes[i] *= pi.input[i]
		}

	}

	return
}

func pcksPhase(params bfv.Parameters, tpk *rlwe.PublicKey, encRes *bfv.Ciphertext, P []*party) (encOut *bfv.Ciphertext) {

	// Collective key switching from the collective secret key to
	// the target public key

	pcks := dbfv.NewPCKSProtocol(params, 3.19)

	for _, pi := range P {
		pi.pcksShare = pcks.AllocateShare()
	}

	l.Println("> PCKS Phase")
	elapsedPCKSParty = runTimedParty(func() {
		for _, pi := range P {
			pcks.GenShare(pi.tsk, tpk, encRes.Ciphertext.Value[1], pi.pcksShare)
		}
	}, len(P))

	pcksCombined := pcks.AllocateShare()
	encOut = bfv.NewCiphertext(params, 1)
	elapsedPCKSCloud = runTimed(func() {
		for _, pi := range P {
			pcks.AggregateShare(pi.pcksShare, pcksCombined, pcksCombined)
		}
		pcks.KeySwitch(encRes, pcksCombined, encOut)

	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedPCKSCloud, elapsedPCKSParty)

	return

}

func rkgphase(params bfv.Parameters, crs drlwe.CRS, P []*party) *rlwe.RelinearizationKey {

	l.Println("> RKG Phase")

	rkg := dbfv.NewRKGProtocol(params) // Relineariation key generation

	for _, pi := range P {
		pi.rlkEphemSk, pi.rkgShareOne, pi.rkgShareTwo = rkg.AllocateShare()
	}

	crp := rkg.SampleCRP(crs)

	elapsedRKGParty = runTimedParty(func() {
		for _, pi := range P {
			rkg.GenShareRoundOne(pi.sk, crp, pi.rlkEphemSk, pi.rkgShareOne)
		}
	}, len(P))

	_, rkgCombined1, rkgCombined2 := rkg.AllocateShare()

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

	rlk := bfv.NewRelinearizationKey(params, 1)
	elapsedRKGCloud += runTimed(func() {
		for _, pi := range P {
			rkg.AggregateShare(pi.rkgShareTwo, rkgCombined2, rkgCombined2)
		}
		rkg.GenRelinearizationKey(rkgCombined1, rkgCombined2, rlk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedRKGCloud, elapsedRKGParty)

	return rlk
}

func ckgphase(params bfv.Parameters, crs drlwe.CRS, P []*party) *rlwe.PublicKey {

	l.Println("> CKG Phase")

	ckg := dbfv.NewCKGProtocol(params) // Public key generation
	crp := ckg.SampleCRP(crs)          // for the public-key

	for _, pi := range P {
		pi.ckgShare = ckg.AllocateShare()
	}

	elapsedCKGParty = runTimedParty(func() {
		for _, pi := range P {
			ckg.GenShare(pi.sk, crp, pi.ckgShare)
		}
	}, len(P))

	ckgCombined := ckg.AllocateShare()

	pk := bfv.NewPublicKey(params)

	elapsedCKGCloud = runTimed(func() {
		for _, pi := range P {
			ckg.AggregateShare(pi.ckgShare, ckgCombined, ckgCombined)
		}
		ckg.GenPublicKey(ckgCombined, crp, pk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKGCloud, elapsedCKGParty)

	return pk
}
