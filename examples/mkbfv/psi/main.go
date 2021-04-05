package main

import (
	"log"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkbfv"
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

type multTask struct {
	wg              *sync.WaitGroup
	op1             *mkbfv.MKCiphertext
	op2             *mkbfv.MKCiphertext
	res             *mkbfv.MKCiphertext
	elapsedmultTask time.Duration
}

var elapsedEncryptParty time.Duration
var elapsedEncryptCloud time.Duration
var elapsedDecryptParty time.Duration
var elapsedEvalCloudCPU time.Duration
var elapsedEvalCloud time.Duration
var elapsedEvalParty time.Duration

func main() {

	l := log.New(os.Stderr, "", 0)

	// retrieve parameters

	NParties := 8 // Default number of parties
	var err error
	if len(os.Args[1:]) >= 1 {

		NParties, err = strconv.Atoi(os.Args[1])
		check(err)
	}

	NGoRoutine := 1 // Default number of Go routines
	if len(os.Args[1:]) >= 2 {
		NGoRoutine, err = strconv.Atoi(os.Args[2])
		check(err)
	}

	// gen parameters and CRS
	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})

	if err != nil {
		panic(err)
	}

	// Use the defaultparams logN=14, logQP=438 with a plaintext modulus T=65537
	params := bfv.DefaultParams[bfv.PN14QP438].WithT(65537)

	sigma := 6.0

	crs := mkbfv.GenCommonPublicParam(params, prng)

	// Setup each participant
	participants := genParticipants(crs, params, sigma, NParties)

	// generate values
	expected, values := genInputs(params, participants)

	// encrypt all values
	encrypted := encPhase(participants, values)

	// gen keys for relin
	evalKeys := getEvaluationKeys(participants)
	pubKeys := getPublicKeys(participants)

	// compute circuit
	encryptedResult := evalPhase(params, NGoRoutine, encrypted, evalKeys, pubKeys)

	var decrypted []uint64

	elapsedDecryptParty = runTimed(func() {

		// gen partial decryption
		partialDecryptions := genPartialDecryption(participants, encryptedResult)
		// decrypt
		decrypted = participants[0].Decrypt(encryptedResult, partialDecryptions)
	})

	// check result

	for i, v := range expected {

		if v != decrypted[i] {
			l.Println("\tincorrect")
			return
		}
	}

	l.Println("\tcorrect")
	l.Printf("> Finished (total cloud: %s, total party: %s)\n",
		elapsedEncryptCloud+elapsedEvalCloud,
		elapsedEncryptParty+elapsedEvalParty+elapsedDecryptParty)

}

func genParticipants(crs *mkbfv.MKDecomposedPoly, params *bfv.Parameters, sigmaSmudging float64, nbrParties int) []mkbfv.MKParticipant {

	res := make([]mkbfv.MKParticipant, nbrParties)

	for i := 0; i < nbrParties; i++ {

		res[i] = mkbfv.NewParticipant(params, sigmaSmudging, crs)

	}

	return res
}

func genPartialDecryption(participants []mkbfv.MKParticipant, cipher *mkbfv.MKCiphertext) []*ring.Poly {

	res := make([]*ring.Poly, len(participants))

	for i, p := range participants {
		res[i] = p.GetPartialDecryption(cipher)
	}

	return res
}

func getEvaluationKeys(participants []mkbfv.MKParticipant) []*mkbfv.MKEvaluationKey {

	res := make([]*mkbfv.MKEvaluationKey, len(participants))

	for i, p := range participants {
		res[i] = p.GetEvaluationKey()
	}

	return res
}

func getPublicKeys(participants []mkbfv.MKParticipant) []*mkbfv.MKPublicKey {

	res := make([]*mkbfv.MKPublicKey, len(participants))

	for i, p := range participants {
		res[i] = p.GetPublicKey()
	}

	return res
}

func encPhase(P []mkbfv.MKParticipant, values [][]uint64) (encInputs []*mkbfv.MKCiphertext) {

	l := log.New(os.Stderr, "", 0)

	encInputs = make([]*mkbfv.MKCiphertext, len(P))

	// Each party encrypts its input vector
	l.Println("> Encrypt Phase")

	elapsedEncryptParty = runTimedParty(func() {
		for i, pi := range P {
			encInputs[i] = pi.Encrypt(values[i])
		}
	}, len(P))

	elapsedEncryptCloud = time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	return
}

func evalPhase(params *bfv.Parameters, NGoRoutine int, encInputs []*mkbfv.MKCiphertext, rlk []*mkbfv.MKEvaluationKey, pubKeys []*mkbfv.MKPublicKey) (encRes *mkbfv.MKCiphertext) {

	l := log.New(os.Stderr, "", 0)

	encLvls := make([][]*mkbfv.MKCiphertext, 0)
	encLvls = append(encLvls, encInputs)
	for nLvl := len(encInputs) / 2; nLvl > 0; nLvl = nLvl >> 1 {
		encLvl := make([]*mkbfv.MKCiphertext, nLvl)
		encLvls = append(encLvls, encLvl)
	}
	encRes = encLvls[len(encLvls)-1][0]

	evaluator := mkbfv.NewMKEvaluator(params)

	// Split the task among the Go routines
	tasks := make(chan *multTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)
	//l.Println("> Spawning", NGoRoutine, "evaluator goroutine")
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {

			for task := range tasks {
				task.elapsedmultTask = runTimed(func() {
					// 1) Multiplication and Relinearization of two input vectors
					task.res = evaluator.MultRelinDynamic(task.op1, task.op2, rlk, pubKeys)
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

func genInputs(params *bfv.Parameters, P []mkbfv.MKParticipant) (expRes []uint64, participantValues [][]uint64) {

	expRes = make([]uint64, params.N())
	for i := range expRes {
		expRes[i] = 1
	}

	participantValues = make([][]uint64, len(P))

	for j := range P {

		input := make([]uint64, params.N())
		for i := range input {
			if utils.RandFloat64(0, 1) > 0.3 || i == 4 {
				input[i] = 1
			}
			expRes[i] *= input[i]
		}

		participantValues[j] = input
	}

	return
}
