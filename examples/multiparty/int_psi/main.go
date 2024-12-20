// This example demonstrates the use of the [multiparty] package to perform a basic N-party private set intersection (PSI) protocol with external receiver.
// This protocol relies on the N-out-of-N threshold variant of the BGV scheme.
// In a nutshell, the parties encode their sets as binary vectors, and a helper server compute (under encryption) the intersection as the compoenent-wise product of all vectors.
// The result is then re-encrypted towards the receiver's key.
// For more details about the PSI circuit example see the paper [Multiparty Homomorphic Encryption from Ring-Learning-With-Errors] by by Christian Mouchet, Juan Troncoso-Pastoriza, Jean-Philippe Bossuat, and Jean-Pierre Hubaux.
//
// To run the example, use the following command:
//
//	go run main.go N NGoRoutines
//
// where N is the number of parties (default:16) and NGoRoutines is the number of Go routines (default: 1) to use during the homomorphic evaluation.
// All parties are run in the same process.
//
// The example demonstrates the following steps:
//
//  1. Setup: The parties generate a collectice public encryption key and a collective relinearization key.
//  2. Inputs: Each party encrypts its input vector and send it to a helper server.
//  3. Evaluation: The helper server computes the multiplication of the input vectors and relinearizes the result.
//  4. Output: The helper server, with the help of the N parties, switches the encryption of the result to the target public key.
//  5. Decryption: The target party decrypts the result with its secret key.
//
// [Multiparty Homomorphic Encryption from Ring-Learning-With-Errors]: https://eprint.iacr.org/2020/304
// TODO: do we want a README version of this docstring ?
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

// party is a type for the parties' state in the protocol, to be kept accross the different phases.
type party struct {
	sk         *rlwe.SecretKey // secret key of the party
	rlkEphemSk *rlwe.SecretKey // ephemeral state to be used in the RKG protocol

	input []uint64 // the input of the party, encoding a set as a binary vector
}

var l = log.New(os.Stderr, "", 0)

func main() {

	// $go run main.go N NGoRoutine
	// N: number of parties
	// NGoRoutines: number of Go routines
	N := 16 // Default number of parties
	var err error
	if len(os.Args[1:]) >= 1 {
		N, err = strconv.Atoi(os.Args[1])
		check(err)

		if N < 2 || N > 128 {
			l.Fatal("N must be in the range [2, ..., 128]")
		}
	}

	NGoRoutine := 1 // Default number of Go routines
	if len(os.Args[1:]) >= 2 {
		NGoRoutine, err = strconv.Atoi(os.Args[2])
		check(err)

		if NGoRoutine < 1 {
			l.Fatal("NGoRoutine must be at least 1")
		}
	}

	// Creating encryption parameters from a default params with logN=14, logQP=438 with a plaintext modulus T=65537
	params, err := bgv.NewParametersFromLiteral(bgv.ParametersLiteral{
		LogN:             14,
		LogQ:             []int{56, 55, 55, 54, 54, 54},
		LogP:             []int{55, 55},
		PlaintextModulus: 65537,
	})
	check(err)

	// Creates a PRNG that will be used to sample the common reference string (crs)
	crs, err := sampling.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	check(err)

	// Generate some keys for the receiver (target party)
	tsk, tpk := rlwe.NewKeyGenerator(params).GenKeyPairNew()

	// Create the N input parties and generate their secret keys
	P := genparties(params, N)

	// Step 1: Setup of the collective public key and relinearization key

	pk := execCKGProtocol(params, crs, P) // generates the collective public key

	rlk := execRKGProtocol(params, crs, P) // generates the collective relinearization key

	evk := rlwe.NewMemEvaluationKeySet(rlk) // creates the evaluation key from the relinearization key

	l.Printf("\tdone (cloud: %s, party: %s)\n",
		elapsedRKGCloud, elapsedRKGParty)
	l.Printf("\tSetup done (cloud: %s, party: %s)\n",
		elapsedRKGCloud+elapsedCKGCloud, elapsedRKGParty+elapsedCKGParty)

	// Step 2: Each party encrypts its input vector
	expRes := genInputs(params, P) // generates the input vectors and the expected result

	encoder := bgv.NewEncoder(params)
	encInputs := inputPhase(params, P, pk, encoder) // encrypts the input vectors

	// Step 3: The helper server computes the multiplication of the input vectors and relinearizes the result
	encRes := evalPhase(params, NGoRoutine, encInputs, evk)

	// Step 4: The helper server switches the encryption of the result to the target public key
	encOut := execPCKSProtocol(params, tpk, encRes, P)

	// Step 5: The target party decrypts the result with its secret key
	l.Println("> ResulPlaintextModulus:")
	decryptor := rlwe.NewDecryptor(params, tsk)
	ptres := bgv.NewPlaintext(params, params.MaxLevel())
	elapsedDecParty := runTimed(func() {
		decryptor.Decrypt(encOut, ptres)
	})

	// Check the result
	res := make([]uint64, params.MaxSlots())
	err = encoder.Decode(ptres, res)
	check(err)

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

func genparties(params bgv.Parameters, N int) []party {

	// Create the parties and generates a secret key for each party
	P := make([]party, N)
	for i := range P {
		P[i].sk = rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	}

	return P
}

func execCKGProtocol(params bgv.Parameters, crs sampling.PRNG, P []party) *rlwe.PublicKey {

	l.Println("> PublicKeyGen Phase")

	// Creates a protocol type for the collective public key generation.
	// The type is stateless and can be used to generate as many public keys as needed.
	ckg := multiparty.NewPublicKeyGenProtocol(params)

	// Allocates the memory for the parties' shares in the protocol
	ckgShares := make([]multiparty.PublicKeyGenShare, len(P))
	for i := range ckgShares {
		ckgShares[i] = ckg.AllocateShare()
	}
	ckgCombined := ckg.AllocateShare() // Allocate the memory for the combined share

	// sample the common reference polynomial (crp) common reference string (crs)
	crp := ckg.SampleCRP(crs)

	// Generate the parties' shares
	elapsedCKGParty = runTimedParty(func() {
		for i, pi := range P {
			ckg.GenShare(pi.sk, crp, &ckgShares[i])
		}
	}, len(P))

	// Aggregate the parties' shares into a collective public key
	pk := rlwe.NewPublicKey(params)
	elapsedCKGCloud = runTimed(func() {
		// Aggregate the parties' shares into a combined share
		for i := range P {
			ckg.AggregateShares(ckgShares[i], ckgCombined, &ckgCombined)
		}

		// Generate the public key from the combined share
		ckg.GenPublicKey(ckgCombined, crp, pk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKGCloud, elapsedCKGParty)

	return pk
}

func execRKGProtocol(params bgv.Parameters, crs sampling.PRNG, P []party) *rlwe.RelinearizationKey {

	l.Println("> RelinearizationKeyGen Phase")

	// Creates a protocol type for the collective relinearization key generation.
	// The type is stateless and can be used to generate as many relinearization keys as needed.
	// The RKG protocol has two rounds.
	rkg := multiparty.NewRelinearizationKeyGenProtocol(params)

	// Allocates the memory for the parties' shares in the protocol
	rkgSharesRoundOne := make([]multiparty.RelinearizationKeyGenShare, len(P))
	rkgSharesRoundTwo := make([]multiparty.RelinearizationKeyGenShare, len(P))
	for i := range P {
		P[i].rlkEphemSk, rkgSharesRoundOne[i], rkgSharesRoundTwo[i] = rkg.AllocateShare()
		// the parties have a private ephemeral secret key in the RKGen protocol
	}
	// Allocate the memory for the combined public shares
	_, rkgCombined1, rkgCombined2 := rkg.AllocateShare()

	// Sample the common reference polynomial (crp) common reference string (crs)
	crp := rkg.SampleCRP(crs)

	// The parties generate their shares for round one
	elapsedRKGParty = runTimedParty(func() {
		for i, pi := range P {
			rkg.GenShareRoundOne(pi.sk, crp, pi.rlkEphemSk, &rkgSharesRoundOne[i])
		}
	}, len(P))

	// the helper aggregates the parties' shares for round one
	elapsedRKGCloud = runTimed(func() {
		for i := range P {
			rkg.AggregateShares(rkgSharesRoundOne[i], rkgCombined1, &rkgCombined1)
		}
	})

	// The parties generate their shares for round two
	elapsedRKGParty += runTimedParty(func() {
		for i, pi := range P {
			rkg.GenShareRoundTwo(pi.rlkEphemSk, pi.sk, rkgCombined1, &rkgSharesRoundTwo[i])
		}
	}, len(P))

	// the helper aggregates the parties' shares for round two and generates the relinearization key
	rlk := rlwe.NewRelinearizationKey(params)
	elapsedRKGCloud += runTimed(func() {
		for i := range P {
			rkg.AggregateShares(rkgSharesRoundTwo[i], rkgCombined2, &rkgCombined2)
		}
		rkg.GenRelinearizationKey(rkgCombined1, rkgCombined2, rlk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedRKGCloud, elapsedRKGParty)

	return rlk
}

func genInputs(params bgv.Parameters, P []party) (expRes []uint64) {

	// generate input vectors for the parties of max size
	expRes = make([]uint64, params.MaxSlots())
	for i := range expRes {
		expRes[i] = 1
	}

	for i := range P {
		P[i].input = make([]uint64, params.MaxSlots())
		for j := range P[i].input {
			if sampling.RandFloat64(0, 1) > 0.3 || j == 4 {
				P[i].input[j] = 1
			}
			expRes[j] *= P[i].input[j]
		}

	}

	return
}

func inputPhase(params bgv.Parameters, P []party, pk *rlwe.PublicKey, encoder *bgv.Encoder) (encInputs []*rlwe.Ciphertext) {

	// Allocate the memory for the encrypted input vectors
	encInputs = make([]*rlwe.Ciphertext, len(P))
	for i := range encInputs {
		encInputs[i] = bgv.NewCiphertext(params, 1, params.MaxLevel())
	}

	// Each party encrypts its input vector
	l.Println("> Encrypt Phase")
	encryptor := rlwe.NewEncryptor(params, pk)

	pt := bgv.NewPlaintext(params, params.MaxLevel())
	elapsedEncryptParty = runTimedParty(func() {
		for i, pi := range P {
			err := encoder.Encode(pi.input, pt)
			check(err)
			err = encryptor.Encrypt(pt, encInputs[i])
			check(err)
		}
	}, len(P))

	elapsedEncryptCloud = time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	return
}

func evalPhase(params bgv.Parameters, NGoRoutine int, encInputs []*rlwe.Ciphertext, evk rlwe.EvaluationKeySet) (encRes *rlwe.Ciphertext) {

	// The eval phase performs the multiplication as a balanced binary tree.
	// For each level of the tree, it performs the multiplications in parallel using at most NGoRoutine Go routines.

	// Allocate the memory for the encrypted result at each level of the tree
	encLvls := make([][]*rlwe.Ciphertext, 0)
	encLvls = append(encLvls, encInputs)
	for nLvl := len(encInputs) / 2; nLvl > 0; nLvl = nLvl >> 1 {
		encLvl := make([]*rlwe.Ciphertext, nLvl)
		for i := range encLvl {
			encLvl[i] = bgv.NewCiphertext(params, 2, params.MaxLevel())
		}
		encLvls = append(encLvls, encLvl)
	}
	encRes = encLvls[len(encLvls)-1][0]

	// Creates a evaluator for the multiplication, with the evaluation key
	evaluator := bgv.NewEvaluator(params, evk, true)

	// Split the task among the Go routines
	// A multTask is a task that multiplies two ciphertexts and relinearizes the result
	type multTask struct {
		wg              *sync.WaitGroup
		op1             *rlwe.Ciphertext
		opOut           *rlwe.Ciphertext
		res             *rlwe.Ciphertext
		elapsedmultTask time.Duration
	}
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
					err := evaluator.Mul(task.op1, task.opOut, task.res)
					check(err)
					// 2) Relinearization
					err = evaluator.Relinearize(task.res, task.res)
					check(err)

				})
				task.wg.Done()
			}
			workers.Done()
		}(i)
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

func execPCKSProtocol(params bgv.Parameters, tpk *rlwe.PublicKey, encRes *rlwe.Ciphertext, P []party) (encOut *rlwe.Ciphertext) {

	// Collective key switching from the collective secret key to
	// the target public key
	l.Println("> PublicKeySwitch Phase")

	// Creates a protocol type for the collective public key switch.
	// The type is stateless and can be used to generate as many public key switches as needed.
	pcks, err := multiparty.NewPublicKeySwitchProtocol(params, ring.DiscreteGaussian{Sigma: 1 << 30, Bound: 6 * (1 << 30)})
	check(err)

	// Allocates the memory for the parties' shares in the protocol
	pcksShares := make([]multiparty.PublicKeySwitchShare, len(P))
	for i := range P {
		pcksShares[i] = pcks.AllocateShare(params.MaxLevel())
	}
	// Allocates the memory for combined share
	pcksCombined := pcks.AllocateShare(params.MaxLevel())

	// Each party generates its share
	elapsedPCKSParty = runTimedParty(func() {
		for i, pi := range P {
			pcks.GenShare(pi.sk, tpk, encRes, &pcksShares[i])
		}
	}, len(P))

	// The helper server aggregates the parties' shares and combutes the output, re-encrpyted, ciphertext
	encOut = bgv.NewCiphertext(params, 1, params.MaxLevel())
	elapsedPCKSCloud = runTimed(func() {
		// Aggregate the parties' shares into a combined share
		for i := range P {
			err := pcks.AggregateShares(pcksShares[i], pcksCombined, &pcksCombined)
			check(err)
		}

		// Generate the output ciphertext
		pcks.KeySwitch(encRes, pcksCombined, encOut)
	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedPCKSCloud, elapsedPCKSParty)

	return
}

var (
	elapsedEncryptParty,
	elapsedEncryptCloud,
	elapsedCKGCloud,
	elapsedCKGParty,
	elapsedRKGCloud,
	elapsedRKGParty,
	elapsedPCKSCloud,
	elapsedPCKSParty,
	elapsedEvalCloudCPU,
	elapsedEvalCloud,
	elapsedEvalParty time.Duration
)

func check(err error) {
	if err != nil {
		l.Fatal(err)
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
