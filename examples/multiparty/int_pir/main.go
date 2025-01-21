// This example demonstrates the use of the [multiparty] package to perform a basic N-party private information retrival (PIR) protocol.
// This protocol relies on the t-out-of-N-threshold variant of the BGV scheme.
//
// In a nutshell, each party uploads an data row to a helper server, as an encrypted integer vector. The result is an encrypted database of N rows. At a later stage,
// a party queries a row in the database to the helper, without revealing the selector query index. The querier does so by encoding its query as
// a binary vector with a single 1 component corresponding to the row it wants to retrive, and send this encrypted vector to the helper.
// The helper can then compute the response (under encryption) by multiplying each row i of the database with the i-th query component, and summing
// the resulting vectors together. Finally, the result can be decrypted by any group of t parties.
//
// Thanks to the t-out-of-N-threshold scheme, only t parties need to be online at query time.
//
// For more details about the PIR circuit example see the paper [Multiparty Homomorphic Encryption from Ring-Learning-With-Errors] by Mouchet, Troncoso-Pastoriza, Bossuat, and Hubaux.
// For more details about the t-out-of-N-threshold scheme, see the paper [An Efficient Threshold Access-Structure for RLWE-Based Multiparty Homomorphic Encryption] by Mouchet, Bertrand and Hubaux.
//
// To run the example, use the following command:
//
//	go run main.go N T NGoRoutines
//
// where N is the number of parties (default:3) T is the threshold (default: 2) and NGoRoutines is the number of Go routines (default: 1) to use during the homomorphic evaluation.
// All parties are run in the same process.
//
// The example demonstrates the following steps:
//
//  1. Setup:
//     a. The parties generate threshold secret key with t-out-of-N access structure.
//     b. The parties generate a collective public encryption key, a relinearization key, and a set of rotation (or, galois) keys.
//  2. Database inputs: Each party encrypts its input row vector and send it to a helper server.
//  3. Query evaluation: A party requests a row in the database. The helper server computes the query output as described above.
//  4. Query output decryption: the helper, with the help of at least t parties, decrypts the query output.
//
// [Multiparty Homomorphic Encryption from Ring-Learning-With-Errors]: https://eprint.iacr.org/2020/304
// [An Efficient Threshold Access-Structure for RLWE-Based Multiparty Homomorphic Encryption]: https://eprint.iacr.org/2022/780
package main

import (
	"errors"
	"log"
	"math/rand"
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
	multiparty.Combiner

	sk         *rlwe.SecretKey // secret key of the party
	tsk        multiparty.ShamirSecretShare
	rlkEphemSk *rlwe.SecretKey // ephemeral state to be used in the RKG protocol

	shamirPt multiparty.ShamirPublicPoint

	input []uint64 // the input of the party, encoding a set as a binary vector
}

// getOnlineParties is a utility function that returns t random parties from a list of parties.
// This simulates a dynamic setting where the system rely on any t parties to be online at query time to
// execute the various protocols.
func getOnlineParties(t int, parties []party) []party {
	if t > len(parties) {
		check(errors.New("t must be less than the number of parties"))
	}
	// randomizes a subset of t parties
	onlineParties := make([]party, len(parties))
	copy(onlineParties, parties)
	rand.Shuffle(len(onlineParties), func(i, j int) { onlineParties[i], onlineParties[j] = onlineParties[j], onlineParties[i] })
	return onlineParties[:t]
}

// getShamirPoints is a utility function that returns the Shamir public points of a group of parties.
func getShamirPoints(parties []party) []multiparty.ShamirPublicPoint {
	shamirPoints := make([]multiparty.ShamirPublicPoint, len(parties))
	for i := range parties {
		shamirPoints[i] = parties[i].shamirPt
	}
	return shamirPoints
}

var l = log.New(os.Stderr, "", 0)

func main() {

	// Parse command line arguments

	N := 3 // Default number of parties
	var err error
	if len(os.Args[1:]) >= 1 {
		N, err = strconv.Atoi(os.Args[1])
		check(err)

		if N < 3 || N > 128 {
			l.Fatal("N must be in the range [3, 128]")
		}
	}

	t := 2 // Default Threshold
	if len(os.Args[1:]) >= 2 {
		t, err = strconv.Atoi(os.Args[2])
		check(err)

		if t < 2 || t >= N {
			l.Fatal("T must be in the range [2, N-1]")
		}
	}

	nGoRoutine := 1 // Default number of Go routines
	if len(os.Args[1:]) >= 3 {
		nGoRoutine, err = strconv.Atoi(os.Args[2])
		check(err)

		if nGoRoutine < 1 {
			l.Fatal("NGoRoutine must be at least 1")
		}
	}

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

	// The circuit relies on rotations/automorphisms for computing an inner-sum-like operation.
	// This obtains the corresponding Galois elements.
	galEls := append(params.GaloisElementsForInnerSum(1, params.N()>>1), params.GaloisElementForRowRotation())

	// Creates a PRNG that will be used to sample the common reference string (crs)
	crs, err := sampling.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}

	// Create the N input parties and generate their secret keys and private inputs
	P := genparties(params, N, t)

	l.Println("> Setup phase") // TODO: improve the logging output

	// Step 1.a: Generation of the threshold secret key

	thresholdizer := multiparty.NewThresholdizer(params.Parameters)

	// Allocates the memory for the parties' shares in the protocol
	tskShares := make([][]multiparty.ShamirSecretShare, N)
	for i := range P {
		tskShares[i] = make([]multiparty.ShamirSecretShare, N)
		for j := range P {
			tskShares[i][j] = thresholdizer.AllocateThresholdSecretShare()
		}
	}

	l.Println("> Threshold Secret Key Generation")
	// The parties generate their shares for the threshold secret key
	elapsedThresholdParty := runTimedParty(func() {
		for i, pi := range P {
			pol, err := thresholdizer.GenShamirPolynomial(t, pi.sk)
			check(err)
			for j, pj := range P {
				thresholdizer.GenShamirSecretShare(pj.shamirPt, pol, &tskShares[i][j])
			}
		}
	}, N)

	// Each party aggregates the shares it has received from the other parties
	elapsedThresholdParty += runTimedParty(func() {
		for i := range P {
			P[i].tsk = thresholdizer.AllocateThresholdSecretShare()
			for j := range P {
				err := thresholdizer.AggregateShares(tskShares[j][i], P[i].tsk, &P[i].tsk) // TODO wierd pointer arguments
				check(err)
			}
		}
	}, N)

	l.Printf("\tdone (cloud: %s, party: %s)\n", time.Duration(0), elapsedThresholdParty)

	// Step 1.b: Setup of the collective public encryption and evaluation keys
	// Note 1: because we are using the t-out-of-N-threshold scheme, generating a key only requires t parties to be online.
	// Note 2: the above note means that the workload could be balanced between the online parties whenever more than t parties are online.

	pk := execCKGProtocol(params, crs, getOnlineParties(t, P)) // Collective public key generation

	rlk := execRKGProtocol(params, crs, getOnlineParties(t, P)) // Collective RelinearizationKey generation

	galKeys := execGTGProtocol(params, crs, galEls, getOnlineParties(t, P)) // Collective GaloisKeys generation

	// Creates the evaluation key set from the rlk and the galKeys
	evk := rlwe.NewMemEvaluationKeySet(rlk, galKeys...)

	l.Printf("\tSetup done (cloud: %s, party: %s)\n",
		elapsedCKGCloud+elapsedRKGCloud+elapsedGKGCloud,
		elapsedCKGParty+elapsedRKGParty+elapsedGKGParty)

	// Step 2: Database inputs

	// Pre-allocates the encrypted database: one ciphertext per party
	encInputs := make([]*rlwe.Ciphertext, N)
	for i := range encInputs {
		encInputs[i] = bgv.NewCiphertext(params, 1, params.MaxLevel())
	}

	// Each party encrypts its input row under the collective public key
	l.Println("> Database input phase")
	encoder := bgv.NewEncoder(params)
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

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	// Step 3: Query evaluation

	l.Println("> Query evaluation phase")

	queryIndex := 2 // Index of the ciphertext to retrieve.
	querier := P[0] // Party performing the query
	others := P[1:] // Other parties

	// Creates a query ciphertext from the query index
	encQuery := genQuery(params, queryIndex, encoder, encryptor)

	// Executes the requests
	encResult := execRequest(params, nGoRoutine, encQuery, encInputs, evk)

	// Step 4: Query output decryption

	// The helper, with the help of at least t parties, performs a re-encryption of the result towards the querier
	participants := getOnlineParties(t-1, others)
	encOut := execCKSProtocol(params, participants, querier, encResult)

	// The querier decrypts the final result
	sk := rlwe.NewSecretKey(params)
	err = querier.Combiner.GenAdditiveShare(append(getShamirPoints(participants), querier.shamirPt), querier.shamirPt, querier.tsk, sk)
	check(err)

	decryptor := rlwe.NewDecryptor(params, sk)
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
		elapsedCKGParty+elapsedRKGParty+elapsedGKGParty+elapsedEncryptParty+elapsedRequestParty+elapsedCKSParty+elapsedDecParty)
}

func genparties(params bgv.Parameters, N, t int) []party {

	P := make([]party, N)
	kgen := rlwe.NewKeyGenerator(params)
	for i := range P {
		P[i].shamirPt = multiparty.ShamirPublicPoint(i + 1)

		P[i].sk = kgen.GenSecretKeyNew()

		P[i].input = make([]uint64, params.N())
		for j := range P[i].input {
			P[i].input[j] = uint64(i)
		}
	}

	shamirPts := getShamirPoints(P)
	for i := range P {
		P[i].Combiner = multiparty.NewCombiner(params.Parameters, P[i].shamirPt, shamirPts, t) // TODO: NewCombiner takes interface
	}

	return P
}

func execCKGProtocol(params bgv.Parameters, crs sampling.PRNG, participants []party) *rlwe.PublicKey {

	l.Println("> PublicKeyGen Phase")

	// Creates a protocol type for the collective public key generation.
	// The type is stateless and can be used to generate as many public keys as needed.
	ckg := multiparty.NewPublicKeyGenProtocol(params)

	// Allocates the memory for the parties' shares in the protocol
	ckgShares := make([]multiparty.PublicKeyGenShare, len(participants))
	tsks := make([]*rlwe.SecretKey, len(participants))
	for i := range ckgShares {
		ckgShares[i] = ckg.AllocateShare()  // the public CKG shares
		tsks[i] = rlwe.NewSecretKey(params) // the t-out-of-t secret keys for group P
	}
	ckgCombined := ckg.AllocateShare() // Allocate the memory for the combined share

	// sample the common reference polynomial (crp) common reference string (crs)
	crp := ckg.SampleCRP(crs)

	// Generate the parties' shares
	elapsedCKGParty = runTimedParty(func() {
		for i, pi := range participants {
			// Generate the t-out-of-t secret key of the party within the group of participants
			err := pi.Combiner.GenAdditiveShare(getShamirPoints(participants), pi.shamirPt, pi.tsk, tsks[i]) // TODO: discuss returning the key directly
			check(err)

			// Generate the public key share of the party from the t-out-of-t secret key
			ckg.GenShare(tsks[i], crp, &ckgShares[i])
		}
	}, len(participants))

	// Aggregate the parties' shares into a collective public key
	pk := rlwe.NewPublicKey(params)
	elapsedCKGCloud = runTimed(func() {
		// Aggregate the parties' shares into a combined share
		for i := range participants {
			ckg.AggregateShares(ckgShares[i], ckgCombined, &ckgCombined)
		}

		// Generate the public key from the combined share
		ckg.GenPublicKey(ckgCombined, crp, pk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKGCloud, elapsedCKGParty)

	return pk
}

func execRKGProtocol(params bgv.Parameters, crs sampling.PRNG, participants []party) *rlwe.RelinearizationKey {

	l.Println("> RelinearizationKeyGen Phase")

	// Creates a protocol type for the collective relinearization key generation.
	// The type is stateless and can be used to generate as many relinearization keys as needed.
	// The RKG protocol has two rounds. Because the ephemeral secret key is not re-shared,
	// the same set of participants must participate to the two rounds.
	rkg := multiparty.NewRelinearizationKeyGenProtocol(params)

	// Allocates the memory for the parties' shares in the protocol
	rkgSharesRoundOne := make([]multiparty.RelinearizationKeyGenShare, len(participants))
	rkgSharesRoundTwo := make([]multiparty.RelinearizationKeyGenShare, len(participants))
	tsks := make([]*rlwe.SecretKey, len(participants))
	for i := range participants {
		// the parties have a private ephemeral secret key in the RKGen protocol
		participants[i].rlkEphemSk, rkgSharesRoundOne[i], rkgSharesRoundTwo[i] = rkg.AllocateShare()
		tsks[i] = rlwe.NewSecretKey(params)
	}
	// Allocate the memory for the combined public shares
	_, rkgCombined1, rkgCombined2 := rkg.AllocateShare()

	// Sample the common reference polynomial (crp) common reference string (crs)
	crp := rkg.SampleCRP(crs)

	// The parties generate their shares for round one
	elapsedRKGParty = runTimedParty(func() {
		for i, pi := range participants {

			// Generate the t-out-of-t secret key of the party within the group of participants
			err := pi.Combiner.GenAdditiveShare(getShamirPoints(participants), pi.shamirPt, pi.tsk, tsks[i])
			check(err)

			// Generate the shares for round one from the t-out-of-t secret key
			rkg.GenShareRoundOne(tsks[i], crp, pi.rlkEphemSk, &rkgSharesRoundOne[i])
		}
	}, len(participants))

	// the helper aggregates the parties' shares for round one
	elapsedRKGCloud = runTimed(func() {
		for i := range participants {
			rkg.AggregateShares(rkgSharesRoundOne[i], rkgCombined1, &rkgCombined1)
		}
	})

	// The parties generate their shares for round two
	elapsedRKGParty += runTimedParty(func() {
		for i, pi := range participants {
			// Generate the shares for round two from the t-out-of-t secret key
			// Note: the same set of participants must participate to the two rounds, so the same tsk is used.
			rkg.GenShareRoundTwo(pi.rlkEphemSk, tsks[i], rkgCombined1, &rkgSharesRoundTwo[i])
		}
	}, len(participants))

	// the helper aggregates the parties' shares for round two and generates the relinearization key
	rlk := rlwe.NewRelinearizationKey(params)
	elapsedRKGCloud += runTimed(func() {
		for i := range participants {
			rkg.AggregateShares(rkgSharesRoundTwo[i], rkgCombined2, &rkgCombined2)
		}
		rkg.GenRelinearizationKey(rkgCombined1, rkgCombined2, rlk)
	})

	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedRKGCloud, elapsedRKGParty)

	return rlk
}

func execGTGProtocol(params bgv.Parameters, crs sampling.PRNG, galEls []uint64, participants []party) (galKeys []*rlwe.GaloisKey) {

	l.Println("> GKG Phase")

	// Creates a protocol type for the collective galois key generation.
	// The type is stateless and can be used to generate as many galois keys as needed.
	gkg := multiparty.NewGaloisKeyGenProtocol(params) // Rotation keys generation

	// Allocates the memory for the parties' shares in the protocol
	gkgShares := make([]multiparty.GaloisKeyGenShare, len(participants))
	tsks := make([]*rlwe.SecretKey, len(participants))
	for i := range participants {
		gkgShares[i] = gkg.AllocateShare()
		tsks[i] = rlwe.NewSecretKey(params)
	}

	// Allocate a slice for storing the output keys
	galKeys = make([]*rlwe.GaloisKey, len(galEls))

	// Runs the GKG protocol for each required Galois key
	// Note: this demo re-uses the allocated shares for each execution.
	for j, galEl := range galEls {

		// Sample the common reference polynomial (crp) common reference string (crs)
		crp := gkg.SampleCRP(crs)

		// The parties generate their shares for the Galois key generation protocol
		elapsedGKGParty += runTimedParty(func() {
			for i, pi := range participants {
				// Generate the t-out-of-t secret key of the party within the group of participants
				err := pi.Combiner.GenAdditiveShare(getShamirPoints(participants), pi.shamirPt, pi.tsk, tsks[i])
				check(err)

				// Generate the shares for the Galois key generation protocol from the t-out-of-t secret key
				err = gkg.GenShare(tsks[i], galEl, crp, &gkgShares[i])
				check(err)
			}

		}, len(participants))

		// The helper aggregates the parties' shares and generates the Galois key
		elapsedGKGCloud += runTimed(func() {

			gkgShareCombined := gkg.AllocateShare() // Allocate the memory for the combined share
			gkgShareCombined.GaloisElement = galEl
			for i := range participants {
				err := gkg.AggregateShares(gkgShares[i], gkgShareCombined, &gkgShareCombined)
				check(err)
			}

			galKeys[j] = rlwe.NewGaloisKey(params)

			if err := gkg.GenGaloisKey(gkgShareCombined, crp, galKeys[j]); err != nil {
				panic(err)
			}
		})
	}
	l.Printf("\tdone (cloud: %s, party %s)\n", elapsedGKGCloud, elapsedGKGParty)

	return
}

func genQuery(params bgv.Parameters, queryIndex int, encoder *bgv.Encoder, encryptor *rlwe.Encryptor) *rlwe.Ciphertext {

	l.Println("> QueryGen Phase")

	// Creates a query vector from the query index
	queryCoeffs := make([]uint64, params.N())
	queryCoeffs[queryIndex] = 1

	// Encrypts the query vector
	query := bgv.NewPlaintext(params, params.MaxLevel())
	var encQuery *rlwe.Ciphertext
	elapsedRequestParty += runTimed(func() {

		err := encoder.Encode(queryCoeffs, query)
		check(err)

		encQuery, err = encryptor.EncryptNew(query)
		check(err)
	})

	return encQuery
}

func execRequest(params bgv.Parameters, NGoRoutine int, encQuery *rlwe.Ciphertext, encInputs []*rlwe.Ciphertext, evk rlwe.EvaluationKeySet) *rlwe.Ciphertext {

	// First, pre-compute the some plaintext masks for the query evaluation as:
	// plainmask[i] = encode([0, ..., 0, 1, 0, ..., 0])  (zero with a 1 at the i-th position).
	// In practice, the masks are pre-computed and reused accross queries.
	encoder := bgv.NewEncoder(params)
	plainMask := make([]*rlwe.Plaintext, len(encInputs))
	for i := range plainMask {
		maskCoeffs := make([]uint64, params.N())
		maskCoeffs[i] = 1
		plainMask[i] = bgv.NewPlaintext(params, params.MaxLevel())
		if err := encoder.Encode(maskCoeffs, plainMask[i]); err != nil {
			panic(err)
		}
	}

	l.Println("> Request Phase")

	// Buffer for the intermediate computation done by the helper
	encPartial := make([]*rlwe.Ciphertext, len(encInputs))
	for i := range encPartial {
		encPartial[i] = bgv.NewCiphertext(params, 2, params.MaxLevel())
	}

	// Creates an evaluator for the homomorphic evaluation
	evaluator := bgv.NewEvaluator(params, evk)

	// Split the task among the Go routines

	// maskTask is a type for the task to be executed by the Go routines
	// The task computes the multiplication of the query with a mask, and the multiplication of the result with a row of the database.
	type maskTask struct {
		query           *rlwe.Ciphertext
		mask            *rlwe.Plaintext
		row             *rlwe.Ciphertext
		res             *rlwe.Ciphertext
		elapsedmaskTask time.Duration
	}
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
			workers.Done()
		}(i)
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
		workers.Wait() // Wait for all the workers to finish
	})

	// collects the elapsed time for each task
	for _, t := range taskList {
		elapsedRequestCloudCPU += t.elapsedmaskTask
	}

	// Creates ciphertexts to store the final result
	resultDeg2 := bgv.NewCiphertext(params, 2, params.MaxLevel()) // to receive the sum of the partial results.
	result := bgv.NewCiphertext(params, 1, params.MaxLevel())     // to receive the relinearized final result

	// Summation of all the partial result among the different Go routines
	// The sum is computed over the degree-2 ciphertexts from the previous step. Then, the result is relinearized.
	// This avoids performing N relinearizations.
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

func execCKSProtocol(params bgv.Parameters, participants []party, receiver party, ctIn *rlwe.Ciphertext) *rlwe.Ciphertext {

	l.Println("> KeySwitch Phase")

	// Creates a protocol type for the collective key-switching protocol, with smudging distribution parameter of 2^30.
	// The type is stateless and can be used to generate as many key-switching keys as needed.
	cks, err := multiparty.NewKeySwitchProtocol(params, ring.DiscreteGaussian{Sigma: 1 << 30, Bound: 6 * (1 << 30)})
	check(err)

	// Allocates the memory for the parties' shares in the protocol
	cksShares := make([]multiparty.KeySwitchShare, len(participants))
	tsks := make([]*rlwe.SecretKey, len(participants))
	for i := range participants {
		cksShares[i] = cks.AllocateShare(params.MaxLevel()) // Allocate the memory for the public share
		tsks[i] = rlwe.NewSecretKey(params)                 // Allocate the memory for the t-out-of-t secret key
	}
	cksCombined := cks.AllocateShare(params.MaxLevel()) // Allocate the memory for the combined share

	// To generate a re-encryption of the result ciphertexts towards the querier,
	// each party except for the receiver generates a key-switching share towards
	// secret-key zero (i.e., a decryption share).
	zero := rlwe.NewSecretKey(params)

	// The parties (except the receiver) generate their shares for the key-switching protocol
	elapsedCKSParty = runTimedParty(func() {
		for i, pi := range participants {

			// Generate the t-out-of-t secret key with the reciever and t-1 other parties
			err := pi.Combiner.GenAdditiveShare(append(getShamirPoints(participants), receiver.shamirPt), pi.shamirPt, pi.tsk, tsks[i])
			check(err)

			// Generate the key-switching share with the t-out-of-t secret key
			cks.GenShare(tsks[i], zero, ctIn, &cksShares[i])
		}
	}, len(participants))

	// The helper aggregates the parties' shares and generates the key-switching key
	ctOut := bgv.NewCiphertext(params, 1, params.MaxLevel())
	elapsedCKSCloud = runTimed(func() {
		// Aggregate the parties' shares into a combined share
		for i := range participants {
			err := cks.AggregateShares(cksShares[i], cksCombined, &cksCombined)
			check(err)
		}
		// Generate the re-encryption from the combined share
		cks.KeySwitch(ctIn, cksCombined, ctOut)
	})
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKSCloud, elapsedCKSParty)

	return ctOut
}

var (
	elapsedCKGCloud        time.Duration
	elapsedCKGParty        time.Duration
	elapsedRKGCloud        time.Duration
	elapsedRKGParty        time.Duration
	elapsedGKGCloud        time.Duration
	elapsedGKGParty        time.Duration
	elapsedCKSCloud        time.Duration
	elapsedCKSParty        time.Duration
	elapsedEncryptCloud    time.Duration
	elapsedRequestParty    time.Duration
	elapsedRequestCloud    time.Duration
	elapsedRequestCloudCPU time.Duration
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
