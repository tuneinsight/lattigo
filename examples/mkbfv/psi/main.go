package main

import (
	"log"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkbfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
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

	paramsDef := bfv.PN14QP438
	paramsDef.T = 65537
	p, err := bfv.NewParametersFromLiteral(paramsDef)
	if err != nil {
		panic(err)
	}
	params := &p

	sigma := 6.0

	crs := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)

	// Setup each participant
	participants := genParticipants(crs, params, sigma, NParties)

	// generate values
	expected, values := genInputs(params, participants)

	// encrypt all values
	encrypted := encPhase(participants, values)

	// gen keys for relin
	evalKeys := getEvaluationKeys(participants)
	pubKeys := getPublicKeys(participants)

	// create Evaluator and IDs
	evaluator := mkbfv.NewMKEvaluator(params)

	ids := make([]uint64, NParties)
	for i := uint64(0); i < uint64(NParties); i++ {
		ids[i] = i
	}

	ciphers := evaluator.ConvertToMKCiphertext(encrypted, ids)

	// compute circuit
	encryptedResult := evalPhase(params, NGoRoutine, ciphers, evalKeys, pubKeys, evaluator)

	resultBFV := evaluator.ConvertToBFVCiphertext(encryptedResult)

	var decrypted []uint64

	elapsedDecryptParty = runTimed(func() {

		// gen partial decryption
		partialDecryptions := genPartialDecryption(participants, resultBFV)
		// decrypt
		decrypted = participants[0].Decrypt(resultBFV[0], partialDecryptions)
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

func genParticipants(crs *mkrlwe.MKDecomposedPoly, params *bfv.Parameters, sigmaSmudging float64, nbrParties int) []participant {

	res := make([]participant, nbrParties)

	for i := 0; i < nbrParties; i++ {

		res[i] = newParticipant(params, sigmaSmudging, crs)

	}

	return res
}

func genPartialDecryption(participants []participant, cipher []*bfv.Ciphertext) []*ring.Poly {

	res := make([]*ring.Poly, len(participants))

	for i, p := range participants {
		res[i] = p.GetPartialDecryption(cipher[i])
	}

	return res
}

func getEvaluationKeys(participants []participant) []*mkrlwe.MKEvaluationKey {

	res := make([]*mkrlwe.MKEvaluationKey, len(participants))

	for i, p := range participants {
		res[i] = p.GetEvaluationKey()
		res[i].PeerID = uint64(i)
	}

	return res
}

func getPublicKeys(participants []participant) []*mkrlwe.MKPublicKey {

	res := make([]*mkrlwe.MKPublicKey, len(participants))

	for i, p := range participants {
		res[i] = p.GetPublicKey()
		res[i].PeerID = uint64(i)
	}

	return res
}

func encPhase(P []participant, values [][]uint64) (encInputs []*bfv.Ciphertext) {

	l := log.New(os.Stderr, "", 0)

	encInputs = make([]*bfv.Ciphertext, len(P))

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

func evalPhase(params *bfv.Parameters, NGoRoutine int, encInputs []*mkbfv.MKCiphertext, rlk []*mkrlwe.MKEvaluationKey, pubKeys []*mkrlwe.MKPublicKey, evaluator mkbfv.MKEvaluator) (encRes *mkbfv.MKCiphertext) {

	l := log.New(os.Stderr, "", 0)

	encLvls := make([][]*mkbfv.MKCiphertext, 0)
	encLvls = append(encLvls, encInputs)
	for nLvl := len(encInputs) / 2; nLvl > 0; nLvl = nLvl >> 1 {
		encLvl := make([]*mkbfv.MKCiphertext, nLvl)
		for i := range encLvl {
			encLvl[i] = new(mkbfv.MKCiphertext)
		}
		encLvls = append(encLvls, encLvl)
	}
	encRes = encLvls[len(encLvls)-1][0]

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
					tmpRes := evaluator.Mul(task.op1, task.op2)
					evaluator.RelinInPlace(tmpRes, getRelinKeyForParticipants(rlk, tmpRes.PeerID), getPublicKeyForParticipants(pubKeys, tmpRes.PeerID))
					task.res.PeerID = tmpRes.PeerID
					task.res.Ciphertexts = tmpRes.Ciphertexts
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

func genInputs(params *bfv.Parameters, P []participant) (expRes []uint64, participantValues [][]uint64) {

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

func getRelinKeyForParticipants(rlk []*mkrlwe.MKEvaluationKey, peerID []uint64) []*mkrlwe.MKEvaluationKey {

	res := make([]*mkrlwe.MKEvaluationKey, 0)

	for _, v := range rlk {
		if mkrlwe.Contains(peerID, v.PeerID) >= 0 {
			res = append(res, v)
		}
	}

	return res
}

func getPublicKeyForParticipants(pk []*mkrlwe.MKPublicKey, peerID []uint64) []*mkrlwe.MKPublicKey {

	res := make([]*mkrlwe.MKPublicKey, 0)

	for _, v := range pk {
		if mkrlwe.Contains(peerID, v.PeerID) >= 0 {
			res = append(res, v)
		}
	}

	return res
}

//--------------------------------Participants interface for tests------------------------------------------

// MKParticipant is a type for participants in a multy key bfv scheme
type participant interface {
	GetEvaluationKey() *mkrlwe.MKEvaluationKey
	GetPublicKey() *mkrlwe.MKPublicKey
	Encrypt(values []uint64) *bfv.Ciphertext
	Decrypt(cipher *bfv.Ciphertext, partialDecryptions []*ring.Poly) []uint64
	GetPartialDecryption(ciphertext *bfv.Ciphertext) *ring.Poly
}

type mkParticipant struct {
	encryptor     mkbfv.MKEncryptor
	decryptor     mkrlwe.MKDecryptor
	params        *bfv.Parameters
	keys          *mkrlwe.MKKeys
	encoder       bfv.Encoder
	ringQ         *ring.Ring
	sigmaSmudging float64
}

// GetEvaluationKey returns the evaluation key of the participant
func (participant *mkParticipant) GetEvaluationKey() *mkrlwe.MKEvaluationKey {
	return participant.keys.EvalKey
}

// GetPublicKey returns the publik key of the participant
func (participant *mkParticipant) GetPublicKey() *mkrlwe.MKPublicKey {
	return participant.keys.PublicKey
}

// Encrypt constructs a ciphertext from the given values
func (participant *mkParticipant) Encrypt(values []uint64) *bfv.Ciphertext {
	if values == nil || len(values) <= 0 {
		panic("Cannot encrypt uninitialized or empty values")
	}

	pt := newPlaintext(values, participant.encoder, participant.params)

	return participant.encryptor.Encrypt(pt)
}

// Decrypt returns the decryption of the ciphertext given the partial decryption
func (participant *mkParticipant) Decrypt(cipher *bfv.Ciphertext, partialDecryptions []*ring.Poly) []uint64 {

	if cipher == nil || cipher.Degree() != 1 {
		panic("Cannot decrypt uninitialized ciphertext or cipher of degree greater than 1")
	}

	if partialDecryptions == nil || len(partialDecryptions) < 1 {
		panic("Decryption necessitates at least one partialy decrypted ciphertext")
	}

	decrypted := participant.decryptor.MergeDec(cipher.Element, uint64(len(participant.ringQ.Modulus)-1), partialDecryptions)

	pt := bfv.NewPlaintext(*participant.params)
	pt.SetValue(decrypted)

	return participant.encoder.DecodeUintNew(pt)
}

// GetPartialDecryption returns the partial decryption of an element in the ciphertext
// this function should only be used by participants that were involved in the given ciphertext
func (participant *mkParticipant) GetPartialDecryption(cipher *bfv.Ciphertext) *ring.Poly {

	return participant.decryptor.PartDec(cipher.Element, uint64(len(participant.ringQ.Modulus)-1), participant.keys.SecretKey, participant.sigmaSmudging)
}

// newParticipant creates a participant for the multi key bfv scheme
// the bfv parameters as well as the standard deviation used for partial decryption must be provided
func newParticipant(params *bfv.Parameters, sigmaSmudging float64, crs *mkrlwe.MKDecomposedPoly) participant {

	if crs == nil || params == nil {
		panic("Uninitialized parameters. Cannot create new participant")
	}

	if sigmaSmudging < params.Sigma() {
		panic("Sigma must be at least greater than the standard deviation of the gaussian distribution")
	}

	if len(crs.Poly) != int(params.Beta()) {
		panic("CRS must be the same dimention as returned by the function bfv.Parameters.Beta()")
	}

	keys := mkrlwe.KeyGen(&params.Parameters, mkrlwe.CopyNewDecomposed(crs))

	encryptor := mkbfv.NewMKEncryptor(keys.PublicKey, params)
	decryptor := mkrlwe.NewMKDecryptor(&params.Parameters)
	encoder := bfv.NewEncoder(*params)
	ringQ := mkrlwe.GetRingQ(&params.Parameters)

	return &mkParticipant{
		encryptor:     encryptor,
		decryptor:     decryptor,
		params:        params,
		keys:          keys,
		encoder:       encoder,
		ringQ:         ringQ,
		sigmaSmudging: sigmaSmudging,
	}
}

// newPlaintext initializes a new bfv Plaintext with an encoded slice of uint64
func newPlaintext(value []uint64, encoder bfv.Encoder, params *bfv.Parameters) *bfv.Plaintext {

	plaintext := bfv.NewPlaintext(*params)

	// Encode
	encoder.EncodeUint(value, plaintext)

	return plaintext
}
