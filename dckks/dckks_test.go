package dckks

import (
	"fmt"
	"github.com/ldsec/lattigo/utils"
	"math"
	"sort"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

func testString(opname string, parties uint64, params *ckks.Parameters) string {
	return fmt.Sprintf("%sparties=%d/logN=%d/logQ=%d/levels=%d/a=%d/b=%d",
		opname,
		parties,
		params.LogN,
		params.LogQP,
		params.MaxLevel+1,
		params.Alpha,
		params.Beta)
}

type dckksTestContext struct {
	params       *ckks.Parameters
	dckksContext *dckksContext
	encoder      ckks.Encoder
	evaluator    ckks.Evaluator

	encryptorPk0 ckks.Encryptor
	decryptorSk0 ckks.Decryptor
	decryptorSk1 ckks.Decryptor

	pk0 *ckks.PublicKey
	pk1 *ckks.PublicKey

	sk0 *ckks.SecretKey
	sk1 *ckks.SecretKey

	sk0Shards []*ckks.SecretKey
	sk1Shards []*ckks.SecretKey
}

type dckksTestParameters struct {
	parties    uint64
	verbose    bool
	medianprec float64

	ckksParameters []*ckks.Parameters
}

var err error
var testParams = new(dckksTestParameters)

func init() {

	testParams.parties = 10

	testParams.medianprec = 15
	testParams.verbose = false

	testParams.ckksParameters = ckks.DefaultParams[ckks.PN13QP218 : ckks.PN14QP438+2]
}

func TestDCKKS(t *testing.T) {
	t.Run("PublicKeyGen", testPublicKeyGen)
	t.Run("RelinKeyGen", testRelinKeyGen)
	t.Run("RelinKeyGenNaive", testRelinKeyGenNaive)
	t.Run("KeySwitching", testKeyswitching)
	t.Run("PublicKeySwitching", testPublicKeySwitching)
	t.Run("RotKeyGenConjugate", testRotKeyGenConjugate)
	t.Run("RotKeyGenCols", testRotKeyGenCols)
	t.Run("Refresh", testRefresh)
	t.Run("RefreshAndPermute", testRefreshAndPermute)
}

func gendckksTestContext(contextParameters *ckks.Parameters) (params *dckksTestContext) {

	params = new(dckksTestContext)

	params.params = contextParameters.Copy()

	dckksContext := newDckksContext(contextParameters)

	params.dckksContext = dckksContext

	params.encoder = ckks.NewEncoder(contextParameters)
	params.evaluator = ckks.NewEvaluator(contextParameters)

	kgen := ckks.NewKeyGenerator(contextParameters)

	// SecretKeys
	params.sk0Shards = make([]*ckks.SecretKey, testParams.parties)
	params.sk1Shards = make([]*ckks.SecretKey, testParams.parties)
	tmp0 := params.dckksContext.contextQP.NewPoly()
	tmp1 := params.dckksContext.contextQP.NewPoly()

	for j := uint64(0); j < testParams.parties; j++ {
		params.sk0Shards[j] = kgen.GenSecretKey()
		params.sk1Shards[j] = kgen.GenSecretKey()
		params.dckksContext.contextQP.Add(tmp0, params.sk0Shards[j].Get(), tmp0)
		params.dckksContext.contextQP.Add(tmp1, params.sk1Shards[j].Get(), tmp1)
	}

	params.sk0 = new(ckks.SecretKey)
	params.sk1 = new(ckks.SecretKey)

	params.sk0.Set(tmp0)
	params.sk1.Set(tmp1)

	// Publickeys
	params.pk0 = kgen.GenPublicKey(params.sk0)
	params.pk1 = kgen.GenPublicKey(params.sk1)

	params.encryptorPk0 = ckks.NewEncryptorFromPk(contextParameters, params.pk0)

	params.decryptorSk0 = ckks.NewDecryptor(contextParameters, params.sk0)

	params.decryptorSk1 = ckks.NewDecryptor(contextParameters, params.sk1)

	return
}

func testPublicKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQP)

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			crp := crpGenerator.ReadNew()

			type Party struct {
				*CKGProtocol
				s  *ring.Poly
				s1 CKGShare
			}

			ckgParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.CKGProtocol = NewCKGProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.s1 = p.AllocateShares()
				ckgParties[i] = p
			}
			P0 := ckgParties[0]

			// Each party creates a new CKGProtocol instance
			for i, p := range ckgParties {
				p.GenShare(p.s, crp, p.s1)
				if i > 0 {
					P0.AggregateShares(p.s1, P0.s1, P0.s1)
				}
			}

			pk := &ckks.PublicKey{}
			P0.GenPublicKey(P0.s1, crp, pk)

			// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
			encryptorTest := ckks.NewEncryptorFromPk(parameters, pk)

			coeffs, _, ciphertext := newTestVectors(params, encryptorTest, 1, t)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertext, t)
		})
	}
}

func testRelinKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RKGProtocol
				u      *ring.Poly
				s      *ring.Poly
				share1 RKGShareRoundOne
				share2 RKGShareRoundTwo
				share3 RKGShareRoundThree
			}

			rkgParties := make([]*Party, parties)

			for i := range rkgParties {
				p := new(Party)
				p.RKGProtocol = NewEkgProtocol(parameters)
				p.u = p.NewEphemeralKey()
				p.s = sk0Shards[i].Get()
				p.share1, p.share2, p.share3 = p.AllocateShares()
				rkgParties[i] = p
			}

			P0 := rkgParties[0]
			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQP)
			crp := make([]*ring.Poly, parameters.Beta)

			for i := uint64(0); i < parameters.Beta; i++ {
				crp[i] = crpGenerator.ReadNew()
			}

			// ROUND 1
			for i, p := range rkgParties {
				p.GenShareRoundOne(p.u, p.s, crp, p.share1)
				if i > 0 {
					P0.AggregateShareRoundOne(p.share1, P0.share1, P0.share1)
				}
			}

			//ROUND 2
			for i, p := range rkgParties {
				p.GenShareRoundTwo(P0.share1, p.s, crp, p.share2)
				if i > 0 {
					P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
				}
			}

			// ROUND 3
			for i, p := range rkgParties {
				p.GenShareRoundThree(P0.share2, p.u, p.s, p.share3)
				if i > 0 {
					P0.AggregateShareRoundThree(p.share3, P0.share3, P0.share3)
				}
			}

			evk := ckks.NewRelinKey(parameters)
			P0.GenRelinearizationKey(P0.share2, P0.share3, evk)

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, 1, t)

			for i := range coeffs {
				coeffs[i] *= coeffs[i]
			}

			evaluator.MulRelin(ciphertext, ciphertext, evk, ciphertext)

			evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)

			require.Equal(t, ciphertext.Degree(), uint64(1))

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertext, t)

		})
	}
}

func testRelinKeyGenNaive(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		evaluator := params.evaluator
		pk0 := params.pk0
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RKGProtocolNaive
				u      *ring.Poly
				s      *ring.Poly
				share1 RKGNaiveShareRoundOne
				share2 RKGNaiveShareRoundTwo
			}

			rkgParties := make([]*Party, parties)

			for i := range rkgParties {
				p := new(Party)
				p.RKGProtocolNaive = NewRKGProtocolNaive(parameters)
				p.s = sk0Shards[i].Get()
				p.share1, p.share2 = p.AllocateShares()
				rkgParties[i] = p
			}

			P0 := rkgParties[0]

			// ROUND 1
			for i, p := range rkgParties {
				rkgParties[i].GenShareRoundOne(p.s, pk0.Get(), p.share1)
				if i > 0 {
					P0.AggregateShareRoundOne(p.share1, P0.share1, P0.share1)
				}
			}

			// ROUND 2
			for i, p := range rkgParties {
				rkgParties[i].GenShareRoundTwo(P0.share1, p.s, pk0.Get(), p.share2)
				if i > 0 {
					P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
				}
			}

			evk := ckks.NewRelinKey(parameters)
			P0.GenRelinearizationKey(P0.share2, evk)

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, 1, t)

			for i := range coeffs {
				coeffs[i] *= coeffs[i]
			}

			evaluator.MulRelin(ciphertext, ciphertext, evk, ciphertext)

			require.Equal(t, ciphertext.Degree(), uint64(1))

			evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertext, t)
		})
	}
}

func testKeyswitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		encryptorPk0 := params.encryptorPk0
		decryptorSk1 := params.decryptorSk1
		sk0Shards := params.sk0Shards
		sk1Shards := params.sk1Shards

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*CKSProtocol
				s0    *ring.Poly
				s1    *ring.Poly
				share CKSShare
			}

			cksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.CKSProtocol = NewCKSProtocol(parameters, 6.36)
				p.s0 = sk0Shards[i].Get()
				p.s1 = sk1Shards[i].Get()
				p.share = p.AllocateShare()
				cksParties[i] = p
			}
			P0 := cksParties[0]

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, 1, t)

			// Each party creates its CKSProtocol instance with tmp = si-si'
			for i, p := range cksParties {
				p.GenShare(p.s0, p.s1, ciphertext, p.share)
				if i > 0 {
					P0.AggregateShares(p.share, P0.share, P0.share)
				}
			}

			ksCiphertext := ckks.NewCiphertext(parameters, 1, ciphertext.Level(), ciphertext.Scale())

			P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

			verifyTestVectors(params, decryptorSk1, coeffs, ksCiphertext, t)

			P0.KeySwitch(P0.share, ciphertext, ciphertext)

			verifyTestVectors(params, decryptorSk1, coeffs, ksCiphertext, t)

		})
	}
}

func testPublicKeySwitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		encryptorPk0 := params.encryptorPk0
		decryptorSk1 := params.decryptorSk1
		sk0Shards := params.sk0Shards
		pk1 := params.pk1

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, 1, t)

			params.evaluator.DropLevel(ciphertext, 1)

			type Party struct {
				*PCKSProtocol
				s     *ring.Poly
				share PCKSShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.PCKSProtocol = NewPCKSProtocol(parameters, 6.36)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShares(ciphertext.Level())
				pcksParties[i] = p
			}
			P0 := pcksParties[0]

			ciphertextSwitched := ckks.NewCiphertext(parameters, 1, ciphertext.Level(), ciphertext.Scale())

			for i, p := range pcksParties {
				p.GenShare(p.s, pk1, ciphertext, p.share)
				if i > 0 {
					P0.AggregateShares(p.share, P0.share, P0.share)
				}
			}

			P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

			verifyTestVectors(params, decryptorSk1, coeffs, ciphertextSwitched, t)
		})
	}
}

func testRotKeyGenConjugate(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		contextKeys := params.dckksContext.contextQP
		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RTGProtocol
				s     *ring.Poly
				share RTGShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RTGProtocol = NewRotKGProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShare()
				pcksParties[i] = p
			}
			P0 := pcksParties[0]
			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQP)
			crp := make([]*ring.Poly, parameters.Beta)

			for i := uint64(0); i < parameters.Beta; i++ {
				crp[i] = crpGenerator.ReadNew()
			}

			for i, p := range pcksParties {
				p.GenShare(ckks.Conjugate, 0, p.s, crp, &p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			rotkey := ckks.NewRotationKeys()
			P0.Finalize(parameters, P0.share, crp, rotkey)

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, 1, t)

			evaluator.Conjugate(ciphertext, rotkey, ciphertext)

			coeffsWant := make([]complex128, contextKeys.N>>1)

			for i := uint64(0); i < contextKeys.N>>1; i++ {
				coeffsWant[i] = complex(real(coeffs[i]), -imag(coeffs[i]))
			}

			verifyTestVectors(params, decryptorSk0, coeffsWant, ciphertext, t)

		})
	}
}

func testRotKeyGenCols(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		contextKeys := params.dckksContext.contextQP
		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RTGProtocol
				s     *ring.Poly
				share RTGShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RTGProtocol = NewRotKGProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShare()
				pcksParties[i] = p
			}

			P0 := pcksParties[0]
			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, contextKeys)
			crp := make([]*ring.Poly, parameters.Beta)

			for i := uint64(0); i < parameters.Beta; i++ {
				crp[i] = crpGenerator.ReadNew()
			}

			mask := (contextKeys.N >> 1) - 1

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, 1, t)

			receiver := ckks.NewCiphertext(parameters, ciphertext.Degree(), ciphertext.Level(), ciphertext.Scale())

			for k := uint64(1); k < contextKeys.N>>1; k <<= 1 {

				for i, p := range pcksParties {
					p.GenShare(ckks.RotationLeft, k, p.s, crp, &p.share)
					if i > 0 {
						P0.Aggregate(p.share, P0.share, P0.share)
					}
				}

				rotkey := ckks.NewRotationKeys()
				P0.Finalize(parameters, P0.share, crp, rotkey)

				evaluator.RotateColumns(ciphertext, k, rotkey, receiver)

				coeffsWant := make([]complex128, contextKeys.N>>1)

				for i := uint64(0); i < contextKeys.N>>1; i++ {
					coeffsWant[i] = coeffs[(i+k)&mask]
				}

				verifyTestVectors(params, decryptorSk0, coeffsWant, receiver, t)
			}
		})
	}
}

func testRefresh(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		levelStart := uint64(3)

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RefreshProtocol
				s      *ring.Poly
				share1 RefreshShareDecrypt
				share2 RefreshShareRecrypt
			}

			RefreshParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RefreshProtocol = NewRefreshProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share1, p.share2 = p.AllocateShares(levelStart)
				RefreshParties[i] = p
			}

			P0 := RefreshParties[0]
			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQ)
			crp := crpGenerator.ReadNew()

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, 1.0, t)

			for ciphertext.Level() != levelStart {
				evaluator.DropLevel(ciphertext, 1)
			}

			for i, p := range RefreshParties {
				p.GenShares(p.s, levelStart, parties, ciphertext, crp, p.share1, p.share2)
				if i > 0 {
					P0.Aggregate(p.share1, P0.share1, P0.share1)
					P0.Aggregate(p.share2, P0.share2, P0.share2)
				}
			}

			// We refresh the ciphertext with the simulated error
			P0.Decrypt(ciphertext, P0.share1)      // Masked decryption
			P0.Recode(ciphertext)                  // Masked re-encoding
			P0.Recrypt(ciphertext, crp, P0.share2) // Masked re-encryption

			require.Equal(t, ciphertext.Level(), parameters.MaxLevel)

			verifyTestVectors(params, decryptorSk0, coeffs, ciphertext, t)

		})
	}
}

func testRefreshAndPermute(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.ckksParameters {

		params := gendckksTestContext(parameters)

		evaluator := params.evaluator
		encryptorPk0 := params.encryptorPk0
		decryptorSk0 := params.decryptorSk0
		sk0Shards := params.sk0Shards

		levelStart := uint64(3)

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*PermuteProtocol
				s      *ring.Poly
				share1 RefreshShareDecrypt
				share2 RefreshShareRecrypt
			}

			RefreshParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.PermuteProtocol = NewPermuteProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share1, p.share2 = p.AllocateShares(levelStart)
				RefreshParties[i] = p
			}

			P0 := RefreshParties[0]

			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQ)
			crp := crpGenerator.ReadNew()

			coeffs, _, ciphertext := newTestVectors(params, encryptorPk0, 1.0, t)

			for ciphertext.Level() != levelStart {
				evaluator.DropLevel(ciphertext, 1)
			}

			permutation := make([]uint64, parameters.Slots)

			for i := range permutation {
				permutation[i] = ring.RandUniform(prng, parameters.Slots, parameters.Slots-1)
			}

			for i, p := range RefreshParties {
				p.GenShares(p.s, levelStart, parties, ciphertext, crp, parameters.Slots, permutation, p.share1, p.share2)
				if i > 0 {
					P0.Aggregate(p.share1, P0.share1, P0.share1)
					P0.Aggregate(p.share2, P0.share2, P0.share2)
				}
			}

			// We refresh the ciphertext with the simulated error
			P0.Decrypt(ciphertext, P0.share1)                     // Masked decryption
			P0.Permute(ciphertext, permutation, parameters.Slots) // Masked re-encoding
			P0.Recrypt(ciphertext, crp, P0.share2)                // Masked re-encryption

			coeffsPermute := make([]complex128, len(coeffs))

			for i := range coeffs {
				coeffsPermute[i] = coeffs[permutation[i]]
			}

			require.Equal(t, ciphertext.Level(), parameters.MaxLevel)

			verifyTestVectors(params, decryptorSk0, coeffsPermute, ciphertext, t)

		})
	}
}

func newTestVectors(contextParams *dckksTestContext, encryptor ckks.Encryptor, a float64, t *testing.T) (values []complex128, plaintext *ckks.Plaintext, ciphertext *ckks.Ciphertext) {

	slots := uint64(1 << contextParams.params.LogSlots)

	values = make([]complex128, slots)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	for i := uint64(0); i < slots; i++ {
		values[i] = randomComplex(prng, a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = ckks.NewPlaintext(contextParams.params, contextParams.params.MaxLevel, contextParams.params.Scale)

	contextParams.encoder.Encode(plaintext, values, slots)

	ciphertext = encryptor.EncryptNew(plaintext)

	return values, plaintext, ciphertext
}

func verifyTestVectors(contextParams *dckksTestContext, decryptor ckks.Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	var plaintextTest *ckks.Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *ckks.Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*ckks.Ciphertext))
	case *ckks.Plaintext:
		plaintextTest = element.(*ckks.Plaintext)
	}

	slots := uint64(1 << contextParams.params.LogSlots)

	valuesTest = contextParams.encoder.Decode(plaintextTest, slots)

	var deltaReal, deltaImag float64

	var minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, slots)

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	distribReal := make(map[uint64]uint64)
	distribImag := make(map[uint64]uint64)

	for i := range valuesWant {

		deltaReal = math.Abs(real(valuesTest[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesTest[i]) - imag(valuesWant[i]))

		diff[i] += complex(deltaReal, 0)
		diff[i] += complex(0, deltaImag)

		meanprec += diff[i]

		if real(diff[i]) > real(minprec) {
			minprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) > imag(minprec) {
			minprec = complex(real(minprec), imag(diff[i]))
		}

		if real(diff[i]) < real(maxprec) {
			maxprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) < imag(maxprec) {
			maxprec = complex(real(maxprec), imag(diff[i]))
		}

		distribReal[uint64(math.Floor(math.Log2(1/real(diff[i]))))]++
		distribImag[uint64(math.Floor(math.Log2(1/imag(diff[i]))))]++
	}

	meanprec /= complex(float64(slots), 0)
	medianprec = calcmedian(diff)

	if testParams.verbose {
		t.Logf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
		t.Logf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
		t.Logf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
		t.Logf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
		t.Log()
	}

	require.GreaterOrEqual(t, math.Log2(1/real(medianprec)), testParams.medianprec)
	require.GreaterOrEqual(t, math.Log2(1/imag(medianprec)), testParams.medianprec)
}

func calcmedian(values []complex128) (median complex128) {

	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = real(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(tmp[i], imag(values[i]))
	}

	for i := range values {
		tmp[i] = imag(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(real(values[i]), tmp[i])
	}

	index := len(values) / 2

	if len(values)&1 == 1 {
		return values[index]
	}

	if index+1 == len(values) {
		return values[index]
	}

	return (values[index] + values[index+1]) / 2
}
