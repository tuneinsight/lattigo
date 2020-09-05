package dckks

import (
	"flag"
	"fmt"
	"math"
	"sort"
	"testing"

	"github.com/ldsec/lattigo/utils"

	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

var err error
var params = new(testParams)
var defaultParams = ckks.DefaultParams[ckks.PN13QP218 : ckks.PN14QP438+1]

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func testString(opname string, parties uint64, params *ckks.Parameters) string {
	return fmt.Sprintf("%sparties=%d/logN=%d/logQ=%d/levels=%d/a=%d/b=%d",
		opname,
		parties,
		params.LogN(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.Alpha(),
		params.Beta())
}

type testParams struct {
	params     *ckks.Parameters
	parties    uint64
	medianprec float64

	dckksContext *dckksContext

	prng utils.PRNG

	encoder   ckks.Encoder
	evaluator ckks.Evaluator

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

func TestDCKKS(t *testing.T) {

	params.medianprec = 15
	params.parties = 3

	for _, p := range defaultParams {

		if err = genTestParams(p); err != nil {
			panic(err)
		}

		t.Run("PublicKeyGen", testPublicKeyGen)
		t.Run("RelinKeyGen", testRelinKeyGen)
		t.Run("RelinKeyGenNaive", testRelinKeyGenNaive)
		t.Run("KeySwitching", testKeyswitching)
		t.Run("PublicKeySwitching", testPublicKeySwitching)
		t.Run("RotKeyGenConjugate", testRotKeyGenConjugate)
		t.Run("RotKeyGenCols", testRotKeyGenCols)
		t.Run("Refresh", testRefresh)
	}
}

func genTestParams(defaultParams *ckks.DefaultParam) (err error) {

	if params.params, err = ckks.NewParametersFromLogModuli(defaultParams.LogN, defaultParams.LogModuli); err != nil {
		return err
	}

	params.params.SetLogSlots(defaultParams.LogN - 1)
	params.params.SetScale(defaultParams.Scale)

	dckksContext := newDckksContext(params.params)

	params.dckksContext = dckksContext

	params.prng, err = utils.NewPRNG()
	if err != nil {
		return err
	}

	params.encoder = ckks.NewEncoder(params.params)
	params.evaluator = ckks.NewEvaluator(params.params)

	kgen := ckks.NewKeyGenerator(params.params)

	// SecretKeys
	params.sk0Shards = make([]*ckks.SecretKey, params.parties)
	params.sk1Shards = make([]*ckks.SecretKey, params.parties)
	tmp0 := params.dckksContext.contextQP.NewPoly()
	tmp1 := params.dckksContext.contextQP.NewPoly()

	for j := uint64(0); j < params.parties; j++ {
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

	params.encryptorPk0 = ckks.NewEncryptorFromPk(params.params, params.pk0)

	params.decryptorSk0 = ckks.NewDecryptor(params.params, params.sk0)

	params.decryptorSk1 = ckks.NewDecryptor(params.params, params.sk1)

	return nil
}

func testPublicKeyGen(t *testing.T) {

	parties := params.parties
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards
	prng, err := utils.NewKeyedPRNG(nil)
	if err != nil {
		panic(err)
	}

	crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQP)

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		crp := crpGenerator.ReadNew()

		type Party struct {
			*CKGProtocol
			s  *ring.Poly
			s1 CKGShare
		}

		ckgParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.CKGProtocol = NewCKGProtocol(params.params)
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
		encryptorTest := ckks.NewEncryptorFromPk(params.params, pk)

		coeffs, _, ciphertext := newTestVectors(encryptorTest, 1, t)

		verifyTestVectors(decryptorSk0, coeffs, ciphertext, t)
	})

}

func testRelinKeyGen(t *testing.T) {

	parties := params.parties
	evaluator := params.evaluator
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*RKGProtocol
			u      *ring.Poly
			s      *ring.Poly
			share1 RKGShare
			share2 RKGShare
		}

		rkgParties := make([]*Party, parties)

		for i := range rkgParties {
			p := new(Party)
			p.RKGProtocol = NewEkgProtocol(params.params)
			p.u = p.NewEphemeralKey()
			p.s = sk0Shards[i].Get()
			p.share1, p.share2 = p.AllocateShares()
			rkgParties[i] = p
		}

		P0 := rkgParties[0]
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		crpGenerator := ring.NewUniformSampler(prng, params.dckksContext.contextQP)
		crp := make([]*ring.Poly, params.params.Beta())

		for i := uint64(0); i < params.params.Beta(); i++ {
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
			p.GenShareRoundTwo(P0.share1, p.u, p.s, crp, p.share2)
			if i > 0 {
				P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
			}
		}

		evk := ckks.NewRelinKey(params.params)
		P0.GenRelinearizationKey(P0.share1, P0.share2, evk)

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, 1, t)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
		}

		evaluator.MulRelin(ciphertext, ciphertext, evk, ciphertext)

		evaluator.Rescale(ciphertext, params.params.Scale(), ciphertext)

		require.Equal(t, ciphertext.Degree(), uint64(1))

		verifyTestVectors(decryptorSk0, coeffs, ciphertext, t)

	})

}

func testRelinKeyGenNaive(t *testing.T) {

	parties := params.parties
	evaluator := params.evaluator
	pk0 := params.pk0
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards

	t.Run(testString("", parties, params.params), func(t *testing.T) {

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
			p.RKGProtocolNaive = NewRKGProtocolNaive(params.params)
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

		evk := ckks.NewRelinKey(params.params)
		P0.GenRelinearizationKey(P0.share2, evk)

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, 1, t)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
		}

		evaluator.MulRelin(ciphertext, ciphertext, evk, ciphertext)

		require.Equal(t, ciphertext.Degree(), uint64(1))

		evaluator.Rescale(ciphertext, params.params.Scale(), ciphertext)

		verifyTestVectors(decryptorSk0, coeffs, ciphertext, t)
	})
}

func testKeyswitching(t *testing.T) {

	parties := params.parties
	encryptorPk0 := params.encryptorPk0
	decryptorSk1 := params.decryptorSk1
	sk0Shards := params.sk0Shards
	sk1Shards := params.sk1Shards

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*CKSProtocol
			s0    *ring.Poly
			s1    *ring.Poly
			share CKSShare
		}

		cksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.CKSProtocol = NewCKSProtocol(params.params, 6.36)
			p.s0 = sk0Shards[i].Get()
			p.s1 = sk1Shards[i].Get()
			p.share = p.AllocateShare()
			cksParties[i] = p
		}
		P0 := cksParties[0]

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, 1, t)

		// Each party creates its CKSProtocol instance with tmp = si-si'
		for i, p := range cksParties {
			p.GenShare(p.s0, p.s1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		ksCiphertext := ckks.NewCiphertext(params.params, 1, ciphertext.Level(), ciphertext.Scale())

		P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

		verifyTestVectors(decryptorSk1, coeffs, ksCiphertext, t)

		P0.KeySwitch(P0.share, ciphertext, ciphertext)

		verifyTestVectors(decryptorSk1, coeffs, ksCiphertext, t)

	})
}

func testPublicKeySwitching(t *testing.T) {

	parties := params.parties
	encryptorPk0 := params.encryptorPk0
	decryptorSk1 := params.decryptorSk1
	sk0Shards := params.sk0Shards
	pk1 := params.pk1

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, 1, t)

		params.evaluator.DropLevel(ciphertext, 1)

		type Party struct {
			*PCKSProtocol
			s     *ring.Poly
			share PCKSShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.PCKSProtocol = NewPCKSProtocol(params.params, 6.36)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShares(ciphertext.Level())
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		ciphertextSwitched := ckks.NewCiphertext(params.params, 1, ciphertext.Level(), ciphertext.Scale())

		for i, p := range pcksParties {
			p.GenShare(p.s, pk1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

		verifyTestVectors(decryptorSk1, coeffs, ciphertextSwitched, t)
	})
}

func testRotKeyGenConjugate(t *testing.T) {

	parties := params.parties
	contextKeys := params.dckksContext.contextQP
	evaluator := params.evaluator
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *ring.Poly
			share RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(params.params)
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
		crp := make([]*ring.Poly, params.params.Beta())

		for i := uint64(0); i < params.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		for i, p := range pcksParties {
			p.GenShare(ckks.Conjugate, 0, p.s, crp, &p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		rotkey := ckks.NewRotationKeys()
		P0.Finalize(params.params, P0.share, crp, rotkey)

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, 1, t)

		evaluator.Conjugate(ciphertext, rotkey, ciphertext)

		coeffsWant := make([]complex128, contextKeys.N>>1)

		for i := uint64(0); i < contextKeys.N>>1; i++ {
			coeffsWant[i] = complex(real(coeffs[i]), -imag(coeffs[i]))
		}

		verifyTestVectors(decryptorSk0, coeffsWant, ciphertext, t)

	})
}

func testRotKeyGenCols(t *testing.T) {

	parties := params.parties
	contextKeys := params.dckksContext.contextQP
	evaluator := params.evaluator
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *ring.Poly
			share RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(params.params)
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
		crp := make([]*ring.Poly, params.params.Beta())

		for i := uint64(0); i < params.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		mask := (contextKeys.N >> 1) - 1

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, 1, t)

		receiver := ckks.NewCiphertext(params.params, ciphertext.Degree(), ciphertext.Level(), ciphertext.Scale())

		for k := uint64(1); k < contextKeys.N>>1; k <<= 1 {

			for i, p := range pcksParties {
				p.GenShare(ckks.RotationLeft, k, p.s, crp, &p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			rotkey := ckks.NewRotationKeys()
			P0.Finalize(params.params, P0.share, crp, rotkey)

			evaluator.RotateColumns(ciphertext, k, rotkey, receiver)

			coeffsWant := make([]complex128, contextKeys.N>>1)

			for i := uint64(0); i < contextKeys.N>>1; i++ {
				coeffsWant[i] = coeffs[(i+k)&mask]
			}

			verifyTestVectors(decryptorSk0, coeffsWant, receiver, t)
		}
	})
}

func testRefresh(t *testing.T) {

	parties := params.parties
	evaluator := params.evaluator
	encryptorPk0 := params.encryptorPk0
	decryptorSk0 := params.decryptorSk0
	sk0Shards := params.sk0Shards

	levelStart := uint64(3)

	t.Run(testString("", parties, params.params), func(t *testing.T) {

		type Party struct {
			*RefreshProtocol
			s      *ring.Poly
			share1 RefreshShareDecrypt
			share2 RefreshShareRecrypt
		}

		RefreshParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RefreshProtocol = NewRefreshProtocol(params.params)
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

		coeffs, _, ciphertext := newTestVectors(encryptorPk0, 1.0, t)

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

		require.Equal(t, ciphertext.Level(), params.params.MaxLevel())

		verifyTestVectors(decryptorSk0, coeffs, ciphertext, t)

	})
}

func newTestVectors(encryptor ckks.Encryptor, a float64, t *testing.T) (values []complex128, plaintext *ckks.Plaintext, ciphertext *ckks.Ciphertext) {

	slots := params.params.Slots()

	values = make([]complex128, slots)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	for i := uint64(0); i < slots; i++ {
		values[i] = randomComplex(prng, a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = ckks.NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

	params.encoder.Encode(plaintext, values, slots)

	ciphertext = encryptor.EncryptNew(plaintext)

	return values, plaintext, ciphertext
}

func verifyTestVectors(decryptor ckks.Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	var plaintextTest *ckks.Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *ckks.Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*ckks.Ciphertext))
	case *ckks.Plaintext:
		plaintextTest = element.(*ckks.Plaintext)
	}

	slots := params.params.Slots()

	valuesTest = params.encoder.Decode(plaintextTest, slots)

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

	if *printPrecisionStats {
		t.Logf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
		t.Logf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
		t.Logf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
		t.Logf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
		t.Log()
	}

	require.GreaterOrEqual(t, math.Log2(1/real(medianprec)), params.medianprec)
	require.GreaterOrEqual(t, math.Log2(1/imag(medianprec)), params.medianprec)
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
