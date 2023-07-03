package main

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var (
	NbParties = 3
	HSkDense  = 16384
	HSkSparse = 32
)

type Party struct {
	SkDenseLSSS  *rlwe.SecretKey
	SkSparseLSSS *rlwe.SecretKey
}

func main() {

	// First we define the residual CKKS parameters. This is only a template that will be given
	// to the constructor along with the specificities of the bootstrapping circuit we choose, to
	// enable it to create the appropriate ckks.ParametersLiteral that enable the evaluation of the
	// bootstrapping circuit on top of the residual moduli that we defined.
	ckksParamsResidualLit := ckks.ParametersLiteral{
		LogN:              15,            // Log2 of the ringdegree
		LogQ:              []int{59, 48}, // Log2 of the ciphertext prime moduli
		LogP:              []int{60, 60}, // Log2 of the key-switch auxiliary prime moduli
		LogPlaintextScale: 48,            // Log2 of the scale
		Xs:                ring.Ternary{H: HSkDense},
	}

	btpParametersLit := bootstrapping.ParametersLiteral{
		SlotsToCoeffsFactorizationDepthAndLogPlaintextScales: [][]int{{30, 30}},
		CoeffsToSlotsFactorizationDepthAndLogPlaintextScales: [][]int{{52}, {52}},
		EvalModLogPlaintextScale:                             utils.PointyInt(59),
		LogMessageRatio:                                      utils.PointyInt(11),
		IterationsParameters:                                 &bootstrapping.IterationsParameters{BootstrappingPrecision: []float64{16}, ReservedPrimeBitSize: 16},
		DoubleAngle:                                          utils.PointyInt(3),
		EphemeralSecretWeight:                                utils.PointyInt(HSkSparse),
	}

	ckksParamsLitBoot, bootParams, err := bootstrapping.NewParametersFromLiteral(ckksParamsResidualLit, btpParametersLit)
	if err != nil {
		panic(err)
	}

	// This generate ckks.Parameters, with the NTT tables and other pre-computations from the ckks.ParametersLiteral (which is only a template).
	ckksParamsBoot, err := ckks.NewParametersFromLiteral(ckksParamsLitBoot)
	if err != nil {
		panic(err)
	}

	fmt.Printf("CKKS parameters: logN=%d, logSlots=%d, H(%d; %d), logQP=%f, levels=%d, scale=2^%f\n", ckksParamsBoot.LogN(), ckksParamsBoot.PlaintextLogSlots(), ckksParamsBoot.XsHammingWeight(), bootParams.EphemeralSecretWeight, ckksParamsBoot.LogQP(), ckksParamsBoot.QCount(), math.Log2(ckksParamsBoot.PlaintextScale().Float64()))

	skDenseLSSS := GenerateSecretSharedSecretKey(ckksParamsBoot, HSkDense, NbParties)
	skSparseLSSS := GenerateSecretSharedSecretKey(ckksParamsBoot, HSkSparse, NbParties)

	P := make([]Party, NbParties)
	for i := range P {
		P[i].SkDenseLSSS = skDenseLSSS[i]
		P[i].SkSparseLSSS = skSparseLSSS[i]
	}

	fmt.Printf("Generating Collective Public Key\n")
	pk := GenCollectivePublickey(ckksParamsBoot, P)

	galEls := bootParams.GaloisElements(ckksParamsBoot)
	fmt.Printf("Generating Collective Galois Keys (%d)\n", len(galEls))
	gks := GenCollectiveGaloisKeys(ckksParamsBoot, P, galEls)

	fmt.Printf("Generating Collective RelinearizationKey\n")
	rlk := GenCollectiveRelinearizationKey(ckksParamsBoot, P)

	fmt.Printf("Generating Collective EncapsulationKeys\n")
	evk := bootstrapping.EvaluationKeySet{
		MemEvaluationKeySet: rlwe.NewMemEvaluationKeySet(rlk, gks...),
		EvkDtS:              GenCollectiveEvaluationKey(ckksParamsBoot, skDenseLSSS, skSparseLSSS, rlwe.EvaluationKeyParameters{LevelQ: 0, LevelP: 0, BaseTwoDecomposition: 0}),
		EvkStD:              GenCollectiveEvaluationKey(ckksParamsBoot, skSparseLSSS, skDenseLSSS),
	}

	fmt.Printf("Generating Bootstrapper\n")
	var btp *bootstrapping.Bootstrapper
	if btp, err = bootstrapping.NewBootstrapper(ckksParamsBoot, bootParams, &evk); err != nil {
		panic(err)
	}

	encoder := ckks.NewEncoder(ckksParamsBoot)
	encryptor := ckks.NewEncryptor(ckksParamsBoot, pk)

	// Generate a random plaintext with values uniformely distributed in [-1, 1] for the real and imaginary part.
	valuesWant := make([]complex128, ckksParamsBoot.PlaintextSlots())
	for i := range valuesWant {
		valuesWant[i] = sampling.RandComplex128(-1, 1)
	}

	plaintext := ckks.NewPlaintext(ckksParamsBoot, 0)
	if err := encoder.Encode(valuesWant, plaintext); err != nil {
		panic(err)
	}

	// Encrypt
	ciphertext1 := encryptor.EncryptNew(plaintext)

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	valuesTest1 := printDebug(ckksParamsBoot, valuesWant, CollectiveDecryption(ckksParamsBoot, P, ciphertext1), encoder)

	fmt.Println()
	fmt.Println("Bootstrapping...")
	ciphertext2, err := btp.Bootstrap(ciphertext1)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(ckksParamsBoot, valuesTest1, CollectiveDecryption(ckksParamsBoot, P, ciphertext2), encoder)
}

func CollectiveDecryption(params ckks.Parameters, P []Party, ct *rlwe.Ciphertext) (pt *rlwe.Plaintext) {

	cks := make([]drlwe.KeySwitchProtocol, len(P))

	for i := range cks {
		if i == 0 {
			cks[i] = drlwe.NewKeySwitchProtocol(params.Parameters, params.Xe())
		} else {
			cks[i] = cks[0].ShallowCopy()
		}
	}

	shares := make([]drlwe.KeySwitchShare, len(P))
	for i := range shares {
		shares[i] = cks[i].AllocateShare(ct.Level())
	}

	zero := rlwe.NewSecretKey(params.Parameters)

	for i := range shares {
		cks[i].GenShare(P[i].SkDenseLSSS, zero, ct, &shares[i])
		if i > 0 {
			cks[0].AggregateShares(shares[0], shares[i], &shares[0])
		}
	}

	pt = rlwe.NewPlaintext(params, ct.Level())

	cks[0].KeySwitch(ct, shares[0], pt)

	return
}

func GenCollectivePublickey(params ckks.Parameters, P []Party) (pk *rlwe.PublicKey) {
	ckg := make([]drlwe.PublicKeyGenProtocol, len(P))
	for i := range ckg {
		if i == 0 {
			ckg[i] = drlwe.NewPublicKeyGenProtocol(params.Parameters)
		} else {
			ckg[i] = ckg[0].ShallowCopy()
		}
	}

	shares := make([]drlwe.PublicKeyGenShare, len(P))
	for i := range shares {
		shares[i] = ckg[i].AllocateShare()
	}

	prng, err := sampling.NewPRNG()

	if err != nil {
		panic(err)
	}

	crp := ckg[0].SampleCRP(prng)

	for i := range shares {
		ckg[i].GenShare(P[i].SkDenseLSSS, crp, &shares[i])
	}

	for i := 1; i < len(P); i++ {
		ckg[0].AggregateShares(shares[0], shares[i], &shares[0])
	}

	pk = rlwe.NewPublicKey(params)
	ckg[0].GenPublicKey(shares[0], crp, pk)

	return
}

func GenCollectiveEvaluationKey(params ckks.Parameters, skInLSSS, skOutLSSS []*rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (evk *rlwe.EvaluationKey) {

	var evkParamsCpy rlwe.EvaluationKeyParameters

	if len(evkParams) != 0 {
		evkParamsCpy = evkParams[0]
	} else {
		evkParamsCpy = rlwe.EvaluationKeyParameters{LevelQ: params.MaxLevelQ(), LevelP: params.MaxLevelP(), BaseTwoDecomposition: 0}
	}

	evkg := make([]drlwe.EvaluationKeyGenProtocol, NbParties)
	for i := range evkg {
		if i == 0 {
			evkg[i] = drlwe.NewEvaluationKeyGenProtocol(params.Parameters)
		} else {
			evkg[i] = evkg[0].ShallowCopy()
		}
	}

	shares := make([]drlwe.EvaluationKeyGenShare, NbParties)
	for i := range shares {
		shares[i] = evkg[i].AllocateShare(evkParamsCpy)
	}

	prng, err := sampling.NewPRNG()

	if err != nil {
		panic(err)
	}

	crp := evkg[0].SampleCRP(prng, evkParamsCpy)

	for i := range shares {
		evkg[i].GenShare(skInLSSS[i], skOutLSSS[i], crp, &shares[i])
	}

	for i := 1; i < NbParties; i++ {
		evkg[0].AggregateShares(shares[0], shares[i], &shares[0])
	}

	evk = rlwe.NewEvaluationKey(params, evkParamsCpy)
	evkg[0].GenEvaluationKey(shares[0], crp, evk)

	return
}

func GenCollectiveRelinearizationKey(params ckks.Parameters, P []Party) (rlk *rlwe.RelinearizationKey) {

	rkg := make([]drlwe.RelinKeyGenProtocol, len(P))

	for i := range rkg {
		if i == 0 {
			rkg[i] = drlwe.NewRelinKeyGenProtocol(params.Parameters)
		} else {
			rkg[i] = rkg[0].ShallowCopy()
		}
	}

	ephSk := make([]*rlwe.SecretKey, len(P))
	share1 := make([]drlwe.RelinKeyGenShare, len(P))
	share2 := make([]drlwe.RelinKeyGenShare, len(P))

	for i := range rkg {
		ephSk[i], share1[i], share2[i] = rkg[i].AllocateShare()
	}

	prng, err := sampling.NewPRNG()

	if err != nil {
		panic(err)
	}

	crp := rkg[0].SampleCRP(prng)

	for i := range rkg {
		rkg[i].GenShareRoundOne(P[i].SkDenseLSSS, crp, ephSk[i], &share1[i])
	}

	for i := 1; i < len(P); i++ {
		rkg[0].AggregateShares(share1[0], share1[i], &share1[0])
	}

	for i := range rkg {
		rkg[i].GenShareRoundTwo(ephSk[i], P[i].SkDenseLSSS, share1[0], &share2[i])
	}

	for i := 1; i < len(P); i++ {
		rkg[0].AggregateShares(share2[0], share2[i], &share2[0])
	}

	rlk = rlwe.NewRelinearizationKey(params.Parameters)

	rkg[0].GenRelinearizationKey(share1[0], share2[0], rlk)

	return

}

func GenCollectiveGaloisKeys(params ckks.Parameters, P []Party, galEls []uint64) (GaloisKeys []*rlwe.GaloisKey) {

	gkg := make([]drlwe.GaloisKeyGenProtocol, len(P))
	for i := range gkg {
		if i == 0 {
			gkg[i] = drlwe.NewGaloisKeyGenProtocol(params.Parameters)
		} else {
			gkg[i] = gkg[0].ShallowCopy()
		}
	}

	shares := make([]drlwe.GaloisKeyGenShare, len(P))
	for i := range shares {
		shares[i] = gkg[i].AllocateShare()
	}

	GaloisKeys = make([]*rlwe.GaloisKey, len(galEls))

	prng, err := sampling.NewPRNG()

	if err != nil {
		panic(err)
	}

	fmt.Printf("GaloisKeys: ")
	for i, galEl := range galEls {

		fmt.Printf("%d ", galEl)

		crp := gkg[0].SampleCRP(prng)

		for j := range shares {
			gkg[j].GenShare(P[j].SkDenseLSSS, galEl, crp, &shares[j])
		}

		for j := 1; j < len(P); j++ {
			gkg[0].AggregateShares(shares[0], shares[j], &shares[0])
		}

		galoisKey := rlwe.NewGaloisKey(params)

		gkg[0].GenGaloisKey(shares[0], crp, galoisKey)

		GaloisKeys[i] = galoisKey
	}
	fmt.Println()

	return
}

// GenerateSecretSharedSecretKey generates a secret shared secret key mod QP with the specified Hamming weight.
func GenerateSecretSharedSecretKey(params ckks.Parameters, HammingWeight, NbParties int) (shares []*rlwe.SecretKey) {
	kgen := ckks.NewKeyGenerator(params)

	prng, err := sampling.NewPRNG()

	if err != nil {
		panic(err)
	}

	r := params.RingQP()

	sampler := ringqp.NewUniformSampler(prng, *r)

	shares = make([]*rlwe.SecretKey, NbParties)

	shares[0] = kgen.GenSecretKeyWithHammingWeightNew(HammingWeight)

	for i := 1; i < len(shares); i++ {
		shares[i] = &rlwe.SecretKey{
			Value: sampler.ReadNew(),
		}

		r.Sub(shares[0].Value, shares[i].Value, shares[0].Value)
	}

	return
}

func printDebug(params ckks.Parameters, valuesWant []complex128, pt *rlwe.Plaintext, encoder *ckks.Encoder) (valuesTest []complex128) {

	valuesTest = make([]complex128, pt.PlaintextSlots())

	if err := encoder.Decode(pt, valuesTest); err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", pt.Level(), params.LogQLvl(pt.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(pt.PlaintextScale.Float64()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, nil, false)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
