package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	// "unsafe"
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/dckks"
	"github.com/ldsec/lattigo/ring"
	"time"
)

/*use in sigmoid*/
var degree3 = []float64{-0.5, 0.15012, -0.001593}

type party struct {
	sk         *ckks.SecretKey
	rlkEphemSk *ring.Poly

	ckgShare      dckks.CKGShare
	rkgShareOne   dckks.RKGShareRoundOne
	rkgShareTwo   dckks.RKGShareRoundTwo
	rkgShareThree dckks.RKGShareRoundThree
	rtgShare      dckks.RTGShare
	cksShare      dckks.CKSShare

	input       [][]float64
	Dcnum       uint64
	DfactorDim  uint64
	DsamplesDim uint64
	DfdimBits   uint64
	DsdimBits   uint64
	DbBits      uint64
	Dbatch      uint64
	DsBits      uint64

	encrypted []*ckks.Ciphertext
}

func main() {

	startT := time.Now()
	var N uint64 = 1 // number of users

	/*dckks*/
	/*set parameters*/
	var logmoduli ckks.LogModuli
	logmoduli.LogQi = []uint64{45, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30}
	logmoduli.LogPi = []uint64{45, 45, 45}
	params := ckks.NewParametersFromLogModuli(15, 14, 1<<30, logmoduli, 3.2)

	fmt.Println(params.MaxLevel() + 1)

	crsGen := dckks.NewCRPGenerator(params, []byte{'n', 'a', 'n'})
	crsGen.Seed([]byte{'d', 'c', 'k', 'k', 's'})

	crs := crsGen.ClockNew()                    // for the public-key
	crp := make([]*ring.Poly, params.Beta())    // for the relinearization keys
	crpRot := make([]*ring.Poly, params.Beta()) // for the rotation keys
	for i := uint64(0); i < params.Beta(); i++ {
		crp[i] = crsGen.ClockNew()
	}
	for i := uint64(0); i < params.Beta(); i++ {
		crpRot[i] = crsGen.ClockNew()
	}

	colSk := ckks.NewSecretKey(params)
	ckg := dckks.NewCKGProtocol(params)       // public key generation
	rkg := dckks.NewEkgProtocol(params)       // relineariation key generation
	rtg := dckks.NewRotKGProtocol(params)     // rotation keys generation
	cks := dckks.NewCKSProtocol(params, 3.19) // collective public-key re-encryption

	kgen := ckks.NewKeyGenerator(params)

	contextKeys, _ := ring.NewContextWithParams(1<<params.LogN, append(params.Qi, params.Pi...))

	path := "/home/cyx/data/pcs.txt"
	path1 := "/home/cyx/data/pcs.txt"
	var isfirst bool = false
	var factorDimTest, samplesDimTest uint64 = 0, 0
	var factorDim, samplesDim uint64 = 0, 0
	zDataTest := zDataFromFile(path, &factorDimTest, &samplesDimTest, isfirst)
	normalizeZData(zDataTest, factorDimTest, samplesDimTest)

	P := make([]*party, N, N) //users
	for i := range P {
		var factorDimTrain, samplesDimTrain uint64
		pi := &party{}
		pi.sk = kgen.GenSecretKey()

		pi.rlkEphemSk = contextKeys.SampleTernaryMontgomeryNTTNew(1.0 / 3)

		/*Read the file*/
		if i == 0 {
			pi.input = zDataFromFile(path1, &factorDimTrain, &samplesDimTrain, isfirst)
			normalizeZData(pi.input, factorDimTrain, samplesDimTrain)
		}
		// else if i == 1{
		// 	pi.input = zDataFromFile(path2, &factorDimTrain, &samplesDimTrain, isfirst)
		// 	normalizeZData(pi.input, factorDimTrain, samplesDimTrain)
		// }else if i == 2{
		// 	pi.input = zDataFromFile(path3, &factorDimTrain, &samplesDimTrain, isfirst)
		// 	normalizeZData(pi.input, factorDimTrain, samplesDimTrain)
		// }
		pi.DfactorDim = factorDimTrain
		pi.DsamplesDim = samplesDimTrain
		pi.DfdimBits = uint64(math.Ceil(math.Log2(float64(pi.DfactorDim))))
		pi.DsdimBits = uint64(math.Ceil(math.Log2(float64(pi.DsamplesDim))))
		pi.DbBits = uint64(math.Min(float64(15-1-pi.DsdimBits), float64(pi.DfdimBits)))
		pi.DsBits = pi.DsdimBits + pi.DbBits

		factorDim = factorDimTrain
		samplesDim += samplesDimTrain

		contextKeys.Add(colSk.Get(), pi.sk.Get(), colSk.Get()) //TODO: doc says "return"

		pi.ckgShare = ckg.AllocateShares()
		pi.rkgShareOne, pi.rkgShareTwo, pi.rkgShareThree = rkg.AllocateShares()
		pi.rtgShare = rtg.AllocateShare()
		pi.cksShare = cks.AllocateShare()

		P[i] = pi
	}

	var gammaUp, gammaDown float64
	var numIter uint64
	var kdeg uint64

	fmt.Printf("number of iterations = ")
	fmt.Scanf("%d\n", &numIter)
	fmt.Printf("gammaUp = ")
	fmt.Scanf("%f\n", &gammaUp)
	fmt.Printf("gammaDown = ")
	fmt.Scanf("%f\n", &gammaDown)
	fmt.Printf("g_k = ")
	fmt.Scanf("%d\n", &kdeg)

	/*print paramwters*/
	fmt.Printf("CKKS parameters : logN = %d, logQ = %d, levels = %d, scale= %f, sigma = %f \n",
		params.LogN, params.LogQP(), params.MaxLevel()+1, math.Log2(params.Scale), params.Sigma)

	var fdimBits uint64 = uint64(math.Ceil(math.Log2(float64(factorDim))))
	var sdimBits uint64 = uint64(math.Ceil(math.Log2(float64(samplesDim))))
	var bBits uint64 = uint64(math.Min(float64(15-1-sdimBits), float64(fdimBits)))
	var batch uint64 = 1 << bBits
	var sBits uint64 = sdimBits + bBits
	var cnum uint64 = 0

	/*Determine the cnum of each party*/
	for i := range P {
		P[i].Dcnum = uint64(math.Ceil(float64(P[i].DfactorDim) / float64(batch)))
		// P[i].Dcnum = 1
		cnum += P[i].Dcnum
	}
	slots := uint64(1 << (params.LogN - 1))
	// var slots = params

	fmt.Printf("isfirst = %t, number of iterations = %d, g_k = %d \n", isfirst, numIter, kdeg)
	fmt.Printf("factor = %d, train samples = %d, gammaUp = %f \n", factorDim, samplesDim, gammaUp)
	fmt.Printf("gammaDown = %f, batch = %d, slots = %d, cnum = %d \n", gammaDown, batch, slots, cnum)

	// 1) Collective public key generation
	pk := ckks.NewPublicKey(params)

	for _, pi := range P {
		ckg.GenShare(pi.sk.Get(), crs, pi.ckgShare)
	}
	ckgCombined := ckg.AllocateShares()
	for _, pi := range P {
		ckg.AggregateShares(pi.ckgShare, ckgCombined, ckgCombined)
	}
	ckg.GenPublicKey(ckgCombined, crs, pk)

	// 2) Collective relinearization key generation
	for _, pi := range P {
		rkg.GenShareRoundOne(pi.rlkEphemSk, pi.sk.Get(), crp, pi.rkgShareOne)
	}
	rkgCombined1, rkgCombined2, rkgCombined3 := rkg.AllocateShares()
	for _, pi := range P {
		rkg.AggregateShareRoundOne(pi.rkgShareOne, rkgCombined1, rkgCombined1)
	}
	for _, pi := range P {
		rkg.GenShareRoundTwo(rkgCombined1, pi.sk.Get(), crp, pi.rkgShareTwo)
	}
	for _, pi := range P {
		rkg.AggregateShareRoundTwo(pi.rkgShareTwo, rkgCombined2, rkgCombined2)
	}
	for _, pi := range P {
		rkg.GenShareRoundThree(rkgCombined2, pi.rlkEphemSk, pi.sk.Get(), pi.rkgShareThree)
	}
	rlk := ckks.NewRelinKey(params)
	for _, pi := range P {
		rkg.AggregateShareRoundThree(pi.rkgShareThree, rkgCombined3, rkgCombined3)
	}
	rkg.GenRelinearizationKey(rkgCombined2, rkgCombined3, rlk)

	// 3) Collective rotation keys geneneration
	rtk := ckks.NewRotationKeys()

	for _, rot := range []ckks.Rotation{ckks.RotationRight, ckks.RotationLeft, ckks.Conjugate} {
		for k := uint64(1); (rot == ckks.Conjugate && k == 1) || (rot != ckks.Conjugate && k < 1<<(params.LogN-1)); k <<= 1 {

			rtgShareCombined := rtg.AllocateShare()
			rtgShareCombined.Type = rot
			rtgShareCombined.K = k

			for _, pi := range P {
				rtg.GenShare(rot, k, pi.sk.Get(), crpRot, &pi.rtgShare)
			}
			for _, pi := range P {
				rtg.Aggregate(pi.rtgShare, rtgShareCombined, rtgShareCombined)
			}
			rtg.Finalize(params, rtgShareCombined, crpRot, rtk)
		}
	}

	var pwData []float64 = make([]float64, factorDim, factorDim)
	var pvData []float64 = make([]float64, factorDim, factorDim)

	var twData []float64 = make([]float64, factorDim, factorDim)
	var tvData []float64 = make([]float64, factorDim, factorDim)

	var cwData []float64 = make([]float64, factorDim, factorDim)

	encZData := make([]*ckks.Ciphertext, cnum, cnum)
	encWData := make([]*ckks.Ciphertext, cnum, cnum)
	encVData := make([]*ckks.Ciphertext, cnum, cnum)
	for i := range encWData {
		encWData[i] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale)
	}
	for i := range encVData {
		encVData[i] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale)
	}
	for i := range encZData {
		encZData[i] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale)
	}

	/*generate Encoder and Encryptor*/
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptorFromPk(params, pk)

	var j uint64
	repoly := make([]complex128, slots, slots)
	for j = 0; j < slots; j += batch {
		repoly[j] = complex(1.0, 0)
	}
	prepoly := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale)
	var crepoly *ckks.Ciphertext
	encoder.Encode(prepoly, repoly, slots)
	crepoly = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale)
	encryptor.Encrypt(prepoly, crepoly)

	/*encrypt the data*/
	startenc := time.Now()
	myencZData(P, slots, factorDim, samplesDim, batch, cnum, encoder, encryptor, params)
	elapsedenc := time.Now().Sub(startenc)
	fmt.Println("enc Time", elapsedenc)

	CentralizeEncData(P, encZData, N, cnum)

	/*encWVDataZero*/
	encWVDataZero(encWData, encVData, cnum, slots, encoder, encryptor, params)

	/*init wData and vData*/
	initialWDataVDataZero(pwData, pvData, factorDim)
	initialWDataVDataZero(twData, tvData, factorDim)

	decryptor := ckks.NewDecryptor(params, P[0].sk)
	evaluator := ckks.NewEvaluator(params)
	var eta, gamma float64
	var enccor, encauc float64
	// var  trueauc, truecor float64

	var alpha0 float64 = 0.01
	var alpha1 float64 = (1.0 + math.Sqrt(1.0+4.0*alpha0*alpha0)) / 2.0

	var iter uint64
	startlearn := time.Now()
	for iter = 0; iter < numIter; iter++ {
		fmt.Printf("!!! START %d ITERATION !!!\n", iter+1)
		eta = (1.0 - alpha0) / alpha1
		if gammaDown > 0 {
			gamma = gammaUp / gammaDown / float64(samplesDim)
		} else {
			gamma = gammaUp / (float64(iter) - gammaDown) / float64(samplesDim)
		}

		encGrad := make([]*ckks.Ciphertext, cnum, cnum)
		var encIP *ckks.Ciphertext
		encIP = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale)
		startlearnOne := time.Now()
		encIP = encInnerProduct(encIP, crepoly, encZData, encWData, cnum, bBits, slots, batch, params, evaluator, rlk, rtk)
		encSigmoid(encZData, encGrad, encIP, gamma, cnum, sBits, bBits, encoder, encryptor, params, evaluator, rlk, rtk)
		encNLGDstep(encWData, encVData, encGrad, eta, cnum, iter, encoder, encryptor, params, evaluator, rlk)
		elapsedlearnOne := time.Now().Sub(startlearnOne)
		fmt.Println("learn Time Once", elapsedlearnOne)
		startDec := time.Now()
		decWData(cwData, encWData, batch, slots, factorDim, cnum, decryptor, cks, params, encoder, P)
		elapsedDec := time.Now().Sub(startDec)
		fmt.Println("dec Time Once", elapsedDec)
		caculateAUC(zDataTest, cwData, factorDim, samplesDimTest, enccor, encauc)

		alpha0 = alpha1
		alpha1 = (1.0 + math.Sqrt(1.0+4.0*alpha0*alpha0)) / 2.0
		fmt.Printf("!!! STOP %d ITERATION !!!\n", iter+1)
	}
	elapsedlearn := time.Now().Sub(startlearn)
	fmt.Println("learn Time", elapsedlearn)
	elapsed := time.Now().Sub(startT)
	fmt.Println("ALL Time:", elapsed)
}

func normalizeZData(zData [][]float64, factorDim uint64, samplesDim uint64) {
	var m float64
	for i := 0; i < int(factorDim); i++ {
		m = 0.0
		for j := 0; j < int(samplesDim); j++ {
			m = math.Max(m, math.Abs(zData[j][i]))
		}

		if m < 1e-10 {
			continue
		}

		for j := 0; j < int(samplesDim); j++ {
			zData[j][i] /= m
		}
	}
}

/*Read file from path*/
func zDataFromFile(path string, factorDim *uint64, samplesDim *uint64, isfirst bool) (data [][]float64) {
	zline := [][]float64{}
	*factorDim = uint64(0)  //dimension fo x
	*samplesDim = uint64(0) //number of samples
	start, end := 0, 0
	flag := 0
	/*read the file*/
	file, _ := os.Open(path)
	defer file.Close()
	scanner := bufio.NewScanner(file)
	/*According to the line read. And save it in zData*/
	for scanner.Scan() {
		vecline := []float64{}
		line := scanner.Text()
		strline := string(line)
		if flag == 0 {
			for i := 0; i < len(strline); i++ {
				if strline[i] == ',' {
					*factorDim++
				}
			}
			flag++
		} else {
			for i := 0; i < int(*factorDim); i++ {
				start = 0
				if i != int(*factorDim-1) {
					end = strings.Index(strline, ",")
				}
				float, _ := strconv.ParseFloat(strline[start:end], 64)
				end = strings.Index(strline, ",")
				start = end + 1
				if i != int(*factorDim-1) {
					strline = strline[start:]
				}
				vecline = append(vecline, float)
			}
			*samplesDim++
			zline = append(zline, vecline)
		}
	}

	/*init 2D slice*/
	zData := make([][]float64, *samplesDim)
	for i := 0; i < int(*samplesDim); i++ {
		zData[i] = make([]float64, *factorDim)
	}
	if isfirst {
		for j := 0; j < int(*samplesDim); j++ {
			zj := make([]float64, *factorDim)
			zj[0] = 2*zline[j][0] - 1
			for i := 1; i < int(*factorDim); i++ {
				zj[i] = zj[0] * zline[j][i]
			}
			zData[j] = zj
		}
	} else {
		for j := 0; j < int(*samplesDim); j++ {
			zj := make([]float64, *factorDim)
			zj[0] = 2*zline[j][*factorDim-1] - 1
			for i := 1; i < int(*factorDim); i++ {
				zj[i] = zj[0] * zline[j][i-1]
			}
			zData[j] = zj
		}
	}

	return zData
}

func myencZData(P []*party, slots, factorDim, samplesDim, batch, cnum uint64, encoder ckks.Encoder,
	encryptor ckks.Encryptor, params *ckks.Parameters) {

	for it := range P {
		encZData := make([]*ckks.Ciphertext, P[it].Dcnum, P[it].Dcnum)
		for t := range encZData {
			encZData[t] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale)
		}

		pzData := make([]complex128, slots, slots)
		var i, j, l uint64
		for i = 0; i < (P[it].Dcnum - 1); i++ {
			for j = 0; j < samplesDim; j++ {
				for l = 0; l < batch; l++ {
					pzData[batch*j+l] = complex(P[it].input[j][batch*i+l], 0)
				}
			}
			pt := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale)
			encoder.Encode(pt, pzData, slots)
			encZData[i] = encryptor.EncryptNew(pt)
		}

		var rest uint64 = factorDim - batch*(P[it].Dcnum-1)
		for j = 0; j < P[it].DsamplesDim; j++ {
			for l = 0; l < rest; l++ {
				pzData[batch*j+l] = complex(P[it].input[j][batch*(P[it].Dcnum-1)+l], 0)
			}
			for l = rest; l < batch; l++ {
				pzData[batch*j+l] = complex(0, 0)
			}
		}
		pt := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale)
		encoder.Encode(pt, pzData, slots)
		encZData[P[it].Dcnum-1] = encryptor.EncryptNew(pt)
		P[it].encrypted = encZData
	}
}
func CentralizeEncData(P []*party, encZData []*ckks.Ciphertext, N, cnum uint64) {
	var i, j uint64
	var f uint64 = 0
	for i = 0; i < N; i++ {
		for j = 0; j < P[i].Dcnum; j++ {
			encZData[f] = P[i].encrypted[j]
			f++
		}
	}
}

func initialWDataVDataZero(wData, vData []float64, factorDim uint64) {
	var i uint64
	vData = make([]float64, factorDim, factorDim)
	wData = make([]float64, factorDim, factorDim)
	for i = 0; i < factorDim; i++ {
		wData[i] = 0.0
		vData[i] = 0.0
	}
}

func encWVDataZero(encWData, encVData []*ckks.Ciphertext, cnum, slots uint64,
	encoder ckks.Encoder, encryptor ckks.Encryptor, params *ckks.Parameters) {
	zerow := make([]complex128, slots, slots)
	var i, j uint64
	for j = 0; j < slots; j++ {
		zerow[j] = complex(0, 0)
	}
	for i = 0; i < cnum; i++ {
		pt := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale)
		encoder.Encode(pt, zerow, slots)
		encWData[i] = encryptor.EncryptNew(pt)
		encVData[i] = encWData[i]
	}
}
func encInnerProduct(encIP, crepoly *ckks.Ciphertext, encZData, encWData []*ckks.Ciphertext, cnum, bBits, slots, batch uint64,
	params *ckks.Parameters, evaluator ckks.Evaluator, rlk *ckks.EvaluationKey, rtk *ckks.RotationKeys) (encout *ckks.Ciphertext) {
	var i, l uint64
	encIPvec := make([]*ckks.Ciphertext, cnum, cnum)
	/*for memory*/
	for i := range encIPvec {
		encIPvec[i] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale)
	}
	for i = 0; i < cnum; i++ {
		encIPvec[i] = evaluator.DropLevelNew(encZData[i], (encZData[i].Level() - encWData[i].Level()))
		evaluator.MulRelin(encIPvec[i], encWData[i], rlk, encIPvec[i])
		evaluator.Rescale(encIPvec[i], params.Scale, encIPvec[i])

		for l = 0; l < bBits; l++ {
			rot := evaluator.RotateColumnsNew(encIPvec[i], 1<<l, rtk)
			evaluator.Add(encIPvec[i], rot, encIPvec[i])
		}
	}
	encIP = encIPvec[0]
	for i = 1; i < cnum; i++ {
		evaluator.Add(encIP, encIPvec[i], encIP)
	}

	crepoly = evaluator.DropLevelNew(crepoly, crepoly.Level()-encIP.Level())
	evaluator.MulRelin(encIP, crepoly, rlk, encIP)
	evaluator.Rescale(encIP, params.Scale, encIP)
	for l = 0; l < bBits; l++ {
		tmp := evaluator.RotateColumnsNew(encIP, (slots - (1 << l)), rtk)
		evaluator.Add(encIP, tmp, encIP)
	}
	return encIP
}
func encSigmoid(encZData, encGrad []*ckks.Ciphertext, encIP *ckks.Ciphertext, gamma float64, cnum, sBits, bBits uint64,
	encoder ckks.Encoder, encryptor ckks.Encryptor, params *ckks.Parameters, evaluator ckks.Evaluator, rlk *ckks.EvaluationKey,
	rtk *ckks.RotationKeys) {
	var i, l uint64

	encIP2 := evaluator.MulRelinNew(encIP, encIP, rlk)
	evaluator.Rescale(encIP2, params.Scale, encIP2)

	for i := range encGrad {
		encGrad[i] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.Scale)
	}
	evaluator.AddConst(encIP2, degree3[1]/degree3[2], encIP2)

	for i = 0; i < cnum; i++ {
		evaluator.MultByConst(encZData[i], gamma*degree3[2], encGrad[i])
		evaluator.Rescale(encGrad[i], params.Scale, encGrad[i])

		encGrad[i] = evaluator.DropLevelNew(encGrad[i], encGrad[i].Level()-encIP.Level())
		evaluator.MulRelin(encGrad[i], encIP, rlk, encGrad[i])
		evaluator.Rescale(encGrad[i], params.Scale, encGrad[i])

		evaluator.MulRelin(encGrad[i], encIP2, rlk, encGrad[i])
		evaluator.Rescale(encGrad[i], params.Scale, encGrad[i])

		tmp := evaluator.MultByConstNew(encZData[i], gamma*degree3[0])
		evaluator.Rescale(tmp, params.Scale, tmp)
		tmp = evaluator.DropLevelNew(tmp, tmp.Level()-encGrad[0].Level())
		evaluator.Add(encGrad[i], tmp, encGrad[i])
	}

	for i = 0; i < cnum; i++ {
		for l = bBits; l < sBits; l++ {
			tmp := evaluator.RotateColumnsNew(encGrad[i], (1 << l), rtk)
			evaluator.Add(encGrad[i], tmp, encGrad[i])
		}
	}
}
func encNLGDstep(encWData, encVData, encGrad []*ckks.Ciphertext, eta float64, cnum, iter uint64,
	encoder ckks.Encoder, encryptor ckks.Encryptor, params *ckks.Parameters, evaluator ckks.Evaluator,
	rlk *ckks.EvaluationKey) {
	var i uint64
	for i = 0; i < cnum; i++ {
		encVData[i] = evaluator.DropLevelNew(encVData[i], (encVData[i].Level() - encGrad[i].Level()))
		ctmpw := evaluator.SubNew(encVData[i], encGrad[i])

		if iter != 5 {
			evaluator.MultByConst(ctmpw, 1.0-eta, encVData[i])
			evaluator.Rescale(encVData[i], params.Scale, encVData[i])

			evaluator.MultByConst(encWData[i], eta, encWData[i])
			evaluator.Rescale(encWData[i], params.Scale, encWData[i])

			encWData[i] = evaluator.DropLevelNew(encWData[i], (encWData[i].Level() - encVData[i].Level()))
			evaluator.Add(encVData[i], encWData[i], encVData[i])
		}

		encWData[i] = ctmpw
	}
}
func decWData(wData []float64, encWData []*ckks.Ciphertext, batch, slots, factorDim, cnum uint64, decryptor ckks.Decryptor,
	cks *dckks.CKSProtocol, params *ckks.Parameters, encoder ckks.Encoder, P []*party) {
	var i, j uint64
	for i = 0; i < (cnum - 1); i++ {
		/*CKS*/
		zero := params.NewPolyQ()
		cksCombined := (*cks).AllocateShare()
		for _, pi := range P[1:] {
			cks.GenShare(pi.sk.Get(), zero, encWData[i], pi.cksShare)
		}
		encOut := ckks.NewCiphertext(params, 1, encWData[i].Level(), params.Scale)
		for _, pi := range P {
			cks.AggregateShares(pi.cksShare, cksCombined, cksCombined)
		}
		cks.KeySwitch(cksCombined, encWData[i], encOut)

		ptmp := decryptor.DecryptNew(encOut)
		dcw := encoder.Decode(ptmp, slots)

		for j = 0; j < slots; j++ {
			dcw[j] = round(dcw[j])
		}
		for j = 0; j < batch; j++ {
			wData[batch*i+j] = real(dcw[j])
		}
	}
	/*CKS*/
	zero := params.NewPolyQ()
	cksCombined := cks.AllocateShare()
	for _, pi := range P[1:] {
		cks.GenShare(pi.sk.Get(), zero, encWData[cnum-1], pi.cksShare)
	}
	encOut := ckks.NewCiphertext(params, 1, encWData[0].Level(), params.Scale)
	for _, pi := range P {
		cks.AggregateShares(pi.cksShare, cksCombined, cksCombined)
	}
	cks.KeySwitch(cksCombined, encWData[cnum-1], encOut)

	ptmp := decryptor.DecryptNew(encOut)
	dcw := encoder.Decode(ptmp, slots)

	for j = 0; j < slots; j++ {
		dcw[j] = round(dcw[j])
	}
	var rest uint64 = factorDim - batch*(cnum-1)
	for j = 0; j < rest; j++ {
		wData[batch*(cnum-1)+j] = real(dcw[j])
	}
}
func round(x complex128) complex128 {
	var factor float64
	factor = 100000000
	a := math.Round(real(x)*factor) / factor
	b := math.Round(imag(x)*factor) / factor
	return complex(a, b)
}

func trueIP(a, b []float64, size uint64) (result float64) {
	var res float64 = 0.0
	for i := 0; i < int(size); i++ {
		res += a[i] * b[i]
	}
	return res
}

func caculateAUC(zData [][]float64, wData []float64, factorDim, samplesDim uint64, correctness, auc float64) {
	var i, j uint64
	fmt.Printf("w:")

	for i = 0; i < factorDim; i++ {
		fmt.Printf("%f,", wData[i])
	}
	fmt.Println()

	var TN, FP uint64 = 0, 0
	thetaTN, thetaFP := []float64{}, []float64{}

	for i = 0; i < samplesDim; i++ {
		if zData[i][0] > 0 {
			if trueIP(zData[i], wData, factorDim) < 0 {
				TN++
			}
			thetaTN = append(thetaTN, (zData[i][0] * trueIP(zData[i][1:], wData[1:], factorDim-1)))
		} else {
			if trueIP(zData[i], wData, factorDim) < 0 {
				FP++
			}
			thetaFP = append(thetaFP, (zData[i][0] * trueIP(zData[i][1:], wData[1:], factorDim-1)))
		}
	}
	correctness = 100.0 - (100.0 * float64(FP+TN) / float64(samplesDim))
	fmt.Printf("Correctness: %f %.\n", correctness)

	if len(thetaFP) == 0 || len(thetaTN) == 0 {
		fmt.Println("n_test_yi = 0 : cannot compute AUC")
		auc = 0.0
	} else {
		auc = 0.0
		for i = 0; i < uint64(len(thetaTN)); i++ {
			for j = 0; j < uint64(len(thetaFP)); j++ {
				if thetaFP[j] <= thetaTN[i] {
					auc++
				}
			}
		}
		auc /= float64(len(thetaTN) * len(thetaFP))
		fmt.Printf("AUC: %f \n", auc)
	}
}
