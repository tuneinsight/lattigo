package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"os"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

type genTask struct {
	group     []*party
	galoisEls []uint64
}

type genTaskResult struct {
	galEl uint64

	rtgShare *drlwe.RTGShare
}

type party struct {
	*drlwe.RTGProtocol
	*drlwe.Thresholdizer
	drlwe.Combiner

	i        int
	sk       *rlwe.SecretKey
	tsk      *drlwe.ShamirSecretShare
	ssp      *drlwe.ShamirPolynomial
	shamirPk drlwe.ShamirPublicKey

	genTaskQueue chan genTask
}

func (p *party) String() string {
	return fmt.Sprintf("Party#%d", p.i)
}

type cloud struct {
	*drlwe.RTGProtocol

	aggTaskQueue chan genTaskResult
	finDone      chan struct {
		galEl uint64
		rtk   rlwe.SwitchingKey
	}
}

var crp drlwe.RTGCRP

func (p *party) Run(wg *sync.WaitGroup, params rlwe.Parameters, N int, P []*party, C *cloud) {

	var nShares, nTasks int
	var start time.Time
	var cpuTime time.Duration
	var byteSent int
	for task := range p.genTaskQueue {

		start = time.Now()
		var sk *rlwe.SecretKey
		t := len(task.group)
		if t == N {
			sk = p.sk
		} else {
			activePk := make([]drlwe.ShamirPublicKey, 0)
			for _, pi := range task.group {
				activePk = append(activePk, pi.shamirPk)
			}
			sk = rlwe.NewSecretKey(params)
			p.GenAdditiveShare(activePk, p.shamirPk, p.tsk, sk)
		}

		for _, galEl := range task.galoisEls {
			rtgShare := p.AllocateShare()

			p.GenShare(sk, galEl, crp, rtgShare)
			C.aggTaskQueue <- genTaskResult{galEl: galEl, rtgShare: rtgShare}
			nShares++
			byteSent += len(rtgShare.Value) * len(rtgShare.Value[0]) * rtgShare.Value[0][0].GetDataLen64(false)
		}
		nTasks++
		cpuTime += time.Since(start)
	}
	wg.Done()
	fmt.Printf("\tParty %d finished generating %d shares of %d tasks in %s, sent %s\n", p.i, nShares, nTasks, cpuTime, byteCountSI(byteSent))
}

func (c *cloud) Run(galEls []uint64, params rlwe.Parameters, t int) {

	shares := make(map[uint64]*struct {
		share  *drlwe.RTGShare
		needed int
	}, len(galEls))
	for _, galEl := range galEls {
		shares[galEl] = &struct {
			share  *drlwe.RTGShare
			needed int
		}{c.AllocateShare(), t}
	}

	var i int
	var cpuTime time.Duration
	var byteRecv int
	for task := range c.aggTaskQueue {
		start := time.Now()
		acc := shares[task.galEl]
		c.AggregateShare(acc.share, task.rtgShare, acc.share)
		acc.needed--
		if acc.needed == 0 {
			rtk := rlwe.NewSwitchingKey(params, params.MaxLevel(), params.PCount()-1)
			c.GenRotationKey(acc.share, crp, rtk)
			c.finDone <- struct {
				galEl uint64
				rtk   rlwe.SwitchingKey
			}{galEl: task.galEl, rtk: *rtk}
		}
		i++
		cpuTime += time.Since(start)
		byteRecv += len(acc.share.Value) * len(acc.share.Value[0]) * acc.share.Value[0][0].GetDataLen64(false)
	}
	close(c.finDone)
	fmt.Printf("\tCloud finished aggregating %d shares in %s, received %s\n", i, cpuTime, byteCountSI(byteRecv))

}

func getSubGroups(P []*party, t, k int) [][]*party {
	if t == len(P) {
		return [][]*party{P}
	}
	if t > len(P) {
		panic("t > len(P)")
	}

	groups := [][]*party{}
	for i := 0; i < k; i++ {
		start := (i * t) % len(P)
		end := ((i + 1) * t) % len(P)
		switch {
		case i > 0 && start == 0 && end == t:
			return groups
		case start > end:
			group := make([]*party, t)
			copy(group, P[0:end])
			copy(group[end:], P[start:])
			groups = append(groups, group)
		default:
			groups = append(groups, P[start:end])
		}
	}

	return groups
}

func getTasks(galEls []uint64, groups [][]*party) []genTask {
	tasks := make([]genTask, len(groups))
	for i := range tasks {
		tasks[i].group = groups[i]
	}
	for i, galEl := range galEls {
		tasks[i%len(groups)].galoisEls = append(tasks[i%len(groups)].galoisEls, galEl)
	}
	return tasks
}

var flagN = flag.Int("N", 3, "the number of parties")
var flagT = flag.Int("t", 2, "the threshold")
var flagO = flag.Int("o", 0, "the number of online parties")
var flagK = flag.Int("k", 10, "number of rotation keys to generate")
var flagDefaultParams = flag.Int("params", 3, "default param set to use")
var flagJSONParams = flag.String("json", "", "the JSON encoded parameter set to use")

func main() {

	flag.Parse()

	if *flagDefaultParams >= len(rlwe.DefaultParams) {
		panic("invalid default parameter set")
	}

	paramsDef := rlwe.DefaultParams[*flagDefaultParams]

	if *flagJSONParams != "" {
		if err := json.Unmarshal([]byte(*flagJSONParams), &paramsDef); err != nil {
			panic(err)
		}
	}

	params, err := rlwe.NewParametersFromLiteral(paramsDef)
	if err != nil {
		panic(err)
	}

	if *flagN < 2 {
		panic("-N should be >= 2")
	}
	N := *flagN

	if *flagT > N {
		panic("-t should be <= N")
	}
	t := *flagT

	var o int
	if *flagO <= 0 {
		o = N
	} else {
		o = *flagO
	}

	if *flagK < 1 {
		panic("-k should be >= 1")
	}
	k := *flagK

	galEls := make([]uint64, k)
	for i := range galEls {
		galEls[i] = params.GaloisElementForColumnRotationBy(i + 1)
	}

	fmt.Printf("Starting for N=%d, t=%d\n", N, t)
	fmt.Printf("LogN=%d, LogQP=%d, L=%d, k=%d\n", params.LogN(), params.LogQP(), params.QPCount(), k)

	kg := rlwe.NewKeyGenerator(params)

	crs, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	wg := new(sync.WaitGroup)
	C := &cloud{
		RTGProtocol:  drlwe.NewRTGProtocol(params),
		aggTaskQueue: make(chan genTaskResult, len(galEls)*N),
		finDone: make(chan struct {
			galEl uint64
			rtk   rlwe.SwitchingKey
		}, len(galEls)),
	}
	P := make([]*party, N)
	skIdeal := rlwe.NewSecretKey(params)
	for i := range P {
		pi := new(party)
		pi.RTGProtocol = drlwe.NewRTGProtocol(params)
		pi.Thresholdizer = drlwe.NewThresholdizer(params)
		pi.Combiner = drlwe.NewCombiner(params, t)
		pi.i = i
		pi.sk = kg.GenSecretKey()
		params.RingQP().AddLvl(params.QCount()-1, params.PCount()-1, skIdeal.Value, pi.sk.Value, skIdeal.Value)
		pi.tsk = pi.AllocateThresholdSecretShare()
		pi.ssp, err = pi.GenShamirPolynomial(t, pi.sk)
		if err != nil {
			panic(err)
		}
		pi.shamirPk = drlwe.ShamirPublicKey(i + 1)
		pi.genTaskQueue = make(chan genTask, k)
		P[i] = pi
	}

	if t < N {
		fmt.Println("Performing threshold setup")
		shares := make(map[*party]map[*party]*drlwe.ShamirSecretShare, len(P))
		for _, pi := range P {

			shares[pi] = make(map[*party]*drlwe.ShamirSecretShare)

			for _, pj := range P {
				shares[pi][pj] = pi.AllocateThresholdSecretShare()
				pi.GenShamirSecretShare(pj.shamirPk, pi.ssp, shares[pi][pj])
			}
		}

		for _, pi := range P {
			for _, pj := range P {
				pi.AggregateShares(pi.tsk, shares[pj][pi], pi.tsk)
			}
		}
	}

	P = P[:o]

	groups := getSubGroups(P, t, k)
	fmt.Printf("Generating %d rotation keys with %d parties in %d groups\n", len(galEls), len(P), len(groups))

	crp = P[0].SampleCRP(crs)

	go C.Run(galEls, params, t)
	for _, pi := range P {
		go pi.Run(wg, params, N, P, C)
	}

	wg.Add(len(P))
	start := time.Now()
	tasks := getTasks(galEls, groups)
	for _, task := range tasks {
		for _, p := range task.group {
			p.genTaskQueue <- task
		}
	}
	for _, pi := range P {
		close(pi.genTaskQueue)
	}
	wg.Wait()
	close(C.aggTaskQueue)

	rtks := make(map[uint64]rlwe.SwitchingKey)
	for task := range C.finDone {
		rtks[task.galEl] = task.rtk
	}
	fmt.Printf("Generation of %d keys completed in %s\n", len(rtks), time.Since(start))

	fmt.Printf("Checking the keys... ")
	for galEl, rtk := range rtks {
		if !verifyKey(rtk, galEl, skIdeal, params) {
			fmt.Printf("invalid key for galEl=%d\n", galEl)
			os.Exit(1)
		}
	}
	fmt.Println("done")
}

func verifyKey(swk rlwe.SwitchingKey, galEl uint64, skIdeal *rlwe.SecretKey, params rlwe.Parameters) bool {
	skIn := skIdeal.CopyNew()
	skOut := skIdeal.CopyNew()
	galElInv := ring.ModExp(galEl, uint64(4*params.N()-1), uint64(4*params.N()))
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	ringQ, ringP, ringQP := params.RingQ(), params.RingP(), params.RingQP()
	decompPw2 := params.DecompPw2(levelQ, levelP)

	ringQ.PermuteNTT(skIdeal.Value.Q, galElInv, skOut.Value.Q)
	if ringP != nil {
		ringP.PermuteNTT(skIdeal.Value.P, galElInv, skOut.Value.P)
	}

	// Decrypts
	// [-asIn + w*P*sOut + e, a] + [asIn]
	for i := range swk.Value {
		for j := range swk.Value[i] {
			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, swk.Value[i][j].Value[1], skOut.Value, swk.Value[i][j].Value[0])
		}
	}

	// Sums all basis together (equivalent to multiplying with CRT decomposition of 1)
	// sum([1]_w * [RNS*PW2*P*sOut + e]) = PWw*P*sOut + sum(e)
	for i := range swk.Value { // RNS decomp
		if i > 0 {
			for j := range swk.Value[i] { // PW2 decomp
				ringQP.AddLvl(levelQ, levelP, swk.Value[0][j].Value[0], swk.Value[i][j].Value[0], swk.Value[0][j].Value[0])
			}
		}
	}

	if levelP != -1 {
		// sOut * P
		ringQ.MulScalarBigint(skIn.Value.Q, ringP.ModulusAtLevel[levelP], skIn.Value.Q)
	}

	log2Bound := bits.Len64(uint64(params.N() * len(swk.Value) * len(swk.Value[0]) * (params.N()*3*int(math.Floor(rlwe.DefaultSigma*6)) + 2*3*int(math.Floor(rlwe.DefaultSigma*6)) + params.N()*3)))
	for i := 0; i < decompPw2; i++ {

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(swk.Value[0][i].Value[0].Q, skIn.Value.Q, swk.Value[0][i].Value[0].Q)

		// Checks that the error is below the bound
		// Worst error bound is N * floor(6*sigma) * #Keys
		ringQP.InvNTTLvl(levelQ, levelP, swk.Value[0][i].Value[0], swk.Value[0][i].Value[0])
		ringQP.InvMFormLvl(levelQ, levelP, swk.Value[0][i].Value[0], swk.Value[0][i].Value[0])

		// Worst bound of inner sum
		// N*#Keys*(N * #Parties * floor(sigma*6) + #Parties * floor(sigma*6) + N * #Parties  +  #Parties * floor(6*sigma))

		if log2Bound < log2OfInnerSum(levelQ, ringQ, swk.Value[0][i].Value[0].Q) {
			return false
		}

		if levelP != -1 {
			if log2Bound < log2OfInnerSum(levelP, ringP, swk.Value[0][i].Value[0].P) {
				return false
			}
		}

		// sOut * P * PW2
		ringQ.MulScalar(skIn.Value.Q, 1<<params.Pow2Base(), skIn.Value.Q)
	}

	return true
}

// Returns the ceil(log2) of the sum of the absolute value of all the coefficients
func log2OfInnerSum(level int, ringQ *ring.Ring, poly *ring.Poly) (logSum int) {
	sumRNS := make([]uint64, level+1)
	var sum uint64
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		qiHalf := qi >> 1
		coeffs := poly.Coeffs[i]
		sum = 0

		for j := 0; j < ringQ.N; j++ {

			v := coeffs[j]

			if v >= qiHalf {
				sum = ring.CRed(sum+qi-v, qi)
			} else {
				sum = ring.CRed(sum+v, qi)
			}
		}

		sumRNS[i] = sum
	}

	var smallNorm = true
	for i := 1; i < level+1; i++ {
		smallNorm = smallNorm && (sumRNS[0] == sumRNS[i])
	}

	if !smallNorm {
		var qi uint64
		var crtReconstruction *big.Int

		sumBigInt := ring.NewUint(0)
		QiB := new(big.Int)
		tmp := new(big.Int)
		modulusBigint := ring.NewUint(1)

		for i := 0; i < level+1; i++ {

			qi = ringQ.Modulus[i]
			QiB.SetUint64(qi)

			modulusBigint.Mul(modulusBigint, QiB)

			crtReconstruction = new(big.Int)
			crtReconstruction.Quo(ringQ.ModulusAtLevel[len(ringQ.Modulus)-1], QiB)
			tmp.ModInverse(crtReconstruction, QiB)
			tmp.Mod(tmp, QiB)
			crtReconstruction.Mul(crtReconstruction, tmp)

			sumBigInt.Add(sumBigInt, tmp.Mul(ring.NewUint(sumRNS[i]), crtReconstruction))
		}

		sumBigInt.Mod(sumBigInt, modulusBigint)

		logSum = sumBigInt.BitLen()
	} else {
		logSum = bits.Len64(sumRNS[0])
	}

	return
}

func byteCountSI(b int) string {
	const unit = 1000
	if b < unit {
		return fmt.Sprintf("%d B", b)
	}
	div, exp := int64(unit), 0
	for n := b / unit; n >= unit; n /= unit {
		div *= unit
		exp++
	}
	return fmt.Sprintf("%.1f %cB",
		float64(b)/float64(div), "kMGTPE"[exp])
}
