package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math/bits"
	"os"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v3/drlwe"
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
	*drlwe.Combiner

	i        int
	sk       *rlwe.SecretKey
	tsk      *drlwe.ShamirSecretShare
	ssp      *drlwe.ShamirPolynomial
	shamirPk drlwe.ShamirPublicPoint

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
			activePk := make([]drlwe.ShamirPublicPoint, 0)
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
	shamirPks := make([]drlwe.ShamirPublicPoint, 0)
	for i := range P {
		pi := new(party)
		pi.RTGProtocol = drlwe.NewRTGProtocol(params)
		pi.Thresholdizer = drlwe.NewThresholdizer(params)
		pi.i = i
		pi.sk = kg.GenSecretKey()
		params.RingQP().AddLvl(params.QCount()-1, params.PCount()-1, skIdeal.Value, pi.sk.Value, skIdeal.Value)
		pi.tsk = pi.AllocateThresholdSecretShare()
		pi.ssp, err = pi.GenShamirPolynomial(t, pi.sk)
		if err != nil {
			panic(err)
		}
		pi.shamirPk = drlwe.ShamirPublicPoint(i + 1)
		shamirPks = append(shamirPks, pi.shamirPk)
		pi.genTaskQueue = make(chan genTask, k)
		P[i] = pi
	}

	for _, pi := range P {
		pi.Combiner = drlwe.NewCombiner(params, pi.shamirPk, shamirPks, t)
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

	levelQ, levelP := len(params.RingQ().Modulus)-1, len(params.RingP().Modulus)-1
	decompSize := params.DecompPw2(levelQ, levelP) * params.DecompRNS(levelQ, levelP)
	log2bound := bits.Len64(uint64(params.N() * decompSize * (params.N()*3*int(params.NoiseBound()) + 2*3*int(params.NoiseBound()) + params.N()*3)))
	for galEl, rtk := range rtks {
		if !rlwe.RotationKeyIsCorrect(&rtk, galEl, skIdeal, params, log2bound) {
			fmt.Printf("invalid key for galEl=%d\n", galEl)
			os.Exit(1)
		}
	}
	fmt.Println("done")
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
