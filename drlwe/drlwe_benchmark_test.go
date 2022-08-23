package drlwe

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

func BenchmarkDRLWE(b *testing.B) {

	defaultParams := []rlwe.ParametersLiteral{rlwe.TestPN12QP109, rlwe.TestPN13QP218, rlwe.TestPN14QP438, rlwe.TestPN15QP880}
	thresholdInc := 5

	if testing.Short() {
		defaultParams = defaultParams[:2]
		thresholdInc = 5
	}

	if *flagParamString != "" {
		var jsonParams rlwe.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []rlwe.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {
		params, err := rlwe.NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}

		benchPublicKeyGen(params, b)
		benchRelinKeyGen(params, b)
		benchRotKeyGen(params, b)

		// Varying t
		for t := 2; t <= 19; t += thresholdInc {
			benchThreshold(params, t, b)
		}

	}
}

func benchString(opname string, params rlwe.Parameters) string {
	return fmt.Sprintf("%s/LogN=%d/logQP=%d", opname, params.LogN(), params.LogQP())
}

func benchPublicKeyGen(params rlwe.Parameters, b *testing.B) {

	ckg := NewCKGProtocol(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKey()
	s1 := ckg.AllocateShare()
	crs, _ := utils.NewPRNG()

	crp := ckg.SampleCRP(crs)

	b.Run(benchString("PublicKeyGen/Round1/Gen", params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			ckg.GenShare(sk, crp, s1)
		}
	})

	b.Run(benchString("PublicKeyGen/Round1/Agg", params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			ckg.AggregateShare(s1, s1, s1)
		}
	})

	pk := rlwe.NewPublicKey(params)
	b.Run(benchString("PublicKeyGen/Finalize", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ckg.GenPublicKey(s1, crp, pk)
		}
	})
}

func benchRelinKeyGen(params rlwe.Parameters, b *testing.B) {

	rkg := NewRKGProtocol(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKey()
	ephSk, share1, share2 := rkg.AllocateShare()
	rlk := rlwe.NewRelinKey(params, 2)
	crs, _ := utils.NewPRNG()

	crp := rkg.SampleCRP(crs)

	b.Run(benchString("RelinKeyGen/GenRound1", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.GenShareRoundOne(sk, crp, ephSk, share1)
		}
	})

	b.Run(benchString("RelinKeyGen/GenRound2", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.GenShareRoundTwo(ephSk, sk, share1, share2)
		}
	})

	b.Run(benchString("RelinKeyGen/Agg", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.AggregateShare(share1, share1, share1)
		}
	})

	b.Run(benchString("RelinKeyGen/Finalize", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.GenRelinearizationKey(share1, share2, rlk)
		}
	})
}

func benchRotKeyGen(params rlwe.Parameters, b *testing.B) {

	rtg := NewRTGProtocol(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKey()
	share := rtg.AllocateShare()
	crs, _ := utils.NewPRNG()
	crp := rtg.SampleCRP(crs)

	b.Run(benchString("RotKeyGen/Round1/Gen", params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			rtg.GenShare(sk, params.GaloisElementForRowRotation(), crp, share)
		}
	})

	b.Run(benchString("RotKeyGen/Round1/Agg", params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			rtg.AggregateShare(share, share, share)
		}
	})

	rotKey := rlwe.NewSwitchingKey(params, params.QCount()-1, params.PCount()-1)
	b.Run(benchString("RotKeyGen/Finalize", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rtg.GenRotationKey(share, crp, rotKey)
		}
	})
}

func benchThreshold(params rlwe.Parameters, t int, b *testing.B) {

	type Party struct {
		*Thresholdizer
		Combiner
		*CachedCombiner
		gen *ShamirPolynomial
		s   *rlwe.SecretKey
		sk  *rlwe.SecretKey
		tsk *ShamirSecretShare
	}

	shamirPks := make([]ShamirPublicKey, t)
	for i := range shamirPks {
		shamirPks[i] = ShamirPublicKey(i + 1)
	}

	p := new(Party)
	p.s = rlwe.NewSecretKey(params)
	p.Thresholdizer = NewThresholdizer(params)
	p.tsk = p.Thresholdizer.AllocateThresholdSecretShare()
	p.sk = rlwe.NewSecretKey(params)

	b.Run(benchString("Thresholdizer/GenShamirPolynomial", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.gen, _ = p.Thresholdizer.GenShamirPolynomial(t, p.s)
		}
	})

	shamirShare := p.Thresholdizer.AllocateThresholdSecretShare()

	b.Run(benchString("Thresholdizer/GenShamirSecretShare", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Thresholdizer.GenShamirSecretShare(shamirPks[0], p.gen, shamirShare)
		}
	})

	b.Run(benchString("Thresholdizer/AggregateShares", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Thresholdizer.AggregateShares(shamirShare, shamirShare, shamirShare)
		}
	})

	p.Combiner = NewCombiner(params, t)
	p.CachedCombiner = NewCachedCombiner(params, t)

	p.CachedCombiner.Precompute(shamirPks, shamirPks[0])

	b.Run(benchString("Combiner/GenAdditiveShare", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner.GenAdditiveShare(shamirPks, shamirPks[0], p.tsk, p.sk)
		}
	})

	b.Run(benchString("CombinerCached/GenAdditiveShare", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.CachedCombiner.GenAdditiveShare(shamirPks, shamirPks[0], p.tsk, p.sk)
		}
	})
}
