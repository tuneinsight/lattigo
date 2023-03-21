package drlwe

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func BenchmarkDRLWE(b *testing.B) {

	thresholdInc := 5

	var err error

	defaultParamsLiteral := rlwe.TestParamsLiteral[:]

	if *flagParamString != "" {
		var jsonParams rlwe.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParamsLiteral = []rlwe.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, paramsLit := range defaultParamsLiteral {

		for _, DefaultNTTFlag := range []bool{true, false} {

			for _, RingType := range []ring.Type{ring.Standard, ring.ConjugateInvariant}[:] {

				paramsLit.DefaultNTTFlag = DefaultNTTFlag
				paramsLit.RingType = RingType

				var params rlwe.Parameters
				if params, err = rlwe.NewParametersFromLiteral(paramsLit); err != nil {
					b.Fatal(err)
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
	}
}

func benchString(opname string, params rlwe.Parameters) string {
	return fmt.Sprintf("%s/LogN=%d/logQP=%d", opname, params.LogN(), params.LogQP())
}

func benchPublicKeyGen(params rlwe.Parameters, b *testing.B) {

	ckg := NewCKGProtocol(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	s1 := ckg.AllocateShare()
	crs, _ := sampling.NewPRNG()

	crp := ckg.SampleCRP(crs)

	b.Run(benchString("PublicKeyGen/Round1/Gen", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ckg.GenShare(sk, crp, s1)
		}
	})

	b.Run(benchString("PublicKeyGen/Round1/Agg", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ckg.AggregateShares(s1, s1, s1)
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
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	ephSk, share1, share2 := rkg.AllocateShare()
	rlk := rlwe.NewRelinearizationKey(params)
	crs, _ := sampling.NewPRNG()

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
			rkg.AggregateShares(share1, share1, share1)
		}
	})

	b.Run(benchString("RelinKeyGen/Finalize", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.GenRelinearizationKey(share1, share2, rlk)
		}
	})
}

func benchRotKeyGen(params rlwe.Parameters, b *testing.B) {

	rtg := NewGKGProtocol(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	share := rtg.AllocateShare()
	crs, _ := sampling.NewPRNG()
	crp := rtg.SampleCRP(crs)

	b.Run(benchString("RotKeyGen/Round1/Gen", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rtg.GenShare(sk, params.GaloisElementForColumnRotationBy(1), crp, share)
		}
	})

	b.Run(benchString("RotKeyGen/Round1/Agg", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rtg.AggregateShares(share, share, share)
		}
	})

	gkey := rlwe.NewGaloisKey(params)
	b.Run(benchString("RotKeyGen/Finalize", params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rtg.GenGaloisKey(share, crp, gkey)
		}
	})
}

func benchThreshold(params rlwe.Parameters, t int, b *testing.B) {

	type Party struct {
		*Thresholdizer
		*Combiner
		gen *ShamirPolynomial
		s   *rlwe.SecretKey
		sk  *rlwe.SecretKey
		tsk *ShamirSecretShare
	}

	shamirPks := make([]ShamirPublicPoint, t)
	for i := range shamirPks {
		shamirPks[i] = ShamirPublicPoint(i + 1)
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

	p.Combiner = NewCombiner(params, shamirPks[0], shamirPks, t)

	b.Run(benchString("Combiner/GenAdditiveShare", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner.GenAdditiveShare(shamirPks, shamirPks[0], p.tsk, p.sk)
		}
	})
}
