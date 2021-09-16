package drlwe

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/rlwe"
)

func Benchmark_DRLWE(b *testing.B) {

	defaultParams := []rlwe.ParametersLiteral{ /*rlwe.TestPN12QP109, rlwe.TestPN13QP218, */ rlwe.TestPN14QP438 /* rlwe.TestPN15QP880*/}
	thresholdInc := 1

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

		// Varying t
		N := 20
		for t := 2; t < N; t += thresholdInc {
			benchThreshold(params, t, b)
		}

	}
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

	b.Run(testString("Thresholdizer/GenShamirPolynomial/", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.gen, _ = p.Thresholdizer.GenShamirPolynomial(t, p.s)
		}
	})

	shamirShare := p.Thresholdizer.AllocateThresholdSecretShare()

	b.Run(testString("Thresholdizer/GenShamirSecretShare/", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Thresholdizer.GenShamirSecretShare(shamirPks[0], p.gen, shamirShare)
		}
	})

	b.Run(testString("Thresholdizer/AggregateShares/", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Thresholdizer.AggregateShares(shamirShare, shamirShare, shamirShare)
		}
	})

	p.Combiner = NewCombiner(params, t)
	p.CachedCombiner = NewCachedCombiner(params, t)

	p.CachedCombiner.Precompute(shamirPks, shamirPks[0])

	b.Run(testString("Combiner/GenAdditiveShare/", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner.GenAdditiveShare(shamirPks, shamirPks[0], p.tsk, p.sk)
		}
	})

	b.Run(testString("CombinerCached/GenAdditiveShare/", params)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.CachedCombiner.GenAdditiveShare(shamirPks, shamirPks[0], p.tsk, p.sk)
		}
	})
}
