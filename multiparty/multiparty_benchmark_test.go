package multiparty

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func BenchmarkMultiParty(b *testing.B) {

	thresholdInc := 5

	var err error

	defaultParamsLiteral := testInsecure

	if *flagParamString != "" {
		var jsonParams TestParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParamsLiteral = []TestParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, paramsLit := range defaultParamsLiteral {

		for _, NTTFlag := range []bool{true, false} {

			for _, RingType := range []ring.Type{ring.Standard, ring.ConjugateInvariant}[:] {

				paramsLit.NTTFlag = NTTFlag
				paramsLit.RingType = RingType

				var params rlwe.Parameters
				if params, err = rlwe.NewParametersFromLiteral(paramsLit.ParametersLiteral); err != nil {
					b.Fatal(err)
				}

				levelQ := params.MaxLevelQ()
				levelP := params.MaxLevelP()
				bpw2 := paramsLit.BaseTwoDecomposition

				benchPublicKeyGen(params, levelQ, levelP, bpw2, b)
				benchRelinearizationKeyGen(params, levelQ, levelP, bpw2, b)
				benchRotKeyGen(params, levelQ, levelP, bpw2, b)

				// Varying t
				for t := 2; t <= 19; t += thresholdInc {
					benchThreshold(params, levelQ, levelP, bpw2, t, b)
				}
			}
		}
	}
}

func benchString(params rlwe.Parameters, opname string, levelQ, levelP, bpw2 int) string {
	return fmt.Sprintf("%s/logN=%d/#Qi=%d/#Pi=%d/Pw2=%d/NTT=%t/RingType=%s",
		opname,
		params.LogN(),
		levelQ+1,
		levelP+1,
		bpw2,
		params.NTTFlag(),
		params.RingType())
}

func benchPublicKeyGen(params rlwe.Parameters, levelQ, levelP, bpw2 int, b *testing.B) {

	ckg := NewPublicKeyGenProtocol(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	s1 := ckg.AllocateShare()
	crs, _ := sampling.NewPRNG()

	crp := ckg.SampleCRP(crs)

	b.Run(benchString(params, "PublicKeyGen/Round1/Gen", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ckg.GenShare(sk, crp, &s1)
		}
	})

	b.Run(benchString(params, "PublicKeyGen/Round1/Agg", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ckg.AggregateShares(s1, s1, &s1)
		}
	})

	pk := rlwe.NewPublicKey(params)
	b.Run(benchString(params, "PublicKeyGen/Finalize", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ckg.GenPublicKey(s1, crp, pk)
		}
	})
}

func benchRelinearizationKeyGen(params rlwe.Parameters, levelQ, levelP, bpw2 int, b *testing.B) {

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(levelQ), LevelP: utils.Pointy(levelP), BaseTwoDecomposition: utils.Pointy(bpw2)}

	rkg := NewRelinearizationKeyGenProtocol(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	ephSk, share1, share2 := rkg.AllocateShare(evkParams)
	rlk := rlwe.NewRelinearizationKey(params, evkParams)
	crs, _ := sampling.NewPRNG()

	crp := rkg.SampleCRP(crs, evkParams)

	b.Run(benchString(params, "RelinearizationKeyGen/GenRound1", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.GenShareRoundOne(sk, crp, ephSk, &share1)
		}
	})

	b.Run(benchString(params, "RelinearizationKeyGen/GenRound2", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.GenShareRoundTwo(ephSk, sk, share1, &share2)
		}
	})

	b.Run(benchString(params, "RelinearizationKeyGen/Agg", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.AggregateShares(share1, share1, &share1)
		}
	})

	b.Run(benchString(params, "RelinearizationKeyGen/Finalize", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rkg.GenRelinearizationKey(share1, share2, rlk)
		}
	})
}

func benchRotKeyGen(params rlwe.Parameters, levelQ, levelP, bpw2 int, b *testing.B) {

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(levelQ), LevelP: utils.Pointy(levelP), BaseTwoDecomposition: utils.Pointy(bpw2)}

	rtg := NewGaloisKeyGenProtocol(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	share := rtg.AllocateShare(evkParams)
	crs, _ := sampling.NewPRNG()
	crp := rtg.SampleCRP(crs, evkParams)

	b.Run(benchString(params, "RotKeyGen/Round1/Gen", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rtg.GenShare(sk, params.GaloisElement(1), crp, &share)
		}
	})

	b.Run(benchString(params, "RotKeyGen/Round1/Agg", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rtg.AggregateShares(share, share, &share)
		}
	})

	gkey := rlwe.NewGaloisKey(params, evkParams)
	b.Run(benchString(params, "RotKeyGen/Finalize", levelQ, levelP, bpw2), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rtg.GenGaloisKey(share, crp, gkey)
		}
	})
}

func benchThreshold(params rlwe.Parameters, levelQ, levelP, bpw2 int, t int, b *testing.B) {

	type Party struct {
		Thresholdizer
		Combiner
		gen ShamirPolynomial
		s   *rlwe.SecretKey
		sk  *rlwe.SecretKey
		tsk ShamirSecretShare
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

	b.Run(benchString(params, "Thresholdizer/GenShamirPolynomial", levelQ, levelP, bpw2)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.gen, _ = p.Thresholdizer.GenShamirPolynomial(t, p.s)
		}
	})

	shamirShare := p.Thresholdizer.AllocateThresholdSecretShare()

	b.Run(benchString(params, "Thresholdizer/GenShamirSecretShare", levelQ, levelP, bpw2)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Thresholdizer.GenShamirSecretShare(shamirPks[0], p.gen, &shamirShare)
		}
	})

	b.Run(benchString(params, "Thresholdizer/AggregateShares", levelQ, levelP, bpw2)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Thresholdizer.AggregateShares(shamirShare, shamirShare, &shamirShare)
		}
	})

	p.Combiner = NewCombiner(params, shamirPks[0], shamirPks, t)

	b.Run(benchString(params, "Combiner/GenAdditiveShare", levelQ, levelP, bpw2)+fmt.Sprintf("/threshold=%d", t), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.Combiner.GenAdditiveShare(shamirPks, shamirPks[0], p.tsk, p.sk)
		}
	})
}
