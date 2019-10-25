package dckks

import (
	"fmt"
	"testing"

	//"github.com/ldsec/lattigo/ckks"
	//"github.com/ldsec/lattigo/ring"
)


func Benchmark_CKG(b *testing.B) {

}

func Benchmark_EKG(b *testing.B) {
	
}

func Benchmark_EKGNaive(b *testing.B) {

}

func Benchmark_RKG(b *testing.B) {

}

func Benchmark_CKS(b *testing.B) {

}

func Benchmark_PCKS(b *testing.B) {
	
}

func Benchmark_Refresh(b *testing.B) {

	parties := params.parties
	ckksContext := params.ckksContext
	sk0Shards := params.sk0Shards

	crp := params.ckksContext.ContextQ().NewPoly()
	ciphertext := params.ckksContext.NewRandomCiphertext(1, 2, params.ckksContext.Scale())

	refreshShares := make([]*RefreshShares, params.parties)

	startLevel := uint64(2)

	b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round0", 
		parties, 
		ckksContext.LogN(), 
		ckksContext.LogQ(), 
		ckksContext.Levels(), 
		ckksContext.Scale()),
		func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			refreshShares[0] = GenRefreshShares(sk0Shards[0], startLevel, parties, ckksContext, ciphertext.Value()[0], crp)
		}
	})

	for i := uint64(1); i < parties; i++ {
		refreshShares[i] = GenRefreshShares(sk0Shards[i], startLevel, parties, ckksContext, ciphertext.Value()[0], crp)
	}

	b.Run(fmt.Sprintf("parties=%d/logN=%d/logQ=%d/levels=%d/scale=%f/Round1", 
		parties, 
		ckksContext.LogN(), 
		ckksContext.LogQ(), 
		ckksContext.Levels(), 
		ckksContext.Scale()),
		func(b *testing.B) {
		for i := 0; i < b.N; i++ {

			Refresh(ciphertext, refreshShares, ckksContext, crp)

			b.StopTimer()
			ciphertext.Value()[0].Coeffs = ciphertext.Value()[0].Coeffs[:startLevel+1]
			b.StartTimer()
		}
	})
}
