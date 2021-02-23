package bfv

import "github.com/ldsec/lattigo/v2/rlwe"

type SecretKey struct{ rlwe.SecretKey }

type PublicKey struct{ rlwe.PublicKey }

type SwitchingKey struct{ rlwe.SwitchingKey }

type RelinearizationKey struct{ rlwe.RelinearizationKey }

type RotationKeySet struct{ rlwe.RotationKeySet }
type EvaluationKey struct {
	Rlk  *RelinearizationKey
	Rtks *RotationKeySet
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {
	return &PublicKey{*rlwe.NewPublicKey(params.N(), params.QPiCount())}
}

func NewSwitchingKey(params *Parameters) *SwitchingKey {
	return &SwitchingKey{*rlwe.NewSwitchingKey(params.N(), params.QPiCount(), params.Beta())}
}

func NewRelinKey(params *Parameters, maxRelinDegree uint64) *RelinearizationKey {
	return &RelinearizationKey{*rlwe.NewRelinKey(maxRelinDegree, params.N(), params.QPiCount(), params.Beta())}
}

func NewRotationKeySet(params *Parameters, galoisElements []uint64) *RotationKeySet {
	return &RotationKeySet{*rlwe.NewRotationKeySet(galoisElements, params.N(), params.QPiCount(), params.Beta())}
}
