package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
)

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
// TODO: extract galois parameters type and remove the rotation type from the interface
type RTGProtocol struct {
	drlwe.RTGProtocol

	N           uint64
	galElRotRow uint64
	galElRotCol map[ckks.Rotation][]uint64
}

// NewRotKGProtocol creates a new rotkg object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRotKGProtocol(params *ckks.Parameters) (rtg *RTGProtocol) {

	rtg = new(RTGProtocol)
	rtg.N = params.N()

	rtg.RTGProtocol = *drlwe.NewRTGProtocol(rtg.N, params.Qi(), params.Pi(), params.Sigma())

	rtg.galElRotCol = make(map[ckks.Rotation][]uint64)
	rtg.galElRotCol[ckks.RotationRight] = ring.GenGaloisParams(rtg.N, ckks.GaloisGen)                                     // precompute the galois elements operating left rotation (for left-rotation keys)
	rtg.galElRotCol[ckks.RotationLeft] = ring.GenGaloisParams(rtg.N, ring.ModExp(ckks.GaloisGen, (rtg.N<<1)-1, rtg.N<<1)) // precompute the galois elements operating right rotation (for right-rotation keys)
	rtg.galElRotRow = (rtg.N << 1) - 1

	return rtg
}

// GenShare generates a share in the DBFV rotation key generation protocol
func (rtg *RTGProtocol) GenShare(rotType ckks.Rotation, k uint64, sk *ring.Poly, crp []*ring.Poly, shareOut *drlwe.RTGShare) {
	switch rotType {
	case ckks.RotationRight, ckks.RotationLeft:
		rtg.RTGProtocol.GenShare(sk, rtg.galElRotCol[rotType][k&((rtg.N>>1)-1)], crp, shareOut)
		return
	case ckks.Conjugate:
		rtg.RTGProtocol.GenShare(sk, rtg.galElRotRow, crp, shareOut)
		return
	}
}

// GenCKKSRotationKey populates the input RotationKeys struture with the Switching key computed from the protocol.
func (rtg *RTGProtocol) GenCKKSRotationKey(params *ckks.Parameters, rotType ckks.Rotation, k uint64, share *drlwe.RTGShare, crp []*ring.Poly, rotKey *ckks.RotationKeys) {

	swk := make([][2]*ring.Poly, len(share.Value))
	for i, p := range share.Value {
		swk[i] = [2]*ring.Poly{p, crp[i]}
	}

	rotKey.SetRotKey(params, swk, rotType, k)
}
