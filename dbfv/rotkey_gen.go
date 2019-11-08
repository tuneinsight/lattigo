package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// rotkg is the structure storing the parameters for the collective rotation-keys generation.
type RTGProtocol struct {
	bfvContext *bfv.BfvContext

	gaussianSampler *ring.KYSampler

	galElRotRow uint64
	galElRotCol map[bfv.Rotation][]uint64

	tmpSwitchKey [][2]*ring.Poly
	tmpPoly      *ring.Poly
}

type RTGShare struct {
	Type  bfv.Rotation
	K     uint64
	Value []*ring.Poly
}

func (rtg *RTGProtocol) AllocateShare() (rtgShare RTGShare) {
	rtgShare.Value = make([]*ring.Poly, rtg.bfvContext.Beta())
	for i := uint64(0); i < rtg.bfvContext.Beta(); i++ {
		rtgShare.Value[i] = rtg.bfvContext.ContextKeys().NewPoly()
	}
	return
}

// newrotkg creates a new rotkg object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRotKGProtocol(bfvContext *bfv.BfvContext) (rtg *RTGProtocol) {

	rtg = new(RTGProtocol)
	rtg.bfvContext = bfvContext

	rtg.gaussianSampler = bfvContext.GaussianSampler()

	rtg.tmpSwitchKey = make([][2]*ring.Poly, rtg.bfvContext.Beta())
	for i := range rtg.tmpSwitchKey {
		rtg.tmpSwitchKey[i][0] = bfvContext.ContextKeys().NewPoly()
		rtg.tmpSwitchKey[i][1] = bfvContext.ContextKeys().NewPoly()
	}

	rtg.tmpPoly = bfvContext.ContextKeys().NewPoly()

	N := bfvContext.ContextKeys().N
	mask := (N << 1) - 1

	rtg.galElRotCol = make(map[bfv.Rotation][]uint64)
	for _, rotType := range []bfv.Rotation{bfv.RotationLeft, bfv.RotationRight} {
		rtg.galElRotCol[rotType] = make([]uint64, N>>1)
		rtg.galElRotCol[rotType][0] = 1

		gen := bfv.GaloisGen
		if rotType == bfv.RotationRight {
			gen = ring.ModExp(gen, (N<<1)-1, N<<1)
		}

		for i := uint64(1); i < N>>1; i++ {
			rtg.galElRotCol[rotType][i] = (rtg.galElRotCol[rotType][i-1] * gen) & mask
		}
	}

	rtg.galElRotCol[bfv.RotationRight] = make([]uint64, N>>1)

	rtg.galElRotRow = (N << 1) - 1

	return rtg
}

// GenShare is the first and unique round of the rotkg protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + (pi(s_i) - s_i) + e]
//
// and broadcasts it to the other j-1 parties. The protocol must be repeated for each desired rotation.
func (rtg *RTGProtocol) GenShare(rotType bfv.Rotation, k uint64, sk *ring.Poly, crp []*ring.Poly, shareOut *RTGShare) {
	shareOut.Type = rotType
	shareOut.K = k
	switch rotType {
	case bfv.RotationRight, bfv.RotationLeft:
		rtg.genShare(sk, rtg.galElRotCol[rotType][k&((rtg.bfvContext.N()>>1)-1)], crp, shareOut.Value)
		return
	case bfv.RotationRow:
		rtg.genShare(sk, rtg.galElRotRow, crp, shareOut.Value)
		return
	}
}

// genswitchkey is a generic method to generate the public-share of the collective rotation-key.
func (rtg *RTGProtocol) genShare(sk *ring.Poly, galEl uint64, crp []*ring.Poly, evakey []*ring.Poly) {

	contextKeys := rtg.bfvContext.ContextKeys()

	ring.PermuteNTT(sk, galEl, rtg.tmpPoly)
	contextKeys.Sub(rtg.tmpPoly, sk, rtg.tmpPoly)

	for _, pj := range rtg.bfvContext.KeySwitchPrimes() {
		contextKeys.MulScalar(rtg.tmpPoly, pj, rtg.tmpPoly)
	}

	contextKeys.InvMForm(rtg.tmpPoly, rtg.tmpPoly)

	var index uint64

	for i := uint64(0); i < rtg.bfvContext.Beta(); i++ {

		// e
		evakey[i] = rtg.gaussianSampler.SampleNTTNew()

		// a is the CRP

		// e + sk_in * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0
		for j := uint64(0); j < rtg.bfvContext.Alpha(); j++ {

			index = i*rtg.bfvContext.Alpha() + j

			for w := uint64(0); w < contextKeys.N; w++ {
				evakey[i].Coeffs[index][w] = ring.CRed(evakey[i].Coeffs[index][w]+rtg.tmpPoly.Coeffs[index][w], contextKeys.Modulus[index])
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= uint64(len(rtg.bfvContext.ContextQ().Modulus)-1) {
				break
			}
		}

		// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
		contextKeys.MulCoeffsMontgomeryAndSub(crp[i], sk, evakey[i])
		contextKeys.MForm(evakey[i], evakey[i])

	}

	return
}

// Aggregate is the second part of the unique round of the rotkg protocol. Uppon receiving the j-1 public shares,
// each party computes  :
//
// [sum(a*a_j + (pi(a_j) - a_j) + e_j), a]
func (rtg *RTGProtocol) Aggregate(share1, share2, shareOut RTGShare) {
	contextKeys := rtg.bfvContext.ContextKeys()

	for i := uint64(0); i < rtg.bfvContext.Beta(); i++ {
		contextKeys.Add(share1.Value[i], share2.Value[i], shareOut.Value[i])
	}
}

func (rtg *RTGProtocol) Finalize(share RTGShare, crp []*ring.Poly, rotKey *bfv.RotationKeys) {

	k := share.K & ((rtg.bfvContext.N() >> 1) - 1)

	for i := uint64(0); i < rtg.bfvContext.Beta(); i++ {
		rtg.tmpSwitchKey[i][0].Copy(share.Value[i])
		rtg.bfvContext.ContextKeys().MForm(crp[i], rtg.tmpSwitchKey[i][1])
	}

	switch share.Type {
	case bfv.RotationLeft:
		rotKey.SetRotColLeft(rtg.tmpSwitchKey, k)
	case bfv.RotationRight:
		rotKey.SetRotColRight(rtg.tmpSwitchKey, k)
	case bfv.RotationRow:
		rotKey.SetRotRow(rtg.tmpSwitchKey)
	}
}
