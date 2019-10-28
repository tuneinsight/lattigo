package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// rotkg is the structure storing the parameters for the collective rotation-keys generation.
type RotKGProtocol struct {
	bfvContext *bfv.BfvContext

	gaussianSampler *ring.KYSampler

	gen    uint64
	genInv uint64

	galElRotRow      uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64

	rot_col_L map[uint64][][2]*ring.Poly
	rot_col_R map[uint64][][2]*ring.Poly
	rot_row   [][2]*ring.Poly

	polypool *ring.Poly
}

type RotKGShareRotColLeft []*ring.Poly
type RotKGShareRotColRight []*ring.Poly
type RotKGShareConjugate []*ring.Poly

func (rotkg *RotKGProtocol) AllocateShareRotColLeft() (r RotKGShareRotColLeft) {
	r = make([]*ring.Poly, rotkg.bfvContext.Beta())
	for i := uint64(0); i < rotkg.bfvContext.Beta(); i++ {
		r[i] = rotkg.bfvContext.ContextKeys().NewPoly()
	}
	return
}

func (rotkg *RotKGProtocol) AllocateShareRotColRight() (r RotKGShareRotColRight) {
	r = make([]*ring.Poly, rotkg.bfvContext.Beta())
	for i := uint64(0); i < rotkg.bfvContext.Beta(); i++ {
		r[i] = rotkg.bfvContext.ContextKeys().NewPoly()
	}
	return
}

func (rotkg *RotKGProtocol) AllocateShareConjugate() (r RotKGShareConjugate) {
	r = make([]*ring.Poly, rotkg.bfvContext.Beta())
	for i := uint64(0); i < rotkg.bfvContext.Beta(); i++ {
		r[i] = rotkg.bfvContext.ContextKeys().NewPoly()
	}
	return
}

// newrotkg creates a new rotkg object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRotKGProtocol(bfvContext *bfv.BfvContext) (rotkg *RotKGProtocol) {

	rotkg = new(RotKGProtocol)
	rotkg.bfvContext = bfvContext

	rotkg.gaussianSampler = bfvContext.GaussianSampler()

	rotkg.rot_col_L = make(map[uint64][][2]*ring.Poly)
	rotkg.rot_col_R = make(map[uint64][][2]*ring.Poly)

	rotkg.polypool = bfvContext.ContextKeys().NewPoly()

	N := bfvContext.ContextKeys().N

	rotkg.gen = 5
	rotkg.genInv = ring.ModExp(rotkg.gen, (N<<1)-1, N<<1)

	mask := (N << 1) - 1

	rotkg.galElRotColLeft = make([]uint64, N>>1)
	rotkg.galElRotColRight = make([]uint64, N>>1)

	rotkg.galElRotColRight[0] = 1
	rotkg.galElRotColLeft[0] = 1

	for i := uint64(1); i < N>>1; i++ {
		rotkg.galElRotColLeft[i] = (rotkg.galElRotColLeft[i-1] * rotkg.gen) & mask
		rotkg.galElRotColRight[i] = (rotkg.galElRotColRight[i-1] * rotkg.genInv) & mask
	}

	rotkg.galElRotRow = (N << 1) - 1

	return rotkg
}

// GenShareRotLeft is the first and unique round of the rotkg protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + (pi(s_i) - s_i) + e]
//
// and broadcasts it to the other j-1 parties. The protocol must be repeated for each desired rotation.
func (rotkg *RotKGProtocol) GenShareRotLeft(sk *ring.Poly, k uint64, crp []*ring.Poly, shareOut RotKGShareRotColLeft) {
	rotkg.genShare(sk, rotkg.galElRotColLeft[k&((rotkg.bfvContext.N()>>1)-1)], crp, shareOut)
}

// GenShareRotLeft is the first and unique round of the rotkg protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + (pi(s_i) - s_i) + e]
//
// and broadcasts it to the other j-1 parties. The protocol must be repeated for each desired rotation.
func (rotkg *RotKGProtocol) GenShareRotRight(sk *ring.Poly, k uint64, crp []*ring.Poly, shareOut RotKGShareRotColRight) {
	rotkg.genShare(sk, rotkg.galElRotColRight[k&((rotkg.bfvContext.N()>>1)-1)], crp, shareOut)
}

// GenShareRotLeft is the first and unique round of the rotkg protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + (pi(s_i) - s_i) + e_i]
//
// and broadcasts it to the other j-1 parties.
func (rotkg *RotKGProtocol) GenShareConjugate(sk *ring.Poly, crp []*ring.Poly, shareOut RotKGShareConjugate) {
	rotkg.genShare(sk, rotkg.galElRotRow, crp, shareOut)
}

// genswitchkey is a generic method to generate the public-share of the collective rotation-key.
func (rotkg *RotKGProtocol) genShare(sk *ring.Poly, galEl uint64, crp []*ring.Poly, evakey []*ring.Poly) {

	contextKeys := rotkg.bfvContext.ContextKeys()

	ring.PermuteNTT(sk, galEl, rotkg.polypool)
	contextKeys.Sub(rotkg.polypool, sk, rotkg.polypool)

	for _, pj := range rotkg.bfvContext.KeySwitchPrimes() {
		contextKeys.MulScalar(rotkg.polypool, pj, rotkg.polypool)
	}

	contextKeys.InvMForm(rotkg.polypool, rotkg.polypool)

	var index uint64

	for i := uint64(0); i < rotkg.bfvContext.Beta(); i++ {

		// e
		evakey[i] = rotkg.gaussianSampler.SampleNTTNew()

		// a is the CRP

		// e + sk_in * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0
		for j := uint64(0); j < rotkg.bfvContext.Alpha(); j++ {

			index = i*rotkg.bfvContext.Alpha() + j

			for w := uint64(0); w < contextKeys.N; w++ {
				evakey[i].Coeffs[index][w] = ring.CRed(evakey[i].Coeffs[index][w]+rotkg.polypool.Coeffs[index][w], contextKeys.Modulus[index])
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= uint64(len(rotkg.bfvContext.ContextQ().Modulus)-1) {
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
func (rotkg *RotKGProtocol) Aggregate(share1, share2, shareOut []*ring.Poly) {
	contextKeys := rotkg.bfvContext.ContextKeys()

	for i := uint64(0); i < rotkg.bfvContext.Beta(); i++ {
		contextKeys.Add(share1[i], share2[i], shareOut[i])
	}
}

func (rotkg *RotKGProtocol) StoreRotColLeft(share RotKGShareRotColLeft, k uint64, crp []*ring.Poly) {

	k &= ((rotkg.bfvContext.N() >> 1) - 1)

	rotkg.rot_col_L[k] = make([][2]*ring.Poly, rotkg.bfvContext.Beta())

	for i := uint64(0); i < rotkg.bfvContext.Beta(); i++ {
		rotkg.rot_col_L[k][i][0] = share[i].CopyNew()
		rotkg.rot_col_L[k][i][1] = crp[i].CopyNew()
		rotkg.bfvContext.ContextKeys().MForm(rotkg.rot_col_L[k][i][1], rotkg.rot_col_L[k][i][1])
	}
}

func (rotkg *RotKGProtocol) StoreRotColRight(share RotKGShareRotColRight, k uint64, crp []*ring.Poly) {

	k &= ((rotkg.bfvContext.N() >> 1) - 1)

	rotkg.rot_col_R[k] = make([][2]*ring.Poly, rotkg.bfvContext.Beta())

	for i := uint64(0); i < rotkg.bfvContext.Beta(); i++ {
		rotkg.rot_col_R[k][i][0] = share[i].CopyNew()
		rotkg.rot_col_R[k][i][1] = crp[i].CopyNew()
		rotkg.bfvContext.ContextKeys().MForm(rotkg.rot_col_R[k][i][1], rotkg.rot_col_R[k][i][1])
	}
}

func (rotkg *RotKGProtocol) StoreConjugate(share RotKGShareConjugate, crp []*ring.Poly) {

	rotkg.rot_row = make([][2]*ring.Poly, rotkg.bfvContext.Beta())

	for i := uint64(0); i < rotkg.bfvContext.Beta(); i++ {
		rotkg.rot_row[i][0] = share[i].CopyNew()
		rotkg.rot_row[i][1] = crp[i].CopyNew()
		rotkg.bfvContext.ContextKeys().MForm(rotkg.rot_row[i][1], rotkg.rot_row[i][1])
	}
}

// Finalize retrieves all the aggregated rotation-key, creates a new RotationKeys structur,
// fills it with the collective rotation keys and returns it.
func (rotkg *RotKGProtocol) Finalize() (rotkey *bfv.RotationKeys) {
	rotkey = rotkg.bfvContext.NewRotationKeysEmpty()

	for k := range rotkg.rot_col_L {
		rotkey.SetRotColLeft(rotkg.rot_col_L[k], k)
		delete(rotkg.rot_col_L, k)
	}

	for k := range rotkg.rot_col_R {
		rotkey.SetRotColRight(rotkg.rot_col_R[k], k)
		delete(rotkg.rot_col_R, k)
	}

	if rotkg.rot_row != nil {
		rotkey.SetRotRow(rotkg.rot_row)
		rotkg.rot_row = nil
	}

	return rotkey
}
