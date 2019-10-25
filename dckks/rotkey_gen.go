package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

// RKG is the structure storing the parameters for the collective rotation-keys generation.
type RKG struct {
	ckksContext *ckks.CkksContext

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

// newRKG creates a new RKG object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRKG(ckksContext *ckks.CkksContext) (rkg *RKG) {

	rkg = new(RKG)
	rkg.ckksContext = ckksContext

	rkg.gaussianSampler = ckksContext.GaussianSampler()

	rkg.rot_col_L = make(map[uint64][][2]*ring.Poly)
	rkg.rot_col_R = make(map[uint64][][2]*ring.Poly)

	rkg.polypool = ckksContext.ContextKeys().NewPoly()

	N := ckksContext.ContextKeys().N

	rkg.gen = 5
	rkg.genInv = ring.ModExp(rkg.gen, (N<<1)-1, N<<1)

	mask := (N << 1) - 1

	rkg.galElRotColLeft = make([]uint64, N>>1)
	rkg.galElRotColRight = make([]uint64, N>>1)

	rkg.galElRotColRight[0] = 1
	rkg.galElRotColLeft[0] = 1

	for i := uint64(1); i < N>>1; i++ {
		rkg.galElRotColLeft[i] = (rkg.galElRotColLeft[i-1] * rkg.gen) & mask
		rkg.galElRotColRight[i] = (rkg.galElRotColRight[i-1] * rkg.genInv) & mask

	}

	rkg.galElRotRow = (N << 1) - 1

	return rkg
}

// GenShareRotLeft is the first and unique round of the RKG protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + (pi(s_i) - s_i) + e]
//
// and broadcasts it to the other j-1 parties. The protocol must be repeated for each desired rotation.
func (rkg *RKG) GenShareRotLeft(sk *ring.Poly, k uint64, crp []*ring.Poly) (evakey []*ring.Poly) {

	k &= (rkg.ckksContext.ContextKeys().N >> 1) - 1

	ring.PermuteNTT(sk, rkg.galElRotColLeft[k], rkg.polypool)
	rkg.ckksContext.ContextKeys().Sub(rkg.polypool, sk, rkg.polypool)

	for _, pj := range rkg.ckksContext.KeySwitchPrimes() {
		rkg.ckksContext.ContextKeys().MulScalar(rkg.polypool, pj, rkg.polypool)
	}

	rkg.ckksContext.ContextKeys().InvMForm(rkg.polypool, rkg.polypool)

	return rkg.genswitchkey(rkg.polypool, sk, crp)
}

// GenShareRotLeft is the first and unique round of the RKG protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + (pi(s_i) - s_i) + e]
//
// and broadcasts it to the other j-1 parties. The protocol must be repeated for each desired rotation.
func (rkg *RKG) GenShareRotRight(sk *ring.Poly, k uint64, crp []*ring.Poly) (evakey []*ring.Poly) {

	k &= (rkg.ckksContext.ContextKeys().N >> 1) - 1

	ring.PermuteNTT(sk, rkg.galElRotColRight[k], rkg.polypool)
	rkg.ckksContext.ContextKeys().Sub(rkg.polypool, sk, rkg.polypool)

	for _, pj := range rkg.ckksContext.KeySwitchPrimes() {
		rkg.ckksContext.ContextKeys().MulScalar(rkg.polypool, pj, rkg.polypool)
	}

	rkg.ckksContext.ContextKeys().InvMForm(rkg.polypool, rkg.polypool)

	return rkg.genswitchkey(rkg.polypool, sk, crp)
}

// GenShareRotLeft is the first and unique round of the RKG protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + (pi(s_i) - s_i) + e_i]
//
// and broadcasts it to the other j-1 parties.
func (rkg *RKG) GenShareRotRow(sk *ring.Poly, crp []*ring.Poly) (evakey []*ring.Poly) {

	ring.PermuteNTT(sk, rkg.galElRotRow, rkg.polypool)
	rkg.ckksContext.ContextKeys().Sub(rkg.polypool, sk, rkg.polypool)

	for _, pj := range rkg.ckksContext.KeySwitchPrimes() {
		rkg.ckksContext.ContextKeys().MulScalar(rkg.polypool, pj, rkg.polypool)
	}

	rkg.ckksContext.ContextKeys().InvMForm(rkg.polypool, rkg.polypool)

	return rkg.genswitchkey(rkg.polypool, sk, crp)
}

// genswitchkey is a generic method to generate the public-share of the collective rotation-key.
func (rkg *RKG) genswitchkey(sk_in, sk_out *ring.Poly, crp []*ring.Poly) (evakey []*ring.Poly) {

	contextKeys := rkg.ckksContext.ContextKeys()

	evakey = make([]*ring.Poly, rkg.ckksContext.Beta())

	var index uint64

	for i := uint64(0); i < rkg.ckksContext.Beta(); i++ {

		// e
		evakey[i] = rkg.gaussianSampler.SampleNTTNew()

		// a is the CRP

		// e + sk_in * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0
		for j := uint64(0); j < rkg.ckksContext.Alpha(); j++ {

			index = i*rkg.ckksContext.Alpha() + j

			for w := uint64(0); w < contextKeys.N; w++ {
				evakey[i].Coeffs[index][w] = ring.CRed(evakey[i].Coeffs[index][w]+rkg.polypool.Coeffs[index][w], contextKeys.Modulus[index])
			}

			// Handles the case where nb pj does not divides nb qi
			if index == uint64(len(rkg.ckksContext.ContextQ().Modulus)) {
				break
			}
		}

		// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
		contextKeys.MulCoeffsMontgomeryAndSub(crp[i], sk_out, evakey[i])
		contextKeys.MForm(evakey[i], evakey[i])

	}

	return
}

// AggregateRotColL is the second part of the unique round of the RKG protocol. Uppon receiving the j-1 public shares,
// each party computes  :
//
// [sum(a*a_j + (pi(a_j) - a_j) + e_j), a]
func (rkg *RKG) AggregateRotColL(samples [][]*ring.Poly, k uint64, crp []*ring.Poly) {

	k &= (rkg.ckksContext.ContextKeys().N >> 1) - 1

	rkg.rot_col_L[k] = rkg.aggregate(samples, crp)
}

// AggregateRotColR is the second part of the unique round of the RKG protocol. Uppon receiving the j-1 public shares,
// each party computes  :
//
// [sum(a*a_j + (pi(a_j) - a_j) + e_j), a]
func (rkg *RKG) AggregateRotColR(samples [][]*ring.Poly, k uint64, crp []*ring.Poly) {

	k &= (rkg.ckksContext.ContextKeys().N >> 1) - 1

	rkg.rot_col_R[k] = rkg.aggregate(samples, crp)
}

// AggregateRotRow is the second part of the unique round of the RKG protocol. Uppon receiving the j-1 public shares,
// each party computes  :
//
// [sum(a*a_j + (pi(a_j) - a_j) + e_j), a]
func (rkg *RKG) AggregateRotRow(samples [][]*ring.Poly, crp []*ring.Poly) {

	rkg.rot_row = rkg.aggregate(samples, crp)
}

// aggregate is a generic methode for summing the samples and creating a structur similar to the rotation key from the
// summed public shares and the collective random polynomial.
func (rkg *RKG) aggregate(samples [][]*ring.Poly, crp []*ring.Poly) (receiver [][2]*ring.Poly) {

	contextKeys := rkg.ckksContext.ContextKeys()

	receiver = make([][2]*ring.Poly, rkg.ckksContext.Beta())

	for i := uint64(0); i < rkg.ckksContext.Beta(); i++ {

		receiver[i][0] = samples[0][i].CopyNew()
		receiver[i][1] = crp[i].CopyNew()
		contextKeys.MForm(receiver[i][1], receiver[i][1])

		for j := 1; j < len(samples); j++ {
			contextKeys.AddNoMod(receiver[i][0], samples[j][i], receiver[i][0])

			if j&7 == 7 {
				contextKeys.Reduce(receiver[i][0], receiver[i][0])
			}
		}

		if (len(samples)-1)&7 != 7 {
			contextKeys.Reduce(receiver[i][0], receiver[i][0])
		}

	}

	return
}

// Finalize retrieves all the aggregated rotation-key, creates a new RotationKeys structur,
// fills it with the collective rotation keys and returns it.
func (rkg *RKG) Finalize(ckksContext *ckks.CkksContext) (rotkey *ckks.RotationKey) {
	rotkey = ckksContext.NewRotationKeysEmpty()

	for k := range rkg.rot_col_L {
		rotkey.SetRotColLeft(rkg.rot_col_L[k], k)
		delete(rkg.rot_col_L, k)
	}

	for k := range rkg.rot_col_R {
		rotkey.SetRotColRight(rkg.rot_col_R[k], k)
		delete(rkg.rot_col_R, k)
	}

	if rkg.rot_row != nil {
		rotkey.SetRotRow(rkg.rot_row)
		rkg.rot_row = nil
	}

	return rotkey
}
