package bfv

import "github.com/ldsec/lattigo/v2/ring"

// SecretKey is a structure that stores the SecretKey.
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the PublicKey.
type PublicKey struct {
	pk [2]*ring.Poly
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	key [][2]*ring.Poly
}

// RelinearizationKey is a structure that stores the switching-keys required during the relinearization.
type RelinearizationKey struct {
	keys []*SwitchingKey
}

// RotationKeySet is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeySet struct {
	keys map[uint64]*SwitchingKey
}

// NewRotationKeySet returns a new empty RotationKeys struct.
func NewRotationKeySet(params *Parameters) (rotKey *RotationKeySet) {
	rotKey = new(RotationKeySet)
	rotKey.keys = make(map[uint64]*SwitchingKey, 0)
	return
}

// Set stores a copy of the rotKey SwitchingKey inside the receiver RotationKeySet
func (rtks *RotationKeySet) Set(galEl uint64, rotKey *SwitchingKey) {
	s := new(SwitchingKey)
	s.Copy(rotKey)
	rtks.keys[galEl] = s
}

// Delete empties the set of rotation keys
func (rtks *RotationKeySet) Delete() {
	for k := range rtks.keys {
		delete(rtks.keys, k)
	}
}

// GetRotationKey return the rotation key for the given galois element or nil if such key is not in the set. The
// second argument is true  iff the first one is non-nil.
func (rtks *RotationKeySet) GetRotationKey(galoisEl uint64) (*SwitchingKey, bool) {
	rotKey, inSet := rtks.keys[galoisEl]
	return rotKey, inSet
}

// Get returns the switching key backing slice.
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.key
}

// Copy copies the other SwitchingKey inside the receiver.
func (swk *SwitchingKey) Copy(other *SwitchingKey) {
	if other == nil {
		return
	}
	if len(swk.key) == 0 {
		swk.key = make([][2]*ring.Poly, len(other.key), len(other.key))
		for i, o := range other.key {
			n, q := uint64(o[0].GetDegree()), uint64(o[0].GetLenModuli())
			swk.key[i] = [2]*ring.Poly{ring.NewPoly(n, q), ring.NewPoly(n, q)}
		}
	}
	for i, o := range other.key {
		swk.key[i][0].Copy(o[0])
		swk.key[i][1].Copy(o[1])
	}
}

// NewRelinKey creates a new EvaluationKey with zero values.
func NewRelinKey(params *Parameters, maxDegree uint64) (evakey *RelinearizationKey) {

	evakey = new(RelinearizationKey)

	beta := params.Beta()

	evakey.keys = make([]*SwitchingKey, maxDegree)

	for w := uint64(0); w < maxDegree; w++ {

		evakey.keys[w] = new(SwitchingKey)

		evakey.keys[w].key = make([][2]*ring.Poly, beta)

		for i := uint64(0); i < beta; i++ {

			evakey.keys[w].key[i][0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
			evakey.keys[w].key[i][1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
		}
	}

	return
}

// Get returns the slice of SwitchingKeys of the target EvaluationKey.
func (evk *RelinearizationKey) Get() []*SwitchingKey {
	return evk.keys
}

// Set sets the polynomial of the target EvaluationKey as the input polynomials.
func (evk *RelinearizationKey) Set(rlk [][][2]*ring.Poly) {

	evk.keys = make([]*SwitchingKey, len(rlk))
	for i := range rlk {
		evk.keys[i] = new(SwitchingKey)
		evk.keys[i].key = make([][2]*ring.Poly, len(rlk[i]))
		for j := range rlk[i] {
			evk.keys[i].key[j][0] = rlk[i][j][0].CopyNew()
			evk.keys[i].key[j][1] = rlk[i][j][1].CopyNew()
		}
	}
}
