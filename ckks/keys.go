package ckks

import "github.com/ldsec/lattigo/v2/ring"

// SecretKey is a structure that stores the SecretKey
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the PublicKey
type PublicKey struct {
	pk [2]*ring.Poly
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	key [][2]*ring.Poly
}

// RelinearizationKey is a structure that stores the switching-keys required during the relinearization.
type RelinearizationKey struct {
	evakey *SwitchingKey
}

// RotationKeySet is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeySet struct {
	permuteNTTIndex map[uint64][]uint64
	keys            map[uint64]*SwitchingKey
}

// EvaluationKey is a structure representing the complete set of switching-keys required for evaluation
type EvaluationKey struct {
	Rlk  *RelinearizationKey
	Rtks *RotationKeySet
}

// BootstrappingKey is a structure that stores the switching-keys required during the bootstrapping.
type BootstrappingKey EvaluationKey

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params *Parameters) *SecretKey {

	sk := new(SecretKey)
	sk.sk = params.NewPolyQP()
	return sk
}

// Get returns the value of the SecretKey.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the value of the SecretKey to the provided value.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {

	pk = new(PublicKey)

	pk.pk[0] = params.NewPolyQP()
	pk.pk[1] = params.NewPolyQP()

	return
}

// Get returns the value of the the public key.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the value of the public key to the provided value.
func (pk *PublicKey) Set(poly [2]*ring.Poly) {
	pk.pk[0] = poly[0].CopyNew()
	pk.pk[1] = poly[1].CopyNew()
}

// Get returns the switching key backing slice
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.key
}

// NewRelinKey returns a new EvaluationKey with zero values.
func NewRelinKey(params *Parameters) (evakey *RelinearizationKey) {

	evakey = new(RelinearizationKey)
	evakey.evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey.key = make([][2]*ring.Poly, params.Beta())
	for i := uint64(0); i < params.Beta(); i++ {

		evakey.evakey.key[i][0] = params.NewPolyQP()
		evakey.evakey.key[i][1] = params.NewPolyQP()
	}

	return
}

// Get returns the slice of switching keys of the evaluation-key.
func (evk *RelinearizationKey) Get() *SwitchingKey {
	return evk.evakey
}

// Set sets the target Evaluation key with the input polynomials.
func (evk *RelinearizationKey) Set(rlk [][2]*ring.Poly) {

	evk.evakey = new(SwitchingKey)
	evk.evakey.key = make([][2]*ring.Poly, len(rlk))
	for j := range rlk {
		evk.evakey.key[j][0] = rlk[j][0].CopyNew()
		evk.evakey.key[j][1] = rlk[j][1].CopyNew()
	}
}

// Copy copies other into the receiver.
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

// NewRotationKeySet returns a new empty RotationKeys struct.
func NewRotationKeySet(params *Parameters) (rotKey *RotationKeySet) {
	rotKey = new(RotationKeySet)
	rotKey.keys = make(map[uint64]*SwitchingKey, 0)
	rotKey.permuteNTTIndex = make(map[uint64][]uint64)
	return
}

// GetRotationKey return the rotation key for the given galois element or nil if such key is not in the set. The
// second argument is true  iff the first one is non-nil.
func (rtks *RotationKeySet) GetRotationKey(galoisEl uint64) (*SwitchingKey, bool) {
	rotKey, inSet := rtks.keys[galoisEl]
	return rotKey, inSet
}

// Set stores a copy of the rotKey SwitchingKey inside the receiver RotationKeySet
func (rtks *RotationKeySet) Set(galoisEl uint64, swk *SwitchingKey) {
	s := new(SwitchingKey)
	s.Copy(swk)
	rtks.keys[galoisEl] = s
	N := uint64(len(swk.key[0][0].Coeffs[0]))
	rtks.permuteNTTIndex[galoisEl] = ring.PermuteNTTIndex(galoisEl, N)
}

// delete empties the set of rotation keys
func (rtks RotationKeySet) delete() {
	for k := range rtks.keys {
		delete(rtks.keys, k)
		delete(rtks.permuteNTTIndex, k)
	}
}
