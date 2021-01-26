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

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey *SwitchingKey
}

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	permuteNTTIndex map[uint64][]uint64
	keys            map[uint64]*SwitchingKey
	params          *Parameters
}

// BootstrappingKey is a structure that stores the switching-keys required during the bootstrapping.
type BootstrappingKey struct {
	relinkey *EvaluationKey // Relinearization key
	rotkeys  *RotationKeys  // Rotation and conjugation keys
}

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
func NewRelinKey(params *Parameters) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)
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
func (evk *EvaluationKey) Get() *SwitchingKey {
	return evk.evakey
}

// Set sets the target Evaluation key with the input polynomials.
func (evk *EvaluationKey) Set(rlk [][2]*ring.Poly) {

	evk.evakey = new(SwitchingKey)
	evk.evakey.key = make([][2]*ring.Poly, len(rlk))
	for j := range rlk {
		evk.evakey.key[j][0] = rlk[j][0].CopyNew()
		evk.evakey.key[j][1] = rlk[j][1].CopyNew()
	}
}

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

// NewRotationKeys returns a new empty RotationKeys struct.
func NewRotationKeys(params *Parameters) (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	rotKey.keys = make(map[uint64]*SwitchingKey, 0)
	rotKey.permuteNTTIndex = make(map[uint64][]uint64)
	rotKey.params = params.Copy()
	return
}

// Delete empties the set of rotation keys
func (rtk RotationKeys) Delete() {
	for k := range rtk.keys {
		delete(rtk.keys, k)
	}
}

// SetRotKeyGalEl copies the given switching key in the set
func (rotKeys *RotationKeys) SetRotKeyGalEl(galoisEl uint64, swk *SwitchingKey) {
	rotKey, inSet := rotKeys.keys[galoisEl]
	if !inSet {
		rotKey = new(SwitchingKey)
		rotKeys.keys[galoisEl] = rotKey
	}
	if rotKey != swk {
		rotKey.Copy(swk)
	}
}

func (rotKeys *RotationKeys) GetRotKey(galoisEl uint64) (*SwitchingKey, bool) {
	rotKey, inSet := rotKeys.keys[galoisEl]
	return rotKey, inSet
}

// // SetRotKey sets the target RotationKeys' SwitchingKey for the specified rotation type and amount with the input polynomials.
// func (rotKey *RotationKeys) SetRotKey(params *Parameters, evakey [][2]*ring.Poly, rotType Rotation, k uint64) {

// 	switch rotType {
// 	case RotationLeft:

// 		if rotKey.evakeyRotColLeft == nil {
// 			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
// 		}

// 		if rotKey.permuteNTTLeftIndex == nil {
// 			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
// 		}

// 		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {

// 			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(GaloisGen, k, params.N())

// 			rotKey.evakeyRotColLeft[k] = new(SwitchingKey)
// 			rotKey.evakeyRotColLeft[k].evakey = make([][2]*ring.Poly, len(evakey))
// 			for j := range evakey {
// 				rotKey.evakeyRotColLeft[k].evakey[j][0] = evakey[j][0].CopyNew()
// 				rotKey.evakeyRotColLeft[k].evakey[j][1] = evakey[j][1].CopyNew()
// 			}
// 		}

// 	case RotationRight:

// 		if rotKey.evakeyRotColRight == nil {
// 			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
// 		}

// 		if rotKey.permuteNTTRightIndex == nil {
// 			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
// 		}

// 		if rotKey.evakeyRotColRight[k] == nil && k != 0 {

// 			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(GaloisGen, 2*params.N()-1-k, params.N())

// 			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
// 			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
// 			for j := range evakey {
// 				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
// 				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
// 			}
// 		}

// 	case Conjugate:

// 		if rotKey.evakeyConjugate == nil {

// 			rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(2*params.N()-1, 1, params.N())

// 			rotKey.evakeyConjugate = new(SwitchingKey)
// 			rotKey.evakeyConjugate.evakey = make([][2]*ring.Poly, len(evakey))
// 			for j := range evakey {
// 				rotKey.evakeyConjugate.evakey[j][0] = evakey[j][0].CopyNew()
// 				rotKey.evakeyConjugate.evakey[j][1] = evakey[j][1].CopyNew()
// 			}
// 		}
// 	}
// }
