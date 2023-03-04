package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// SecretKey is a type for generic RLWE secret keys.
// The Value field stores the polynomial in NTT and Montgomery form.
type SecretKey struct {
	Value ringqp.Poly
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params Parameters) *SecretKey {
	return &SecretKey{Value: params.RingQP().NewPoly()}
}

// LevelQ returns the level of the modulus Q of the target.
func (sk *SecretKey) LevelQ() int {
	return sk.Value.Q.Level()
}

// LevelP returns the level of the modulus P of the target.
// Returns -1 if P is absent.
func (sk *SecretKey) LevelP() int {
	if sk.Value.P != nil {
		return sk.Value.P.Level()
	}

	return -1
}

// CopyNew creates a deep copy of the receiver secret key and returns it.
func (sk *SecretKey) CopyNew() *SecretKey {
	if sk == nil {
		return nil
	}
	return &SecretKey{sk.Value.CopyNew()}
}

// PublicKey is a type for generic RLWE public keys.
// The Value field stores the polynomials in NTT and Montgomery form.
type PublicKey struct {
	CiphertextQP
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params Parameters) (pk *PublicKey) {
	return &PublicKey{CiphertextQP{Value: [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()}, MetaData: MetaData{IsNTT: true, IsMontgomery: true}}}
}

// LevelQ returns the level of the modulus Q of the target.
func (pk *PublicKey) LevelQ() int {
	return pk.Value[0].Q.Level()
}

// LevelP returns the level of the modulus P of the target.
// Returns -1 if P is absent.
func (pk *PublicKey) LevelP() int {
	if pk.Value[0].P != nil {
		return pk.Value[0].P.Level()
	}

	return -1
}

// Equals checks two PublicKey struct for equality.
func (pk *PublicKey) Equals(other *PublicKey) bool {
	if pk == other {
		return true
	}
	return pk.Value[0].Equals(other.Value[0]) && pk.Value[1].Equals(other.Value[1])
}

// CopyNew creates a deep copy of the receiver PublicKey and returns it.
func (pk *PublicKey) CopyNew() *PublicKey {
	if pk == nil {
		return nil
	}
	return &PublicKey{*pk.CiphertextQP.CopyNew()}
}

// EvaluationKeySetInterface is an interface implementing methods
// to load the RelinearizationKey and GaloisKeys in the Evaluator.
type EvaluationKeySetInterface interface {
	Add(evk interface{}) (err error)
	GetGaloisKey(galEl uint64) (evk *GaloisKey, err error)
	GetGaloisKeysList() (galEls []uint64)
	GetRelinearizationKey() (evk *RelinearizationKey, err error)
}

// EvaluationKeySet is a generic struct that complies to the `EvaluationKeys` interface.
// This interface can be re-implemented by users to suit application specific requirement,
// notably evaluation keys loading and persistence.
type EvaluationKeySet struct {
	*RelinearizationKey
	GaloisKeys map[uint64]*GaloisKey
}

// Add stores the evaluation key in the EvaluationKeySet.
// Supported types are *rlwe.EvalutionKey and *rlwe.GaloiKey
func (evk *EvaluationKeySet) Add(key interface{}) (err error) {
	switch key := key.(type) {
	case *RelinearizationKey:
		evk.RelinearizationKey = key
	case *GaloisKey:
		evk.GaloisKeys[key.GaloisElement] = key
	default:
		return fmt.Errorf("unsuported type. Supported types are *rlwe.EvalutionKey and *rlwe.GaloiKey, but have %T", key)
	}

	return
}

// NewEvaluationKeySet returns a new EvaluationKeySet with nil RelinearizationKey and empty GaloisKeys map.
func NewEvaluationKeySet() (evk *EvaluationKeySet) {
	return &EvaluationKeySet{
		RelinearizationKey: nil,
		GaloisKeys:         make(map[uint64]*GaloisKey),
	}
}

func (evk *EvaluationKeySet) GetGaloisKey(galEl uint64) (gk *GaloisKey, err error) {
	var ok bool
	if gk, ok = evk.GaloisKeys[galEl]; !ok {
		return nil, fmt.Errorf("GaloiKey[%d] is nil", galEl)
	}

	return
}

func (evk *EvaluationKeySet) GetGaloisKeysList() (galEls []uint64) {

	if evk.GaloisKeys == nil {
		return []uint64{}
	}

	galEls = make([]uint64, len(evk.GaloisKeys))

	var i int
	for galEl := range evk.GaloisKeys {
		galEls[i] = galEl
		i++
	}

	return
}

func (evk *EvaluationKeySet) GetRelinearizationKey() (rk *RelinearizationKey, err error) {
	if evk.RelinearizationKey != nil {
		return evk.RelinearizationKey, nil
	}

	return nil, fmt.Errorf("RelinearizationKey is nil")
}

// EvaluationKey is a public key indended to be used during the evaluation phase of a homomorphic circuit.
// It provides a one way public and non-interactive re-encryption from a ciphertext encrypted under `skIn`
// to a ciphertext encrypted under `skOut`.
//
// Such re-encryption is for example used for:
//
// - Homomorphic relinearization: re-encryption of a quadratic ciphertext (that requires (1, sk sk^2) to be decrypted)
// to a linear ciphertext (that required (1, sk) to be decrypted). In this case skIn = sk^2 an skOut = sk.
//
// - Homomorphic automorphisms: an automorphism in the ring Z[X]/(X^{N}+1) is defined as pi_k: X^{i} -> X^{i^k} with
// k coprime to 2N. Pi_sk is for exampled used during homomorphic slot rotations. Applying pi_k to a ciphertext encrypted
// under sk generates a new ciphertext encrypted under pi_k(sk), and an Evaluationkey skIn = pi_k(sk) to skOut = sk
// is used to bring it back to its original key.
type EvaluationKey struct {
	GadgetCiphertext
}

// NewEvaluationKey returns a new EvaluationKey with pre-allocated zero-value
func NewEvaluationKey(params Parameters, levelQ, levelP int) *EvaluationKey {
	return &EvaluationKey{GadgetCiphertext: *NewGadgetCiphertext(
		params,
		levelQ,
		levelP,
		params.DecompRNS(levelQ, levelP),
		params.DecompPw2(levelQ, levelP),
	)}
}

// Equals checks two EvaluationKeys for equality.
func (evk *EvaluationKey) Equals(other *EvaluationKey) bool {
	return evk.GadgetCiphertext.Equals(&other.GadgetCiphertext)
}

// CopyNew creates a deep copy of the target EvaluationKey and returns it.
func (evk *EvaluationKey) CopyNew() *EvaluationKey {
	return &EvaluationKey{GadgetCiphertext: *evk.GadgetCiphertext.CopyNew()}
}

type RelinearizationKey struct {
	*EvaluationKey
}

func NewRelinearizationKey(params Parameters) *RelinearizationKey {
	return &RelinearizationKey{EvaluationKey: NewEvaluationKey(params, params.MaxLevelQ(), params.MaxLevelP())}
}

func (rlk *RelinearizationKey) Equals(other *RelinearizationKey) bool {
	return rlk.EvaluationKey.Equals(other.EvaluationKey)
}

func (rlk *RelinearizationKey) CopyNew() *RelinearizationKey {
	return &RelinearizationKey{EvaluationKey: rlk.EvaluationKey.CopyNew()}
}

type GaloisKey struct {
	GaloisElement uint64
	NthRoot       uint64
	*EvaluationKey
}

func NewGaloisKey(params Parameters) *GaloisKey {
	return &GaloisKey{EvaluationKey: NewEvaluationKey(params, params.MaxLevelQ(), params.MaxLevelP())}
}

func (gk *GaloisKey) Equals(other *GaloisKey) bool {
	return gk.EvaluationKey.Equals(other.EvaluationKey) && gk.GaloisElement == other.GaloisElement && gk.NthRoot == other.NthRoot
}

func (gk *GaloisKey) CopyNew() *GaloisKey {
	return &GaloisKey{
		GaloisElement: gk.GaloisElement,
		NthRoot:       gk.NthRoot,
		EvaluationKey: gk.EvaluationKey.CopyNew(),
	}
}
