package mkckks

import (
	"encoding/binary"
	"hash/fnv"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MKParticipant is a type for participants in a multi key ckks scheme
type MKParticipant interface {
	GetID() uint64
	GetEvaluationKey() *mkrlwe.MKEvaluationKey
	GetPublicKey() *mkrlwe.MKPublicKey
	Encrypt(values []complex128) *mkrlwe.MKCiphertext
	Decrypt(cipher *mkrlwe.MKCiphertext, partialDecryptions []*ring.Poly) []complex128
	GetPartialDecryption(ciphertext *mkrlwe.MKCiphertext) *ring.Poly
	GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey
	GetSecretKey() *mkrlwe.MKSecretKey
	SetSecretKey(newKey *mkrlwe.MKSecretKey)
}

type mkParticipant struct {
	id        uint64
	encryptor MKEncryptor
	decryptor mkrlwe.MKDecryptor
	keys      *mkrlwe.MKKeys
	encoder   ckks.Encoder
	ringQ     *ring.Ring
	params    *rlwe.Parameters
}

// GetID returns the id of the participant
func (participant *mkParticipant) GetID() uint64 {
	return participant.id
}

// GetEvaluationKey returns the evaluation key of the participant
func (participant *mkParticipant) GetEvaluationKey() *mkrlwe.MKEvaluationKey {
	return participant.keys.EvalKey
}

// GetPublicKey returns the publik key of the participant
func (participant *mkParticipant) GetPublicKey() *mkrlwe.MKPublicKey {
	return participant.keys.PublicKey
}

// Encrypt constructs a ciphertext from the given values
func (participant *mkParticipant) Encrypt(values []complex128) *mkrlwe.MKCiphertext {
	if values == nil || len(values) <= 0 {
		panic("Cannot encrypt uninitialized or empty values")
	}
	return participant.encryptor.EncryptMK(participant.encoder.EncodeNTTAtLvlNew(participant.params.MaxLevel(), values, participant.params.LogSlots()))
}

// Decrypt returns the decryption of the ciphertext given the partial decryption
func (participant *mkParticipant) Decrypt(cipher *mkrlwe.MKCiphertext, partialDecryptions []*ring.Poly) []complex128 {

	if cipher == nil || cipher.Ciphertexts == nil || len(cipher.Ciphertexts.Value) < 2 {
		panic("Cannot decrypt uninitialized ciphertext nor ciphertext containing only one value")
	}

	if partialDecryptions == nil || len(partialDecryptions) < 1 {
		panic("Decryption necessitates at least one partialy decrypted ciphertext")
	}

	decrypted := participant.decryptor.MergeDec(cipher.Ciphertexts.Value[0], cipher.Ciphertexts.Scale(), cipher.Ciphertexts.Level(), partialDecryptions)

	return participant.encoder.Decode(decrypted, participant.params.LogSlots())
}

// GetPartialDecryption returns the partial decryption of an element in the ciphertext
// this function should only be used by participants that were involved in the given ciphertext
func (participant *mkParticipant) GetPartialDecryption(ciphertext *mkrlwe.MKCiphertext) *ring.Poly {

	cipherPart := participant.getCiphertextPart(ciphertext)

	if cipherPart == nil {
		panic("Participant is not involved in the given ciphertext. Partial decryption impossible.")
	}
	return participant.decryptor.PartDec(cipherPart, ciphertext.Ciphertexts.Level(), participant.keys.SecretKey)
}

// NewParticipant creates a participant for the multi key ckks scheme
// the ckks parameters as well as the standard deviation used for partial decryption must be provided
func NewParticipant(params *rlwe.Parameters, sigmaSmudging float64, crs *mkrlwe.MKDecomposedPoly) MKParticipant {

	if crs == nil || params == nil {
		panic("Uninitialized parameters. Cannot create new participant")
	}

	if sigmaSmudging < params.Sigma() {
		panic("Sigma must be at least greater than the standard deviation of the gaussian distribution")
	}

	if len(crs.Poly) != int(params.Beta()) {
		panic("CRS must be the same dimention as returned by the function ckks.Parameters.Beta()")
	}

	keys := mkrlwe.KeyGen(params, mkrlwe.CopyNewDecomposed(crs))

	uid := hashPublicKey(keys.PublicKey.Key)

	keys.PublicKey.PeerID = uid
	keys.SecretKey.PeerID = uid
	keys.EvalKey.PeerID = uid

	encryptor := NewMKEncryptor(keys.PublicKey, params, uid)
	decryptor := mkrlwe.NewMKDecryptor(params, sigmaSmudging)
	encoder := ckks.NewEncoder(params)
	ringQ := mkrlwe.GetRingQ(params)

	return &mkParticipant{
		id:        uid,
		encryptor: encryptor,
		decryptor: decryptor,
		keys:      keys,
		encoder:   encoder,
		ringQ:     ringQ,
		params:    params,
	}
}

// computes the hash of the public key using the FNV hashing algorithm
func hashPublicKey(pk [2]*mkrlwe.MKDecomposedPoly) uint64 {

	coeffs := pk[0].Poly[0].Coeffs // b[0] is the ckks public key
	h64 := fnv.New64()

	for _, v := range coeffs {
		for _, vIn := range v {
			h64.Write(uintToBytes(vIn))
		}
	}

	return h64.Sum64()
}

// converts a uint64 in a slice of bytes
func uintToBytes(i uint64) []byte {

	b := make([]byte, 8)
	binary.LittleEndian.PutUint64(b, uint64(i))

	return b
}

// returns the part of the ciphertext corresponding to the participant
func (participant *mkParticipant) getCiphertextPart(ciphertext *mkrlwe.MKCiphertext) *ring.Poly {

	for i, v := range ciphertext.PeerIDs {

		if v == participant.id {
			return ciphertext.Ciphertexts.Value[i+1]
		}
	}

	return nil

}

// GetRotationKeys returns the rotation key set associated with the given rotation
func (participant *mkParticipant) GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey {

	galEl := participant.params.GaloisElementForColumnRotationBy(rot)

	evalKey := mkrlwe.GaloisEvaluationKeyGen(galEl, participant.keys.SecretKey, participant.params)
	evalKey.PeerID = participant.id

	return evalKey
}

// GetSecretKey returns the secret key of the participant
func (participant *mkParticipant) GetSecretKey() *mkrlwe.MKSecretKey { // TODO: remove these 2 functions and the key switch if not useful else run Keygen with new key
	return participant.keys.SecretKey
}

// SetSecretKey changes the secret key of the participant
func (participant *mkParticipant) SetSecretKey(newKey *mkrlwe.MKSecretKey) {
	participant.keys.SecretKey = newKey
}
