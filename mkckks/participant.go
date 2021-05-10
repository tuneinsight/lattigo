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
	Encrypt(values []complex128) *MKCiphertext
	Decrypt(cipher *MKCiphertext, partialDecryptions []*ring.Poly) []complex128
	GetPartialDecryption(ciphertext *MKCiphertext) *ring.Poly
	GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey
	GetSecretKey() *mkrlwe.MKSecretKey
}

type mkParticipant struct {
	id        uint64
	encryptor MKEncryptor
	decryptor mkrlwe.MKDecryptor
	keys      *mkrlwe.MKKeys
	encoder   ckks.Encoder
	ringQ     *ring.Ring
	params    *ckks.Parameters
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
func (participant *mkParticipant) Encrypt(values []complex128) *MKCiphertext {
	if values == nil || len(values) <= 0 {
		panic("Cannot encrypt uninitialized or empty values")
	}
	return participant.encryptor.EncryptMK(participant.encoder.EncodeNTTAtLvlNew(participant.params.MaxLevel(), values, participant.params.LogSlots()))
}

// Decrypt returns the decryption of the ciphertext given the partial decryption
func (participant *mkParticipant) Decrypt(cipher *MKCiphertext, partialDecryptions []*ring.Poly) []complex128 {

	if cipher == nil || cipher.Ciphertexts == nil || len(cipher.Ciphertexts.Value) < 2 {
		panic("Cannot decrypt uninitialized ciphertext nor ciphertext containing only one value")
	}

	if partialDecryptions == nil || len(partialDecryptions) < 1 {
		panic("Decryption necessitates at least one partialy decrypted ciphertext")
	}

	decrypted := participant.decryptor.MergeDec(cipher.Ciphertexts.Value[0], cipher.Ciphertexts.Level(), partialDecryptions)

	pt := ckks.NewPlaintext(*participant.params, cipher.Ciphertexts.Level(), cipher.Ciphertexts.Scale())

	pt.SetValue(decrypted)

	return participant.encoder.Decode(pt, participant.params.LogSlots())
}

// GetPartialDecryption returns the partial decryption of an element in the ciphertext
// this function should only be used by participants that were involved in the given ciphertext
func (participant *mkParticipant) GetPartialDecryption(ciphertext *MKCiphertext) *ring.Poly {

	cipherPart := participant.getCiphertextPart(ciphertext)

	if cipherPart == nil {
		panic("Participant is not involved in the given ciphertext. Partial decryption impossible.")
	}
	return participant.decryptor.PartDec(cipherPart, ciphertext.Ciphertexts.Level(), participant.keys.SecretKey)
}

// NewParticipant creates a participant for the multi key ckks scheme
// the ckks parameters as well as the standard deviation used for partial decryption must be provided
func NewParticipant(params *ckks.Parameters, sigmaSmudging float64, crs *mkrlwe.MKDecomposedPoly) MKParticipant {

	if crs == nil || params == nil {
		panic("Uninitialized parameters. Cannot create new participant")
	}

	if sigmaSmudging < params.Sigma() {
		panic("Sigma must be at least greater than the standard deviation of the gaussian distribution")
	}

	if len(crs.Poly) != int(params.Beta()) {
		panic("CRS must be the same dimention as returned by the function ckks.Parameters.Beta()")
	}

	keys := mkrlwe.KeyGen(&params.Parameters, mkrlwe.CopyNewDecomposed(crs))

	uid := hashPublicKey(keys.PublicKey.Key)

	keys.PublicKey.PeerID = uid
	keys.SecretKey.PeerID = uid
	keys.EvalKey.PeerID = uid

	encryptor := NewMKEncryptor(keys.PublicKey, params, uid)
	decryptor := mkrlwe.NewMKDecryptor(&params.Parameters, sigmaSmudging)
	encoder := ckks.NewEncoder(*params)
	ringQ := mkrlwe.GetRingQ(&params.Parameters)

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

// NewParticipantFromSecretKey creates a participant for the multi key bfv scheme from a bfv secret key
// the bfv parameters as well as the standard deviation used for partial decryption must be provided
func NewParticipantFromSecretKey(params *ckks.Parameters, sigmaSmudging float64, crs *mkrlwe.MKDecomposedPoly, sk *rlwe.SecretKey) MKParticipant {

	if crs == nil || params == nil || sk == nil {
		panic("Uninitialized parameters. Cannot create new participant")
	}

	if sigmaSmudging < params.Sigma() {
		panic("Sigma must be at least greater than the standard deviation of the gaussian distribution")
	}

	if len(crs.Poly) != int(params.Beta()) {
		panic("CRS must be the same dimention as returned by the function bfv.Parameters.Beta()")
	}

	keys := mkrlwe.KeyGenWithSecretKey(&params.Parameters, mkrlwe.CopyNewDecomposed(crs), sk)

	uid := hashPublicKey(keys.PublicKey.Key)

	keys.PublicKey.PeerID = uid
	keys.SecretKey.PeerID = uid
	keys.EvalKey.PeerID = uid

	encryptor := NewMKEncryptor(keys.PublicKey, params, uid)
	decryptor := mkrlwe.NewMKDecryptor(&params.Parameters, sigmaSmudging)
	encoder := ckks.NewEncoder(*params)
	ringQ := mkrlwe.GetRingQ(&params.Parameters)

	return &mkParticipant{
		id:        uid,
		encryptor: encryptor,
		decryptor: decryptor,
		params:    params,
		keys:      keys,
		encoder:   encoder,
		ringQ:     ringQ,
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
func (participant *mkParticipant) getCiphertextPart(ciphertext *MKCiphertext) *ring.Poly {

	for i, v := range ciphertext.PeerID {

		if v == participant.id {
			return ciphertext.Ciphertexts.Value[i+1]
		}
	}

	return nil

}

// GetRotationKeys returns the rotation key set associated with the given rotation
func (participant *mkParticipant) GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey {

	galEl := participant.params.GaloisElementForColumnRotationBy(rot)

	evalKey := mkrlwe.GaloisEvaluationKeyGen(galEl, participant.keys.SecretKey, &participant.params.Parameters)
	evalKey.PeerID = participant.id

	return evalKey
}

// GetSecretKey returns the secret key of the participant
func (participant *mkParticipant) GetSecretKey() *mkrlwe.MKSecretKey {
	return participant.keys.SecretKey
}
