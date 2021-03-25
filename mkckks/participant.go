package mkckks

import (
	"encoding/binary"
	"hash/fnv"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKParticipant is a type for participants in a multi key ckks scheme
type MKParticipant interface {
	GetID() uint64
	Encrypt(values []complex128) *MKCiphertext
	Decrypt(cipher *MKCiphertext, partialDecryptions []*ring.Poly) []complex128
	GetPartialDecryption(ciphertext *MKCiphertext) *ring.Poly
}

type mkParticipant struct {
	id        uint64
	encryptor MKEncryptor
	decryptor MKDecryptor
	keys      *MKKeys
	encoder   ckks.Encoder
	ringQ     *ring.Ring
	params    *ckks.Parameters
}

// GetID returns the id of the participant
func (participant *mkParticipant) GetID() uint64 {
	return participant.id
}

// Encrypt constructs a ciphertext from the given values
func (participant *mkParticipant) Encrypt(values []complex128) *MKCiphertext {

	return participant.encryptor.EncryptMK(participant.encoder.EncodeNTTAtLvlNew(participant.params.MaxLevel(), values, participant.params.LogSlots()))
}

// Decrypt returns the decryption of the ciphertext given the partial decryption
func (participant *mkParticipant) Decrypt(cipher *MKCiphertext, partialDecryptions []*ring.Poly) []complex128 {

	decrypted := participant.decryptor.MergeDec(cipher.ciphertexts.Value()[0], cipher.ciphertexts.Scale(), cipher.ciphertexts.Level(), partialDecryptions)

	return participant.encoder.Decode(decrypted, participant.params.LogSlots())
}

// GetPartialDecryption returns the partial decryption of an element in the ciphertext
// this function should only be used by participants that were involved in the given ciphertext
func (participant *mkParticipant) GetPartialDecryption(ciphertext *MKCiphertext) *ring.Poly {

	cipherPart := participant.getCiphertextPart(ciphertext)

	if cipherPart == nil {
		panic("Participant is not involved in the given ciphertext. Partial decryption impossible.")
	}
	return participant.decryptor.PartDec(cipherPart, ciphertext.ciphertexts.Level(), participant.keys.secretKey)
}

// NewParticipant creates a participant for the multi key ckks scheme
// the ckks parameters as well as the standard deviation used for partial decryption must be provided
func NewParticipant(params *ckks.Parameters, sigmaSmudging float64) MKParticipant {

	a := GenCommonPublicParam(params) // TODO: check if it should really be common or not in both ckks and bfv
	keys := KeyGen(params, a)

	uid := hashPublicKey(keys.publicKey.key)

	keys.publicKey.peerID = uid
	keys.secretKey.peerID = uid
	keys.evalKey.peerID = uid

	encryptor := NewMKEncryptor(keys.publicKey, params, uid)
	decryptor := NewMKDecryptor(params, sigmaSmudging)
	encoder := ckks.NewEncoder(params)
	ringQ := GetRingQ(params)

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
func hashPublicKey(pk [2]*MKDecomposedPoly) uint64 {

	coeffs := pk[0].poly[0].Coeffs // b[0] is the ckks public key
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

	for i, v := range ciphertext.peerIDs {

		if v == participant.id {
			return ciphertext.ciphertexts.Value()[i+1]
		}
	}

	return nil

}
