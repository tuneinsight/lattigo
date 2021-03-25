package mkbfv

import (
	"encoding/binary"
	"hash/fnv"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKParticipant is a type for participants in a multy key bfv scheme
type MKParticipant interface {
	GetID() uint64
	Encrypt(values []uint64) *MKCiphertext
	Decrypt(c0 *ring.Poly, partialDecryptions []*ring.Poly) []uint64
	GetPartialDecryption(ciphertext *MKCiphertext) *ring.Poly
}

type mkParticipant struct {
	id        uint64
	encryptor MKEncryptor
	decryptor MKDecryptor
	keys      *MKKeys
	encoder   bfv.Encoder
	ringQ     *ring.Ring
}

// GetID returns the id of the participant
func (participant *mkParticipant) GetID() uint64 {
	return participant.id
}

// Encrypt constructs a ciphertext from the given values
func (participant *mkParticipant) Encrypt(values []uint64) *MKCiphertext {
	return participant.encryptor.EncryptMK(newPlaintext(values, participant.ringQ, participant.encoder))
}

// Decrypt returns the decryption of the ciphertext given the partial decryption
func (participant *mkParticipant) Decrypt(c0 *ring.Poly, partialDecryptions []*ring.Poly) []uint64 {

	decrypted := participant.decryptor.MergeDec(c0, partialDecryptions)

	return participant.encoder.DecodeUintNew(decrypted)
}

// GetPartialDecryption returns the partial decryption of an element in the ciphertext
// this function should only be used by participants that were involved in the given ciphertext
func (participant *mkParticipant) GetPartialDecryption(ciphertext *MKCiphertext) *ring.Poly {

	cipherPart := participant.getCiphertextPart(ciphertext)

	if cipherPart == nil {
		panic("Participant is not involved in the given ciphertext. Partial decryption impossible.")
	}
	return participant.decryptor.PartDec(cipherPart, participant.keys.secretKey)
}

// NewParticipant creates a participant for the multi key bfv scheme
// the bfv parameters as well as the standard deviation used for partial decryption must be provided
func NewParticipant(params *bfv.Parameters, sigmaSmudging float64) MKParticipant {

	a := GenCommonPublicParam(params)
	keys := KeyGen(params, a)

	uid := hashPublicKey(keys.publicKey.key)

	keys.publicKey.peerID = uid
	keys.secretKey.peerID = uid
	keys.evalKey.peerID = uid

	encryptor := NewMKEncryptor(keys.publicKey, params, uid)
	decryptor := NewMKDecryptor(params, sigmaSmudging)
	encoder := bfv.NewEncoder(params)
	ringQ := GetRingQ(params)

	return &mkParticipant{
		id:        uid,
		encryptor: encryptor,
		decryptor: decryptor,
		keys:      keys,
		encoder:   encoder,
		ringQ:     ringQ,
	}
}

// computes the hash of the public key using the FNV hashing algorithm
func hashPublicKey(pk [2]*MKDecomposedPoly) uint64 {

	coeffs := pk[0].poly[0].Coeffs // b[0] is the bfv public key
	h64 := fnv.New64()

	for _, v := range coeffs {
		for _, vIn := range v {
			h64.Write(uintToBytes(vIn))
		}
	}

	return h64.Sum64()
}

func uintToBytes(i uint64) []byte {

	b := make([]byte, 8)
	binary.LittleEndian.PutUint64(b, uint64(i))

	return b
}

// newPlaintext initializes a new bfv Plaintext with an encoded slice of uint64
func newPlaintext(value []uint64, ringQ *ring.Ring, encoder bfv.Encoder) *bfv.Plaintext {

	plaintext := new(bfv.Element)

	ptValues := make([]*ring.Poly, 1)

	ptValues[0] = ringQ.NewPoly()
	plaintext.SetValue(ptValues)

	// Encode
	encoder.EncodeUint(value, plaintext.Plaintext())

	return plaintext.Plaintext()
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
