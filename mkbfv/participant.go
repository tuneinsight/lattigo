package mkbfv

import (
	"encoding/binary"
	"hash/fnv"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKParticipant is a type for participants in a multy key bfv scheme
type MKParticipant interface {
	GetID() uint64
	GetEvaluationKey() *mkrlwe.MKEvaluationKey
	GetPublicKey() *mkrlwe.MKPublicKey
	Encrypt(values []uint64) *MKCiphertext
	Decrypt(cipher *MKCiphertext, partialDecryptions []*ring.Poly) []uint64
	GetPartialDecryption(ciphertext *MKCiphertext) *ring.Poly
	GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey
}

type mkParticipant struct {
	id        uint64
	encryptor MKEncryptor
	decryptor mkrlwe.MKDecryptor
	params    *bfv.Parameters
	keys      *mkrlwe.MKKeys
	encoder   bfv.Encoder
	ringQ     *ring.Ring
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
func (participant *mkParticipant) Encrypt(values []uint64) *MKCiphertext {
	if values == nil || len(values) <= 0 {
		panic("Cannot encrypt uninitialized or empty values")
	}

	pt := newPlaintext(values, participant.encoder, participant.params)

	return participant.encryptor.EncryptMK(pt)
}

// Decrypt returns the decryption of the ciphertext given the partial decryption
func (participant *mkParticipant) Decrypt(cipher *MKCiphertext, partialDecryptions []*ring.Poly) []uint64 {

	if cipher == nil || cipher.Ciphertexts == nil {
		panic("Cannot decrypt uninitialized ciphertext")
	}

	if partialDecryptions == nil || len(partialDecryptions) < 1 {
		panic("Decryption necessitates at least one partialy decrypted ciphertext")
	}

	c0 := cipher.Ciphertexts.Value[0]

	//pass c0 in NTT
	participant.ringQ.NTT(c0, c0)

	decrypted := participant.decryptor.MergeDec(c0, uint64(len(participant.ringQ.Modulus)-1), partialDecryptions)

	//pass result out of NTT domain
	participant.ringQ.InvNTT(decrypted, decrypted)

	pt := bfv.NewPlaintext(*participant.params)
	pt.SetValue(decrypted)

	return participant.encoder.DecodeUintNew(pt)
}

// GetPartialDecryption returns the partial decryption of an element in the ciphertext
// this function should only be used by participants that were involved in the given ciphertext
func (participant *mkParticipant) GetPartialDecryption(ciphertext *MKCiphertext) *ring.Poly {

	cipherPart := participant.getCiphertextPart(ciphertext)

	if cipherPart == nil {
		panic("Participant is not involved in the given ciphertext. Partial decryption impossible.")
	}

	//pass ciphertext in NTT before partial decryption
	participant.ringQ.NTT(cipherPart, cipherPart)

	return participant.decryptor.PartDec(cipherPart, uint64(len(participant.ringQ.Modulus)-1), participant.keys.SecretKey)
}

// NewParticipant creates a participant for the multi key bfv scheme
// the bfv parameters as well as the standard deviation used for partial decryption must be provided
func NewParticipant(params *bfv.Parameters, sigmaSmudging float64, crs *mkrlwe.MKDecomposedPoly) MKParticipant {

	if crs == nil || params == nil {
		panic("Uninitialized parameters. Cannot create new participant")
	}

	if sigmaSmudging < params.Sigma() {
		panic("Sigma must be at least greater than the standard deviation of the gaussian distribution")
	}

	if len(crs.Poly) != int(params.Beta()) {
		panic("CRS must be the same dimention as returned by the function bfv.Parameters.Beta()")
	}

	keys := mkrlwe.KeyGen(&params.Parameters, mkrlwe.CopyNewDecomposed(crs))

	uid := hashPublicKey(keys.PublicKey.Key)

	keys.PublicKey.PeerID = uid
	keys.SecretKey.PeerID = uid
	keys.EvalKey.PeerID = uid

	encryptor := NewMKEncryptor(keys.PublicKey, params, uid)
	decryptor := mkrlwe.NewMKDecryptor(&params.Parameters, sigmaSmudging)
	encoder := bfv.NewEncoder(*params)
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

	coeffs := pk[0].Poly[0].Coeffs // b[0] is the bfv public key
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
func newPlaintext(value []uint64, encoder bfv.Encoder, params *bfv.Parameters) *bfv.Plaintext {

	plaintext := bfv.NewPlaintext(*params)

	// Encode
	encoder.EncodeUint(value, plaintext)

	return plaintext
}

// GetRotationKeys returns the rotation key set associated with the given rotation
func (participant *mkParticipant) GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey {

	galEl := participant.params.GaloisElementForColumnRotationBy(rot)

	evalKey := mkrlwe.GaloisEvaluationKeyGen(galEl, participant.keys.SecretKey, &participant.params.Parameters)
	evalKey.PeerID = participant.id

	return evalKey
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
