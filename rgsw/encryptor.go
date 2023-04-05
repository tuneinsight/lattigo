package rgsw

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// Encryptor is a type for encrypting RGSW ciphertexts. It implements the rlwe.Encryptor
// interface overriding the `Encrypt` and `EncryptZero` methods to accept rgsw.Ciphertext
// types in addition to ciphertexts types in the rlwe package.
type Encryptor struct {
	rlwe.Encryptor

	params rlwe.Parameters
	buffQP ringqp.Poly
}

// NewEncryptor creates a new Encryptor type. Note that only secret-key encryption is
// supported at the moment.
func NewEncryptor(params rlwe.Parameters, sk *rlwe.SecretKey) *Encryptor {
	return &Encryptor{rlwe.NewEncryptor(params, sk), params, *params.RingQP().NewPoly()}
}

// Encrypt encrypts a plaintext pt into a ciphertext ct, which can be a rgsw.Ciphertext
// or any of the `rlwe` cipheretxt types.
func (enc *Encryptor) Encrypt(pt *rlwe.Plaintext, ct interface{}) {

	var rgswCt *Ciphertext
	var isRGSW bool
	if rgswCt, isRGSW = ct.(*Ciphertext); !isRGSW {
		enc.Encryptor.Encrypt(pt, ct)
		return
	}

	enc.EncryptZero(rgswCt)

	levelQ := rgswCt.LevelQ()
	ringQ := enc.params.RingQ().AtLevel(levelQ)

	if pt != nil {
		ringQ.MForm(pt.Value, enc.buffQP.Q)
		if !pt.IsNTT {
			ringQ.NTT(enc.buffQP.Q, enc.buffQP.Q)
		}
		rlwe.AddPolyTimesGadgetVectorToGadgetCiphertext(
			enc.buffQP.Q,
			[]rlwe.GadgetCiphertext{rgswCt.Value[0], rgswCt.Value[1]},
			*enc.params.RingQP(),
			enc.params.Pow2Base(),
			enc.buffQP.Q)
	}
}

// EncryptZero generates an encryption of zero into a ciphertext ct, which can be a rgsw.Ciphertext
// or any of the `rlwe` cipheretxt types.
func (enc *Encryptor) EncryptZero(ct interface{}) {

	var rgswCt *Ciphertext
	var isRGSW bool
	if rgswCt, isRGSW = ct.(*Ciphertext); !isRGSW {
		enc.Encryptor.EncryptZero(ct)
		return
	}

	levelQ := rgswCt.LevelQ()
	levelP := rgswCt.LevelP()
	decompRNS := enc.params.DecompRNS(levelQ, levelP)
	decompPw2 := enc.params.DecompPw2(levelQ, levelP)

	for j := 0; j < decompPw2; j++ {
		for i := 0; i < decompRNS; i++ {
			enc.Encryptor.EncryptZero(rgswCt.Value[0].Value[i][j])
			enc.Encryptor.EncryptZero(rgswCt.Value[1].Value[i][j])
		}
	}
}

// ShallowCopy creates a shallow copy of this Encryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *Encryptor) ShallowCopy() *Encryptor {
	return &Encryptor{Encryptor: enc.Encryptor.ShallowCopy(), params: enc.params, buffQP: *enc.params.RingQP().NewPoly()}
}
