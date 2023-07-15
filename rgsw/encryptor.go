package rgsw

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// Encryptor is a type for encrypting RGSW ciphertexts. It implements the rlwe.Encryptor
// interface overriding the `Encrypt` and `EncryptZero` methods to accept rgsw.Ciphertext
// types in addition to ciphertexts types in the rlwe package.
type Encryptor struct {
	rlwe.EncryptorInterface

	params rlwe.Parameters
	buffQP ringqp.Poly
}

// NewEncryptor creates a new Encryptor type. Note that only secret-key encryption is
// supported at the moment.
func NewEncryptor[T *rlwe.SecretKey | *rlwe.PublicKey](params rlwe.Parameters, key T) (*Encryptor, error) {
	enc, err := rlwe.NewEncryptor(params, key)
	return &Encryptor{enc, params, params.RingQP().NewPoly()}, err
}

// Encrypt encrypts a plaintext pt into a ciphertext ct, which can be a rgsw.Ciphertext
// or any of the `rlwe` cipheretxt types.
func (enc Encryptor) Encrypt(pt *rlwe.Plaintext, ct interface{}) (err error) {

	var rgswCt *Ciphertext
	var isRGSW bool
	if rgswCt, isRGSW = ct.(*Ciphertext); !isRGSW {
		return enc.EncryptorInterface.Encrypt(pt, ct)
	}

	if err = enc.EncryptZero(rgswCt); err != nil {
		return
	}

	levelQ := rgswCt.LevelQ()
	ringQ := enc.params.RingQ().AtLevel(levelQ)

	if pt != nil {

		if !pt.IsNTT {
			ringQ.NTT(pt.Value, enc.buffQP.Q)

			if !pt.IsMontgomery {
				ringQ.MForm(enc.buffQP.Q, enc.buffQP.Q)
			}

		} else {
			if !pt.IsMontgomery {
				ringQ.MForm(pt.Value, enc.buffQP.Q)
			} else {
				ring.CopyLvl(levelQ, enc.buffQP.Q, pt.Value)
			}
		}

		return rlwe.AddPolyTimesGadgetVectorToGadgetCiphertext(
			enc.buffQP.Q,
			[]rlwe.GadgetCiphertext{rgswCt.Value[0], rgswCt.Value[1]},
			*enc.params.RingQP(),
			enc.buffQP.Q)
	}

	return
}

// EncryptZero generates an encryption of zero into a ciphertext ct, which can be a rgsw.Ciphertext
// or any of the `rlwe` cipheretxt types.
func (enc Encryptor) EncryptZero(ct interface{}) (err error) {

	var rgswCt *Ciphertext
	var isRGSW bool
	if rgswCt, isRGSW = ct.(*Ciphertext); !isRGSW {
		return enc.EncryptorInterface.EncryptZero(ct)
	}

	decompRNS := rgswCt.Value[0].DecompRNS()
	decompPw2 := rgswCt.Value[0].DecompPw2()

	for j := 0; j < decompPw2; j++ {
		for i := 0; i < decompRNS; i++ {

			if err = enc.EncryptorInterface.EncryptZero(rlwe.OperandQP{MetaData: rlwe.MetaData{IsNTT: true, IsMontgomery: true}, Value: []ringqp.Poly(rgswCt.Value[0].Value[i][j])}); err != nil {
				return
			}

			if err = enc.EncryptorInterface.EncryptZero(rlwe.OperandQP{MetaData: rlwe.MetaData{IsNTT: true, IsMontgomery: true}, Value: []ringqp.Poly(rgswCt.Value[1].Value[i][j])}); err != nil {
				return
			}
		}
	}

	return
}

// ShallowCopy creates a shallow copy of this Encryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc Encryptor) ShallowCopy() *Encryptor {
	return &Encryptor{EncryptorInterface: enc.EncryptorInterface.ShallowCopy(), params: enc.params, buffQP: enc.params.RingQP().NewPoly()}
}
