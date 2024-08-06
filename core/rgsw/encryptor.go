package rgsw

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
)

// Encryptor is a type for encrypting RGSW ciphertexts. It implements the [rlwe.Encryptor]
// interface overriding the [rlwe.Encryptor.Encrypt] and [rlwe.Encryptor.EncryptZero] methods to accept [rgsw.Ciphertext]
// types in addition to ciphertexts types in the rlwe package.
type Encryptor struct {
	*rlwe.Encryptor
	buffQP ringqp.Poly
}

// NewEncryptor creates a new Encryptor type. Note that only secret-key encryption is
// supported at the moment.
func NewEncryptor(params rlwe.ParameterProvider, key rlwe.EncryptionKey) *Encryptor {
	return &Encryptor{rlwe.NewEncryptor(params, key), params.GetRLWEParameters().RingQP().NewPoly()}
}

// Encrypt encrypts a plaintext pt into a ciphertext ct, which can be a [rgsw.Ciphertext]
// or any of the `rlwe` cipheretxt types.
func (enc Encryptor) Encrypt(pt *rlwe.Plaintext, ct interface{}) (err error) {

	var rgswCt *Ciphertext
	var isRGSW bool
	if rgswCt, isRGSW = ct.(*Ciphertext); !isRGSW {
		return enc.Encryptor.Encrypt(pt, ct)
	}

	if err = enc.EncryptZero(rgswCt); err != nil {
		return
	}

	params := enc.GetRLWEParameters()

	levelQ := rgswCt.LevelQ()
	ringQ := params.RingQ().AtLevel(levelQ)

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
				pt.Value.CopyLvl(levelQ, enc.buffQP.Q)
			}
		}

		if err := rlwe.AddPolyTimesGadgetVectorToGadgetCiphertext(
			enc.buffQP.Q,
			[]rlwe.GadgetCiphertext{rgswCt.Value[0], rgswCt.Value[1]},
			*params.RingQP(),
			enc.buffQP.Q); err != nil {
			// Sanity check, this error should not happen.
			panic(err)
		}
	}

	return nil
}

// EncryptZero generates an encryption of zero into a ciphertext ct, which can be a [rgsw.Ciphertext]
// or any of the `rlwe` ciphertext types.
func (enc Encryptor) EncryptZero(ct interface{}) (err error) {

	var rgswCt *Ciphertext
	var isRGSW bool
	if rgswCt, isRGSW = ct.(*Ciphertext); !isRGSW {
		return enc.Encryptor.EncryptZero(rgswCt)
	}

	BaseRNSDecompositionVectorSize := rgswCt.Value[0].BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := rgswCt.Value[0].BaseTwoDecompositionVectorSize()

	metadata := &rlwe.MetaData{}
	metadata.IsMontgomery = true
	metadata.IsNTT = true

	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		for j := 0; j < BaseTwoDecompositionVectorSize[i]; j++ {

			// extract the RingQ polynomial in case not P moduli have been provided
			if rgswCt.LevelP() == -1 {
				q0 := make([]ring.Poly, len(rgswCt.Value[0].Value[i][j]))
				q1 := make([]ring.Poly, len(rgswCt.Value[1].Value[i][j]))
				for k := range q0 {
					q0[k] = rgswCt.Value[0].Value[i][j][k].Q
					q1[k] = rgswCt.Value[1].Value[i][j][k].Q
				}
				if err = enc.Encryptor.EncryptZero(&rlwe.Ciphertext{Element: rlwe.Element[ring.Poly]{MetaData: metadata, Value: q0}}); err != nil {
					return
				}
				if err = enc.Encryptor.EncryptZero(&rlwe.Ciphertext{Element: rlwe.Element[ring.Poly]{MetaData: metadata, Value: q1}}); err != nil {
					return
				}
			} else {
				if err = enc.Encryptor.EncryptZero(rlwe.Element[ringqp.Poly]{MetaData: metadata, Value: []ringqp.Poly(rgswCt.Value[0].Value[i][j])}); err != nil {
					return
				}
				if err = enc.Encryptor.EncryptZero(rlwe.Element[ringqp.Poly]{MetaData: metadata, Value: []ringqp.Poly(rgswCt.Value[1].Value[i][j])}); err != nil {
					return
				}
			}
		}
	}

	return nil
}

// ShallowCopy creates a shallow copy of this [Encryptor] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc Encryptor) ShallowCopy() *Encryptor {
	return &Encryptor{Encryptor: enc.Encryptor.ShallowCopy(), buffQP: enc.GetRLWEParameters().RingQP().NewPoly()}
}
