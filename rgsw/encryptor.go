package rgsw

import (
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

type Encryptor struct {
	rlwe.Encryptor

	params rlwe.Parameters
	buffQP ringqp.Poly
}

func NewEncryptor(params rlwe.Parameters, sk *rlwe.SecretKey) *Encryptor {
	return &Encryptor{rlwe.NewEncryptor(params, sk), params, params.RingQP().NewPoly()}
}

func (enc *Encryptor) Encrypt(pt *rlwe.Plaintext, ct interface{}) {

	var rgswCt *Ciphertext
	var isRGSW bool
	if rgswCt, isRGSW = ct.(*Ciphertext); !isRGSW {
		enc.Encryptor.Encrypt(pt, ct)
	}

	params := enc.params
	ringQ := params.RingQ()
	levelQ := rgswCt.LevelQ()
	levelP := rgswCt.LevelP()

	decompRNS := params.DecompRNS(levelQ, levelP)
	decompBIT := params.DecompBIT(levelQ, levelP)

	for j := 0; j < decompBIT; j++ {
		for i := 0; i < decompRNS; i++ {
			enc.EncryptZero(&rgswCt.Value[0].Value[i][j])
			enc.EncryptZero(&rgswCt.Value[1].Value[i][j])
		}
	}

	if pt != nil {
		ringQ.MFormLvl(levelQ, pt.Value, enc.buffQP.Q)
		if !pt.Value.IsNTT {
			ringQ.NTTLvl(levelQ, enc.buffQP.Q, enc.buffQP.Q)
		}
		rlwe.AddPolyTimesGadgetVectorToGadgetCiphertext(
			enc.buffQP.Q,
			[]rlwe.GadgetCiphertext{rgswCt.Value[0], rgswCt.Value[1]},
			*params.RingQP(),
			params.LogBase2(),
			enc.buffQP.Q)
	}
}
