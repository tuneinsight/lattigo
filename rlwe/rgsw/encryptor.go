package rgsw

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/gadget"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Encryptor a generic RLWE encryption interface.
type Encryptor interface {
	Encrypt(pt *rlwe.Plaintext, ct *Ciphertext)
	ShallowCopy() Encryptor
	WithKey(key interface{}) Encryptor
}

type encryptor struct {
	*encryptorBase
	*encryptorSamplers
	*encryptorBuffers
}

type skEncryptor struct {
	encryptor
	sk *rlwe.SecretKey
}

// NewEncryptor creates a new Encryptor
// Accepts either a secret-key or a public-key.
func NewEncryptor(params rlwe.Parameters, key interface{}) Encryptor {
	enc := newEncryptor(params)
	return enc.setKey(key)
}

func newEncryptor(params rlwe.Parameters) encryptor {
	return encryptor{
		encryptorBase:     newEncryptorBase(params),
		encryptorSamplers: newEncryptorSamplers(params),
		encryptorBuffers:  newEncryptorBuffers(params),
	}
}

// encryptorBase is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptorBase struct {
	params rlwe.Parameters
}

func newEncryptorBase(params rlwe.Parameters) *encryptorBase {
	return &encryptorBase{params}
}

type encryptorSamplers struct {
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSamplerQ *ring.UniformSampler
	uniformSamplerP *ring.UniformSampler
}

func newEncryptorSamplers(params rlwe.Parameters) *encryptorSamplers {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var uniformSamplerP *ring.UniformSampler
	if params.PCount() != 0 {
		uniformSamplerP = ring.NewUniformSampler(prng, params.RingP())
	}

	return &encryptorSamplers{
		gaussianSampler: ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
		ternarySampler:  ring.NewTernarySamplerWithHammingWeight(prng, params.RingQ(), params.HammingWeight(), false),
		uniformSamplerQ: ring.NewUniformSampler(prng, params.RingQ()),
		uniformSamplerP: uniformSamplerP,
	}
}

type encryptorBuffers struct {
	poolQP ringqp.Poly
}

func newEncryptorBuffers(params rlwe.Parameters) *encryptorBuffers {
	return &encryptorBuffers{
		poolQP: params.RingQP().NewPoly(),
	}
}

// ShallowCopy creates a shallow copy of this skEncryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *skEncryptor) ShallowCopy() Encryptor {
	return &skEncryptor{*enc.encryptor.ShallowCopy(), enc.sk}
}

// ShallowCopy creates a shallow copy of this encryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *encryptor) ShallowCopy() *encryptor {
	return &encryptor{
		encryptorBase:     enc.encryptorBase,
		encryptorSamplers: newEncryptorSamplers(enc.params),
		encryptorBuffers:  newEncryptorBuffers(enc.params),
	}
}

// WithKey creates a shallow copy of this encryptor with a new key in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *encryptor) WithKey(key interface{}) Encryptor {
	return enc.ShallowCopy().setKey(key)
}

func (enc *encryptor) setKey(key interface{}) Encryptor {
	switch key := key.(type) {
	case *rlwe.SecretKey:
		if key.Value.Q.Degree() != enc.params.N() {
			panic("cannot setKey: sk ring degree does not match params ring degree")
		}
		return &skEncryptor{*enc, key}
	default:
		panic("cannot setKey: key must be *rlwe.SecretKey")
	}
}

func (enc *skEncryptor) Encrypt(pt *rlwe.Plaintext, ct *Ciphertext) {

	params := enc.params
	ringQ := params.RingQ()
	levelQ := ct.LevelQ()
	levelP := ct.LevelP()

	decompRNS := params.DecompRNS(levelQ, levelP)
	decompBIT := params.DecompBIT(levelQ, levelP)

	for j := 0; j < decompBIT; j++ {
		for i := 0; i < decompRNS; i++ {
			enc.encryptZeroSymetricQP(levelQ, levelP, enc.sk.Value, true, true, true, ct.Value[0].Value[i][j])
			enc.encryptZeroSymetricQP(levelQ, levelP, enc.sk.Value, true, true, true, ct.Value[1].Value[i][j])
		}
	}

	if pt != nil {
		ringQ.MFormLvl(levelQ, pt.Value, enc.poolQP.Q)
		if !pt.Value.IsNTT {
			ringQ.NTTLvl(levelQ, enc.poolQP.Q, enc.poolQP.Q)
		}
		gadget.AddPolyToGadgetMatrix(
			enc.poolQP.Q,
			[]gadget.Ciphertext{ct.Value[0], ct.Value[1]},
			*params.RingQP(),
			params.LogBase2(),
			enc.poolQP.Q)
	}
}

func (enc *encryptor) encryptZeroSymetricQP(levelQ, levelP int, sk ringqp.Poly, sample, montgomery, ntt bool, ct [2]ringqp.Poly) {

	params := enc.params
	ringQP := params.RingQP()

	hasModulusP := ct[0].P != nil

	if ntt {
		enc.gaussianSampler.ReadLvl(levelQ, ct[0].Q)

		if hasModulusP {
			ringQP.ExtendBasisSmallNormAndCenter(ct[0].Q, levelP, nil, ct[0].P)
		}

		ringQP.NTTLvl(levelQ, levelP, ct[0], ct[0])
	}

	if sample {
		enc.uniformSamplerQ.ReadLvl(levelQ, ct[1].Q)

		if hasModulusP {
			enc.uniformSamplerP.ReadLvl(levelP, ct[1].P)
		}
	}

	ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, ct[1], sk, ct[0])

	if !ntt {
		ringQP.InvNTTLvl(levelQ, levelP, ct[0], ct[0])
		ringQP.InvNTTLvl(levelQ, levelP, ct[1], ct[1])

		e := enc.poolQP
		enc.gaussianSampler.ReadLvl(levelQ, e.Q)

		if hasModulusP {
			ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
		}

		ringQP.AddLvl(levelQ, levelP, ct[0], e, ct[0])
	}

	if montgomery {
		ringQP.MFormLvl(levelQ, levelP, ct[0], ct[0])
		ringQP.MFormLvl(levelQ, levelP, ct[1], ct[1])
	}
}
