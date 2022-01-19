package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Encryptor a generic RLWE encryption interface.
type Encryptor interface {
	Encrypt(pt *Plaintext, ct *Ciphertext)
	EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ct *Ciphertext)
	ShallowCopy() Encryptor
	WithKey(key interface{}) Encryptor
}

type pkEncryptor struct {
	encryptor
	pk *PublicKey
}

type skEncryptor struct {
	encryptor
	sk *SecretKey
}

type encryptor struct {
	*encryptorBase
	*encryptorSamplers
	*encryptorBuffers
	basisextender *ring.BasisExtender
}

// NewEncryptor creates a new Encryptor
// Accepts either a secret-key or a public-key
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	switch key := key.(type) {
	case *PublicKey:
		if key.Value[0].Q.Degree() != params.N() || key.Value[1].Q.Degree() != params.N() {
			panic("cannot newEncryptor: pk ring degree does not match params ring degree")
		}
		return &pkEncryptor{newEncryptor(params), key}
	case *SecretKey:
		if key.Value.Q.Degree() != params.N() {
			panic("cannot newEncryptor: sk ring degree does not match params ring degree")
		}
		return &skEncryptor{newEncryptor(params), key}
	case nil:
		enc := newEncryptor(params)
		return &enc
	default:
		panic("key must be either *rlwe.PublicKey, *rlwe.SecretKey or nil")
	}
}

func newEncryptor(params Parameters) encryptor {

	var bc *ring.BasisExtender
	if params.PCount() != 0 {
		bc = ring.NewBasisExtender(params.RingQ(), params.RingP())
	}

	return encryptor{
		encryptorBase:     newEncryptorBase(params),
		encryptorSamplers: newEncryptorSamplers(params),
		encryptorBuffers:  newEncryptorBuffers(params),
		basisextender:     bc,
	}
}

// encryptorBase is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptorBase struct {
	params Parameters
}

func newEncryptorBase(params Parameters) *encryptorBase {
	return &encryptorBase{params}
}

type encryptorSamplers struct {
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSampler  *ring.UniformSampler
}

func newEncryptorSamplers(params Parameters) *encryptorSamplers {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ringQ := params.RingQ()
	return &encryptorSamplers{
		gaussianSampler: ring.NewGaussianSampler(prng, ringQ, params.Sigma(), int(6*params.Sigma())),
		ternarySampler:  ring.NewTernarySampler(prng, ringQ, 0.5, false),
		uniformSampler:  ring.NewUniformSampler(prng, ringQ),
	}
}

type encryptorBuffers struct {
	poolQ [1]*ring.Poly
	poolP [3]*ring.Poly
}

func newEncryptorBuffers(params Parameters) *encryptorBuffers {

	ringQ := params.RingQ()
	ringP := params.RingP()

	var poolP [3]*ring.Poly
	if params.PCount() != 0 {
		poolP = [3]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return &encryptorBuffers{
		poolQ: [1]*ring.Poly{ringQ.NewPoly()},
		poolP: poolP,
	}
}

// Encrypt encrypts the input plaintext using the stored public-key and write the result on ct.
// The encryption procedure first samples an new encryption of zero under the public-key and
// then adds the plaintext.
// The encryption procedures depends on the parameters. If the auxiliary modulus P is defined,
// then the encryption of zero is sampled in QP before being rescaled by P, else it is directly
// samples in Q.
func (enc *pkEncryptor) Encrypt(pt *Plaintext, ct *Ciphertext) {
	enc.encrypt(pt, enc.pk, ct)
}

// EncryptFromCRP is not defined when using a public-key. This method will panic.
func (enc *pkEncryptor) EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ct *Ciphertext) {
	enc.encryptFromCRP(pt, enc.pk, ct)
}

// Encrypt encrypts the input plaintext and write the result on ct.
func (enc *skEncryptor) Encrypt(pt *Plaintext, ct *Ciphertext) {
	enc.encrypt(pt, enc.sk, ct)
}

// EncryptFromCRP encrypts the input plaintext and writes the result on ct.
// The encryption algorithm depends on the implementor.
func (enc *skEncryptor) EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ct *Ciphertext) {
	enc.encryptFromCRP(pt, enc.sk, ct)
}

// Encrypt is not defined when the key is nil. This method will panic.
func (enc *encryptor) Encrypt(pt *Plaintext, ct *Ciphertext) {
	panic("cannot encrypt, key is nil")
}

// EncryptFromCRP is not defined when the key is nil. This method will panic.
func (enc *encryptor) EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ct *Ciphertext) {
	panic("cannot encrypt, key is nil")
}

// ShallowCopy creates a shallow copy of this pkEncryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *pkEncryptor) ShallowCopy() Encryptor {
	return &pkEncryptor{*enc.encryptor.ShallowCopy().(*encryptor), enc.pk}
}

// ShallowCopy creates a shallow copy of this skEncryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *skEncryptor) ShallowCopy() Encryptor {
	return &skEncryptor{*enc.encryptor.ShallowCopy().(*encryptor), enc.sk}
}

// ShallowCopy creates a shallow copy of this encryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *encryptor) ShallowCopy() Encryptor {

	var bc *ring.BasisExtender
	if enc.params.PCount() != 0 {
		bc = enc.basisextender.ShallowCopy()
	}

	return &encryptor{
		encryptorBase:     enc.encryptorBase,
		encryptorSamplers: newEncryptorSamplers(enc.params),
		encryptorBuffers:  newEncryptorBuffers(enc.params),
		basisextender:     bc,
	}
}

// WithKey creates a shallow copy of this pkEncryptor with a new key in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *pkEncryptor) WithKey(key interface{}) Encryptor {
	return enc.ShallowCopy().(*pkEncryptor).setKey(key)
}

// WithKey creates a shallow copy of this skEncryptor with a new key in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *skEncryptor) WithKey(key interface{}) Encryptor {
	return enc.ShallowCopy().(*skEncryptor).setKey(key)
}

// WithKey creates a shallow copy of this encryptor with a new key in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *encryptor) WithKey(key interface{}) Encryptor {
	return enc.ShallowCopy().(*encryptor).setKey(key)
}

func (enc *encryptor) encrypt(plaintext *Plaintext, key interface{}, ciphertext *Ciphertext) {
	switch key := key.(type) {
	case *PublicKey:

		enc.uniformSampler.ReadLvl(utils.MinInt(plaintext.Level(), ciphertext.Level()), ciphertext.Value[1])

		if enc.basisextender != nil {
			enc.encryptPk(plaintext, key, ciphertext)
		} else {
			enc.encryptPkNoP(plaintext, key, ciphertext)
		}

	case *SecretKey:

		enc.uniformSampler.ReadLvl(utils.MinInt(plaintext.Level(), ciphertext.Level()), ciphertext.Value[1])

		enc.encryptSk(plaintext, key, ciphertext)

	case nil:
		panic("key is nil")

	default:
		panic("key must be either rlwe.PublicKey or rlwe.SecretKey")
	}
}

func (enc *encryptor) encryptFromCRP(plaintext *Plaintext, key interface{}, crp *ring.Poly, ciphertext *Ciphertext) {
	switch key := key.(type) {
	case *PublicKey:

		panic("Cannot encrypt with CRP using a public-key")

	case *SecretKey:

		if key.Value.Q.Degree() != enc.params.N() {
			panic("cannot newEncryptor: sk ring degree does not match params ring degree")
		}

		ring.CopyValues(crp, ciphertext.Value[1])

		enc.encryptSk(plaintext, key, ciphertext)

	default:
		panic("key must be either rlwe.PublicKey or rlwe.SecretKey")
	}
}

func (enc *encryptor) encryptPk(plaintext *Plaintext, pk *PublicKey, ciphertext *Ciphertext) {
	ringQ := enc.params.RingQ()
	ringQP := enc.params.RingQP()

	levelQ := utils.MinInt(plaintext.Level(), ciphertext.Level())
	levelP := 0

	poolQ0 := enc.poolQ[0]
	poolP0 := enc.poolP[0]
	poolP1 := enc.poolP[1]
	poolP2 := enc.poolP[2]

	// We sample a R-WLE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)

	ciphertextNTT := ciphertext.Value[0].IsNTT

	u := PolyQP{Q: poolQ0, P: poolP2}

	enc.ternarySampler.ReadLvl(levelQ, u.Q)
	ringQP.ExtendBasisSmallNormAndCenter(u.Q, levelP, nil, u.P)

	// (#Q + #P) NTT
	ringQP.NTTLvl(levelQ, levelP, u, u)
	ringQP.MFormLvl(levelQ, levelP, u, u)

	ct0QP := PolyQP{Q: ciphertext.Value[0], P: poolP0}
	ct1QP := PolyQP{Q: ciphertext.Value[1], P: poolP1}

	// ct0 = u*pk0
	// ct1 = u*pk1
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, u, pk.Value[0], ct0QP)
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, u, pk.Value[1], ct1QP)

	// 2*(#Q + #P) NTT
	ringQP.InvNTTLvl(levelQ, levelP, ct0QP, ct0QP)
	ringQP.InvNTTLvl(levelQ, levelP, ct1QP, ct1QP)

	e := PolyQP{Q: poolQ0, P: poolP2}

	enc.gaussianSampler.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.AddLvl(levelQ, levelP, ct0QP, e, ct0QP)

	enc.gaussianSampler.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.AddLvl(levelQ, levelP, ct1QP, e, ct1QP)

	// ct0 = (u*pk0 + e0)/P
	enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct0QP.Q, ct0QP.P, ct0QP.Q)

	// ct1 = (u*pk1 + e1)/P
	enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct1QP.Q, ct1QP.P, ct1QP.Q)

	if ciphertextNTT {

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		}

		// 2*#Q NTT
		ringQ.NTTLvl(levelQ, ciphertext.Value[0], ciphertext.Value[0])
		ringQ.NTTLvl(levelQ, ciphertext.Value[1], ciphertext.Value[1])

		if plaintext.Value.IsNTT {
			// ct0 = (u*pk0 + e0)/P + m
			ringQ.AddLvl(levelQ, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		}

	} else {

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		} else {
			ringQ.InvNTTLvl(levelQ, plaintext.Value, poolQ0)
			ringQ.AddLvl(levelQ, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
		}
	}

	ciphertext.Value[1].IsNTT = ciphertext.Value[0].IsNTT
	ciphertext.Value[0].Coeffs = ciphertext.Value[0].Coeffs[:levelQ+1]
	ciphertext.Value[1].Coeffs = ciphertext.Value[1].Coeffs[:levelQ+1]
}

func (enc *encryptor) encryptPkNoP(plaintext *Plaintext, pk *PublicKey, ciphertext *Ciphertext) {
	levelQ := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := enc.poolQ[0]

	ringQ := enc.params.RingQ()

	ciphertextNTT := ciphertext.Value[0].IsNTT

	enc.ternarySampler.ReadLvl(levelQ, poolQ0)
	ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
	ringQ.MFormLvl(levelQ, poolQ0, poolQ0)

	// ct0 = u*pk0
	ringQ.MulCoeffsMontgomeryLvl(levelQ, poolQ0, pk.Value[0].Q, ciphertext.Value[0])
	// ct1 = u*pk1
	ringQ.MulCoeffsMontgomeryLvl(levelQ, poolQ0, pk.Value[1].Q, ciphertext.Value[1])

	if ciphertextNTT {

		// ct1 = u*pk1 + e1
		enc.gaussianSampler.ReadLvl(levelQ, poolQ0)
		ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
		ringQ.AddLvl(levelQ, ciphertext.Value[1], poolQ0, ciphertext.Value[1])

		// ct0 = u*pk0 + e0
		enc.gaussianSampler.ReadLvl(levelQ, poolQ0)

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, poolQ0, plaintext.Value, poolQ0)
			ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
			ringQ.AddLvl(levelQ, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
		} else {
			ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
			ringQ.AddLvl(levelQ, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
			ringQ.AddLvl(levelQ, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		}

	} else {

		ringQ.InvNTTLvl(levelQ, ciphertext.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(levelQ, ciphertext.Value[1], ciphertext.Value[1])

		// ct[0] = pk[0]*u + e0
		enc.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), ciphertext.Value[0])

		// ct[1] = pk[1]*u + e1
		enc.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), ciphertext.Value[1])

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		} else {
			ringQ.InvNTTLvl(levelQ, plaintext.Value, poolQ0)
			ringQ.AddLvl(levelQ, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
		}
	}

	ciphertext.Value[1].IsNTT = ciphertext.Value[0].IsNTT

	ciphertext.Value[0].Coeffs = ciphertext.Value[0].Coeffs[:levelQ+1]
	ciphertext.Value[1].Coeffs = ciphertext.Value[1].Coeffs[:levelQ+1]
}

func (enc *encryptor) encryptSk(plaintext *Plaintext, sk *SecretKey, ciphertext *Ciphertext) {

	ringQ := enc.params.RingQ()

	levelQ := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := enc.poolQ[0]

	ciphertextNTT := ciphertext.Value[0].IsNTT

	ringQ.MulCoeffsMontgomeryLvl(levelQ, ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
	ringQ.NegLvl(levelQ, ciphertext.Value[0], ciphertext.Value[0])

	if ciphertextNTT {

		enc.gaussianSampler.ReadLvl(levelQ, poolQ0)

		if plaintext.Value.IsNTT {
			ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
			ringQ.AddLvl(levelQ, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
			ringQ.AddLvl(levelQ, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		} else {
			ringQ.AddLvl(levelQ, poolQ0, plaintext.Value, poolQ0)
			ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
			ringQ.AddLvl(levelQ, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
		}

		ciphertext.Value[0].IsNTT = true
		ciphertext.Value[1].IsNTT = true

	} else {

		if plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
			ringQ.InvNTTLvl(levelQ, ciphertext.Value[0], ciphertext.Value[0])

		} else {
			ringQ.InvNTTLvl(levelQ, ciphertext.Value[0], ciphertext.Value[0])
			ringQ.AddLvl(levelQ, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		}

		enc.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), ciphertext.Value[0])

		ringQ.InvNTTLvl(levelQ, ciphertext.Value[1], ciphertext.Value[1])

		ciphertext.Value[0].IsNTT = false
		ciphertext.Value[1].IsNTT = false

	}

	ciphertext.Value[0].Coeffs = ciphertext.Value[0].Coeffs[:levelQ+1]
	ciphertext.Value[1].Coeffs = ciphertext.Value[1].Coeffs[:levelQ+1]
}

func (enc *pkEncryptor) setKey(key interface{}) Encryptor {
	return enc.encryptor.setKey(key)
}

func (enc *skEncryptor) setKey(key interface{}) Encryptor {
	return enc.encryptor.setKey(key)
}

func (enc *encryptor) setKey(key interface{}) Encryptor {
	switch key := key.(type) {
	case *PublicKey:
		if key.Value[0].Q.Degree() != enc.params.N() || key.Value[1].Q.Degree() != enc.params.N() {
			panic("cannot newEncryptor: pk ring degree does not match params ring degree")
		}
		return &pkEncryptor{*enc, key}
	case *SecretKey:
		if key.Value.Q.Degree() != enc.params.N() {
			panic("cannot newEncryptor: sk ring degree does not match params ring degree")
		}
		return &skEncryptor{*enc, key}
	default:
		panic("key must be either *rlwe.PublicKey or *rlwe.SecretKey")
	}
}
