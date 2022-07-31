package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Encryptor a generic RLWE encryption interface.
type Encryptor interface {
	Encrypt(pt *Plaintext, ct interface{})
	EncryptZero(ct interface{})

	ShallowCopy() Encryptor
	WithKey(key interface{}) Encryptor
}

// PRNGEncryptor is an interface for encrypting RLWE ciphertexts from a secret-key and
// a pre-determined PRNG. An Encryptor constructed from a secret-key complies to this
// interface.
type PRNGEncryptor interface {
	Encryptor
	WithPRNG(prng utils.PRNG) PRNGEncryptor
}

type encryptorBase struct {
	params Parameters
	*encryptorBuffers

	prng            utils.PRNG
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	basisextender   *ring.BasisExtender
}

type pkEncryptor struct {
	*encryptorBase
	pk *PublicKey
}

type skEncryptor struct {
	encryptorBase
	sk *SecretKey

	uniformSampler ringqp.UniformSampler
}

// NewEncryptor creates a new Encryptor
// Accepts either a secret-key or a public-key.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	switch key := key.(type) {
	case *PublicKey, PublicKey:
		return newPkEncryptor(params, key)
	case *SecretKey, SecretKey:
		return newSkEncryptor(params, key)
	default:
		panic("cannot NewEncryptor: key must be either *rlwe.PublicKey or *rlwe.SecretKey")
	}
}

// NewPRNGEncryptor creates a new PRNGEncryptor instance.
func NewPRNGEncryptor(params Parameters, key *SecretKey) PRNGEncryptor {
	return newSkEncryptor(params, key)
}

func newEncryptorBase(params Parameters) *encryptorBase {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var bc *ring.BasisExtender
	if params.PCount() != 0 {
		bc = ring.NewBasisExtender(params.RingQ(), params.RingP())
	}

	return &encryptorBase{
		params:           params,
		prng:             prng,
		gaussianSampler:  ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
		ternarySampler:   ring.NewTernarySamplerWithHammingWeight(prng, params.ringQ, params.h, false),
		encryptorBuffers: newEncryptorBuffers(params),
		basisextender:    bc,
	}
}

func newSkEncryptor(params Parameters, key interface{}) (enc *skEncryptor) {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(fmt.Errorf("cannot newSkEncryptor: could not create PRNG for symmetric encryptor: %s", err))
	}

	enc = &skEncryptor{*newEncryptorBase(params), nil, ringqp.NewUniformSampler(prng, *params.RingQP())}
	if enc.sk, err = enc.checkSk(key); err != nil {
		panic(err)
	}

	return enc
}

func newPkEncryptor(params Parameters, key interface{}) (enc *pkEncryptor) {
	var err error
	enc = &pkEncryptor{newEncryptorBase(params), nil}
	enc.pk, err = enc.checkPk(key)
	if err != nil {
		panic(err)
	}
	return enc
}

type encryptorBuffers struct {
	buffQ  [2]*ring.Poly
	buffP  [3]*ring.Poly
	buffQP ringqp.Poly
}

func newEncryptorBuffers(params Parameters) *encryptorBuffers {

	ringQ := params.RingQ()
	ringP := params.RingP()

	var buffP [3]*ring.Poly
	if params.PCount() != 0 {
		buffP = [3]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return &encryptorBuffers{
		buffQ:  [2]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
		buffP:  buffP,
		buffQP: params.RingQP().NewPoly(),
	}
}

// Encrypt encrypts the input plaintext using the stored public-key and writes the result on ct.
// The encryption procedure first samples a new encryption of zero under the public-key and
// then adds the plaintext.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
func (enc *pkEncryptor) Encrypt(pt *Plaintext, ct interface{}) {
	if pt == nil {
		enc.EncryptZero(ct)
	} else {
		switch el := ct.(type) {
		case *Ciphertext:
			if enc.params.PCount() > 0 {
				enc.encrypt(pt, el)
			} else {
				enc.encryptNoP(pt, el)
			}
		default:
			panic("cannot Encrypt: input ciphertext type unsupported (must be *rlwe.Ciphertext or *rgsw.Ciphertext)")
		}
	}
}

// EncryptZero generates an encryption of zero under the stored public-key and writes the result on ct.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
func (enc *pkEncryptor) EncryptZero(ct interface{}) {
	switch ct := ct.(type) {
	case *Ciphertext:
		if enc.params.PCount() > 0 {
			enc.encryptZero(ct)
		} else {
			enc.encryptZeroNoP(ct)
		}
	default:
		panic("cannot EncryptZero: input ciphertext type unsupported")
	}
}

func (enc *pkEncryptor) encrypt(pt *Plaintext, ct *Ciphertext) {
	ringQ := enc.params.RingQ()
	levelQ := utils.MinInt(pt.Level(), ct.Level())
	c0, c1 := ct.Value[0], ct.Value[1]
	ctNTT := ct.Value[0].IsNTT // Store the input NTT flag
	ptNTT := pt.Value.IsNTT

	buffQ0 := enc.buffQ[0]

	// encryptZero checks the NTT flag, but
	// we want the result to always be outside
	// of the NTT domain, regardless of the NTT
	// flag to be able to save one NTT after this
	// step.
	ct.Value[0].IsNTT = false
	enc.encryptZero(ct)

	switch {
	case ctNTT && ptNTT:
		ringQ.NTTLvl(levelQ, c0, c0)
		ringQ.NTTLvl(levelQ, c1, c1)
		ringQ.AddLvl(levelQ, c0, pt.Value, c0)
	case ctNTT && !ptNTT:
		ringQ.AddLvl(levelQ, c0, pt.Value, c0)
		ringQ.NTTLvl(levelQ, c0, c0)
		ringQ.NTTLvl(levelQ, c1, c1)
	case !ctNTT && ptNTT:
		ringQ.InvNTTLvl(levelQ, pt.Value, buffQ0)
		ringQ.AddLvl(levelQ, c0, buffQ0, c0)
	case !ctNTT && !ptNTT:
		ringQ.AddLvl(levelQ, c0, pt.Value, c0)
	}

	// Sets the correct NTT flat back
	c0.IsNTT = ctNTT
	c1.IsNTT = ctNTT
	ct.Resize(ct.Degree(), levelQ)
}

func (enc *pkEncryptor) encryptNoP(pt *Plaintext, ct *Ciphertext) {

	ringQ := enc.params.RingQ()
	levelQ := utils.MinInt(pt.Level(), ct.Level())

	ctNTT := ct.Value[0].IsNTT

	buffQ0 := enc.buffQ[0]

	enc.ternarySampler.ReadLvl(levelQ, buffQ0)
	ringQ.NTTLvl(levelQ, buffQ0, buffQ0)

	c0, c1 := ct.Value[0], ct.Value[1]

	// ct0 = NTT(u*pk0)
	ringQ.MulCoeffsMontgomeryLvl(levelQ, buffQ0, enc.pk.Value[0].Q, c0)
	// ct1 = NTT(u*pk1)
	ringQ.MulCoeffsMontgomeryLvl(levelQ, buffQ0, enc.pk.Value[1].Q, c1)

	// c0
	enc.addPtErrorC0(pt, c0)

	// c1
	if ctNTT {
		enc.gaussianSampler.ReadLvl(levelQ, buffQ0)
		ringQ.NTTLvl(levelQ, buffQ0, buffQ0)
		ringQ.AddLvl(levelQ, c1, buffQ0, c1) // ct1 = u*pk1 + e1

	} else {
		ringQ.InvNTTLvl(levelQ, c1, c1)
		// ct[1] = pk[1]*u + e1
		enc.gaussianSampler.ReadAndAddLvl(ct.Level(), c1)
	}

	c1.IsNTT = ct.Value[0].IsNTT
	ct.Resize(ct.Degree(), levelQ)
}

func (enc *pkEncryptor) encryptZero(ct *Ciphertext) {
	ringQP := enc.params.RingQP()
	levelQ := ct.Level()
	levelP := 0

	buffQ0 := enc.buffQ[0]
	buffP0 := enc.buffP[0]
	buffP1 := enc.buffP[1]
	buffP2 := enc.buffP[2]

	u := ringqp.Poly{Q: buffQ0, P: buffP2}

	// We sample a R-WLE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)
	enc.ternarySampler.ReadLvl(levelQ, u.Q)
	ringQP.ExtendBasisSmallNormAndCenter(u.Q, levelP, nil, u.P)

	// (#Q + #P) NTT
	ringQP.NTTLvl(levelQ, levelP, u, u)

	ct0QP := ringqp.Poly{Q: ct.Value[0], P: buffP0}
	ct1QP := ringqp.Poly{Q: ct.Value[1], P: buffP1}

	// ct0 = u*pk0
	// ct1 = u*pk1
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, u, enc.pk.Value[0], ct0QP)
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, u, enc.pk.Value[1], ct1QP)

	// 2*(#Q + #P) NTT
	ringQP.InvNTTLvl(levelQ, levelP, ct0QP, ct0QP)
	ringQP.InvNTTLvl(levelQ, levelP, ct1QP, ct1QP)

	e := ringqp.Poly{Q: buffQ0, P: buffP2}

	enc.gaussianSampler.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.AddLvl(levelQ, levelP, ct0QP, e, ct0QP)

	enc.gaussianSampler.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.AddLvl(levelQ, levelP, ct1QP, e, ct1QP)

	// ct0 = (u*pk0 + e0)/P
	enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct0QP.Q, ct0QP.P, ct.Value[0])

	// ct1 = (u*pk1 + e1)/P
	enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct1QP.Q, ct1QP.P, ct.Value[1])

	if ct.Value[0].IsNTT {
		ringQP.RingQ.NTTLvl(levelQ, ct.Value[0], ct.Value[0])
		ringQP.RingQ.NTTLvl(levelQ, ct.Value[1], ct.Value[1])
	}
}

func (enc *pkEncryptor) encryptZeroNoP(ct *Ciphertext) {

	ringQ := enc.params.RingQ()
	levelQ := ct.Level()
	ctNTT := ct.Value[0].IsNTT
	buffQ0 := enc.buffQ[0]

	enc.ternarySampler.ReadLvl(levelQ, buffQ0)
	ringQ.NTTLvl(levelQ, buffQ0, buffQ0)

	c0, c1 := ct.Value[0], ct.Value[1]

	// ct0 = NTT(u*pk0)
	ringQ.MulCoeffsMontgomeryLvl(levelQ, buffQ0, enc.pk.Value[0].Q, c0)
	// ct1 = NTT(u*pk1)
	ringQ.MulCoeffsMontgomeryLvl(levelQ, buffQ0, enc.pk.Value[1].Q, c1)

	// c0
	if ctNTT {
		enc.gaussianSampler.ReadLvl(levelQ, buffQ0)
		ringQ.NTTLvl(levelQ, buffQ0, buffQ0)
		ringQ.AddLvl(levelQ, c0, buffQ0, c0)
	} else {
		ringQ.InvNTTLvl(levelQ, c0, c0)
		enc.gaussianSampler.ReadAndAddLvl(levelQ, c0)
	}

	// c1
	if ctNTT {
		enc.gaussianSampler.ReadLvl(levelQ, buffQ0)
		ringQ.NTTLvl(levelQ, buffQ0, buffQ0)
		ringQ.AddLvl(levelQ, c1, buffQ0, c1)

	} else {
		ringQ.InvNTTLvl(levelQ, c1, c1)
		enc.gaussianSampler.ReadAndAddLvl(levelQ, c1)
	}

	c1.IsNTT = ct.Value[0].IsNTT
	ct.Resize(levelQ, levelQ)
}

// Encrypt encrypts the input plaintext using the stored secret-key and writes the result on ct.
// The method accepts only *rlwe.Ciphertext, *CiphertextC0 or *rgsw.Ciphertext as input and will panic otherwise.
func (enc *skEncryptor) Encrypt(pt *Plaintext, ct interface{}) {
	if pt == nil {
		enc.EncryptZero(ct)
	} else {
		switch el := ct.(type) {
		case *Ciphertext:
			enc.uniformSampler.ReadLvl(utils.MinInt(pt.Level(), el.Level()), -1, ringqp.Poly{Q: el.Value[1]})
			enc.encryptRLWE(pt, el.Value[0], el.Value[1])
		case *ring.Poly:
			enc.uniformSampler.ReadLvl(utils.MinInt(pt.Level(), el.Level()), -1, ringqp.Poly{Q: enc.buffQ[1]})
			enc.encryptRLWE(pt, el, enc.buffQ[1])
		default:
			panic("cannot Encrypt: input ciphertext type unsuported (must be *rlwe.Ciphertext or *rgsw.Ciphertext)")
		}
	}

}

// EncryptZero generates an encryption of zero using the stored secret-key and writes the result on ct.
// The method accepts only *rlwe.Ciphertext, *CiphertextC0 or *rgsw.Ciphertext as input and will panic otherwise.
func (enc *skEncryptor) EncryptZero(ct interface{}) {
	switch el := ct.(type) {
	case *Ciphertext:
		enc.uniformSampler.ReadLvl(el.Level(), -1, ringqp.Poly{Q: el.Value[1]})
		enc.encryptZero(el.Value[0], el.Value[1])
	case *ring.Poly:
		enc.uniformSampler.ReadLvl(el.Level(), -1, ringqp.Poly{Q: enc.buffQ[1]})
		enc.encryptZero(el, enc.buffQ[1])
	case *CiphertextQP:
		enc.encryptZeroQP(*el)
	default:
		panic("cannot EncryptZero: input ciphertext type unsupported")
	}
}

func (enc *skEncryptor) encryptZero(c0, c1 *ring.Poly) {

	ringQ := enc.params.RingQ()
	levelQ := c0.Level()

	ringQ.MulCoeffsMontgomeryLvl(levelQ, c1, enc.sk.Value.Q, c0) // c0 = NTT(sc1)
	ringQ.NegLvl(levelQ, c0, c0)                                 // c0 = NTT(-sc1)

	if c0.IsNTT {
		enc.gaussianSampler.ReadLvl(levelQ, enc.buffQ[0]) // e
		ringQ.NTTLvl(levelQ, enc.buffQ[0], enc.buffQ[0])  // NTT(e)
		ringQ.AddLvl(levelQ, c0, enc.buffQ[0], c0)        // c0 = NTT(-sc1 + e)
	} else {
		ringQ.InvNTTLvl(levelQ, c0, c0)               // c0 = -sc1
		ringQ.InvNTTLvl(levelQ, c1, c1)               // c1 = c1
		enc.gaussianSampler.ReadAndAddLvl(levelQ, c0) // c0 = -sc1 + e
	}
}

// EncryptZeroSeeded generates en encryption of zero under sk.
// levelQ : level of the modulus Q
// levelP : level of the modulus P
// sk     : secret key
// sampler: uniform sampler; if `sampler` is nil, then will sample using the internal sampler.
// montgomery: returns the result in the Montgomery domain.
func (enc *skEncryptor) encryptZeroQP(ct CiphertextQP) {

	c0, c1 := ct.Value[0], ct.Value[1]

	levelQ, levelP := c0.LevelQ(), c1.LevelP()
	ringQP := enc.params.RingQP()

	// ct = (e, 0)
	enc.gaussianSampler.ReadLvl(levelQ, c0.Q)
	if levelP != -1 {
		ringQP.ExtendBasisSmallNormAndCenter(c0.Q, levelP, nil, c0.P)
	}

	ringQP.NTTLvl(levelQ, levelP, c0, c0)
	// ct[1] is assumed to be sampled in of the Montgomery domain,
	// thus -as will also be in the Montgomery domain (s is by default), therefore 'e'
	// must be switched to the Montgomery domain.
	ringQP.MFormLvl(levelQ, levelP, c0, c0)

	// ct = (e, a)
	enc.uniformSampler.ReadLvl(levelQ, levelP, c1)

	// (-a*sk + e, a)
	ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, c1, enc.sk.Value, c0)

	if !ct.Value[0].Q.IsNTT {

		ringQP.InvNTTLvl(levelQ, levelP, c0, c0)
		ringQP.InvNTTLvl(levelQ, levelP, c1, c1)
	}
}

func (enc *skEncryptor) encryptRLWE(pt *Plaintext, c0, c1 *ring.Poly) {

	ringQ := enc.params.RingQ()
	levelQ := utils.MinInt(pt.Level(), c0.Level())
	ctNTT := c0.IsNTT

	ringQ.MulCoeffsMontgomeryLvl(levelQ, c1, enc.sk.Value.Q, c0)
	ringQ.NegLvl(levelQ, c0, c0)

	enc.addPtErrorC0(pt, c0)

	if !ctNTT && c1 != enc.buffQ[0] {
		ringQ.InvNTTLvl(levelQ, c1, c1)
	}
}

func (enc *encryptorBase) addPtErrorC0(pt *Plaintext, c0 *ring.Poly) {
	ringQ := enc.params.RingQ()
	levelQ := utils.MinInt(pt.Level(), c0.Level())
	ctNTT := c0.IsNTT
	ptNTT := pt.Value.IsNTT
	buffQ0 := enc.buffQ[0]

	switch {
	case ctNTT && ptNTT:
		enc.gaussianSampler.ReadLvl(levelQ, buffQ0)
		ringQ.NTTLvl(levelQ, buffQ0, buffQ0)
		ringQ.AddLvl(levelQ, c0, buffQ0, c0)
		ringQ.AddLvl(levelQ, c0, pt.Value, c0)
	case ctNTT && !ptNTT:
		enc.gaussianSampler.ReadLvl(levelQ, buffQ0)
		ringQ.AddLvl(levelQ, buffQ0, pt.Value, buffQ0)
		ringQ.NTTLvl(levelQ, buffQ0, buffQ0)
		ringQ.AddLvl(levelQ, c0, buffQ0, c0)
	case !ctNTT && ptNTT:
		ringQ.AddLvl(levelQ, c0, pt.Value, c0)
		ringQ.InvNTTLvl(levelQ, c0, c0)
		enc.gaussianSampler.ReadAndAddLvl(c0.Level(), c0)
	case !ctNTT && !ptNTT:
		ringQ.InvNTTLvl(levelQ, c0, c0)
		ringQ.AddLvl(levelQ, c0, pt.Value, c0)
		enc.gaussianSampler.ReadAndAddLvl(c0.Level(), c0)
	}
}

// ShallowCopy creates a shallow copy of this skEncryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *pkEncryptor) ShallowCopy() Encryptor {
	return NewEncryptor(enc.params, enc.pk)
}

// ShallowCopy creates a shallow copy of this skEncryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *skEncryptor) ShallowCopy() Encryptor {
	return NewEncryptor(enc.params, enc.sk)
}

// WithKey returns this encryptor with a new key.
func (enc *skEncryptor) WithKey(key interface{}) Encryptor {
	skPtr, err := enc.checkSk(key)
	if err != nil {
		panic(err)
	}
	return &skEncryptor{enc.encryptorBase, skPtr, enc.uniformSampler}
}

// WithKey returns this encryptor with a new key.
func (enc *pkEncryptor) WithKey(key interface{}) Encryptor {
	pkPtr, err := enc.checkPk(key)
	if err != nil {
		panic(err)
	}
	return &pkEncryptor{enc.encryptorBase, pkPtr}
}

// WithPRNG returns this encrpytor with prng as its source of randomness for the uniform
// element c1.
func (enc skEncryptor) WithPRNG(prng utils.PRNG) PRNGEncryptor {
	return &skEncryptor{enc.encryptorBase, enc.sk, enc.uniformSampler.WithPRNG(prng)}
}

// checkPk checks that a given pk is correct for the parameters.
func (enc encryptorBase) checkPk(key interface{}) (pk *PublicKey, err error) {

	switch key := key.(type) {
	case PublicKey:
		pk = &key
	case *PublicKey:
		pk = key
	default:
		return nil, fmt.Errorf("key is not a valid public key type %T", key)
	}

	if pk.Value[0].Q.N() != enc.params.N() || pk.Value[1].Q.N() != enc.params.N() {
		return nil, fmt.Errorf("pk ring degree does not match params ring degree")
	}

	return pk, nil
}

// checkPk checks that a given pk is correct for the parameters.
func (enc encryptorBase) checkSk(key interface{}) (sk *SecretKey, err error) {

	switch key := key.(type) {
	case SecretKey:
		sk = &key
	case *SecretKey:
		sk = key
	default:
		return nil, fmt.Errorf("key is not a valid public key type %T", key)
	}

	if sk.Value.Q.N() != enc.params.N() {
		panic("cannot checkSk: sk ring degree does not match params ring degree")
	}

	return sk, nil
}
