package rlwe

import (
	"fmt"
	"reflect"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// EncryptionKey is an interface for encryption keys. Valid encryption
// keys are the SecretKey and PublicKey types.
type EncryptionKey interface {
	isEncryptionKey()
}

// NewEncryptor creates a new Encryptor.
// Accepts either a secret-key or a public-key.
func NewEncryptor(params ParametersInterface, key EncryptionKey) (EncryptorInterface, error) {
	switch key := key.(type) {
	case *PublicKey:
		return NewEncryptorPublicKey(params, key)
	case *SecretKey:
		return NewEncryptorSecretKey(params, key)
	case nil:
		return newEncryptorBase(params), nil
	default:
		return nil, fmt.Errorf("cannot NewEncryptor: key must be either *rlwe.PublicKey, *rlwe.SecretKey or nil but have %T", key)
	}
}

// NewPRNGEncryptor creates a new PRNGEncryptor instance.
func NewPRNGEncryptor(params ParametersInterface, key *SecretKey) (PRNGEncryptorInterface, error) {
	return NewEncryptorSecretKey(params, key)
}

type encryptorBase struct {
	params ParametersInterface
	*encryptorBuffers

	prng           sampling.PRNG
	xeSampler      ring.Sampler
	xsSampler      ring.Sampler
	basisextender  *ring.BasisExtender
	uniformSampler ringqp.UniformSampler
}

func newEncryptorBase(params ParametersInterface) *encryptorBase {

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	var bc *ring.BasisExtender
	if params.PCount() != 0 {
		bc = ring.NewBasisExtender(params.RingQ(), params.RingP())
	}

	xeSampler, err := ring.NewSampler(prng, params.RingQ(), params.Xe(), false)

	if err != nil {
		panic(fmt.Errorf("newEncryptorBase: %w", err))
	}

	xsSampler, err := ring.NewSampler(prng, params.RingQ(), params.Xs(), false)

	if err != nil {
		panic(fmt.Errorf("newEncryptorBase: %w", err))
	}

	return &encryptorBase{
		params:           params,
		prng:             prng,
		xeSampler:        xeSampler,
		xsSampler:        xsSampler,
		encryptorBuffers: newEncryptorBuffers(params),
		uniformSampler:   ringqp.NewUniformSampler(prng, *params.RingQP()),
		basisextender:    bc,
	}
}

// EncryptorSecretKey is an encryptor using an `rlwe.SecretKey` to encrypt.
type EncryptorSecretKey struct {
	encryptorBase
	sk *SecretKey
}

// NewEncryptorSecretKey creates a new EncryptorSecretKey from the provided parameters and secret key.
func NewEncryptorSecretKey(params ParametersInterface, sk *SecretKey) (enc *EncryptorSecretKey, err error) {

	enc = &EncryptorSecretKey{*newEncryptorBase(params), nil}

	if err = enc.checkSk(sk); err != nil {
		return nil, fmt.Errorf("cannot NewEncryptorSecretKey: %w", err)
	}

	enc.sk = sk

	return
}

// EncryptorPublicKey is an encryptor using an `rlwe.PublicKey` to encrypt.
type EncryptorPublicKey struct {
	encryptorBase
	pk *PublicKey
}

// NewEncryptorPublicKey creates a new EncryptorPublicKey from the provided parameters and secret key.
func NewEncryptorPublicKey(params ParametersInterface, pk *PublicKey) (enc *EncryptorPublicKey, err error) {

	enc = &EncryptorPublicKey{*newEncryptorBase(params), nil}

	if err = enc.checkPk(pk); err != nil {
		return nil, fmt.Errorf("cannot NewEncryptorPublicKey: %w", err)
	}

	enc.pk = pk

	return
}

type encryptorBuffers struct {
	buffQ  [2]ring.Poly
	buffP  [3]ring.Poly
	buffQP ringqp.Poly
}

func newEncryptorBuffers(params ParametersInterface) *encryptorBuffers {

	ringQ := params.RingQ()
	ringP := params.RingP()

	var buffP [3]ring.Poly
	if params.PCount() != 0 {
		buffP = [3]ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return &encryptorBuffers{
		buffQ:  [2]ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
		buffP:  buffP,
		buffQP: params.RingQP().NewPoly(),
	}
}

// Encrypt encrypts the input plaintext using the stored public-key and writes the result on ct.
// The encryption procedure first samples a new encryption of zero under the public-key and
// then adds the Plaintext.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
// If a Plaintext is given, then the output Ciphertext MetaData will match the Plaintext MetaData.
func (enc EncryptorPublicKey) Encrypt(pt *Plaintext, ct interface{}) (err error) {

	if pt == nil {
		return enc.EncryptZero(ct)
	} else {
		switch ct := ct.(type) {
		case *Ciphertext:

			*ct.MetaData = *pt.MetaData

			level := utils.Min(pt.Level(), ct.Level())

			ct.Resize(ct.Degree(), level)

			if err = enc.EncryptZero(ct); err != nil {
				return fmt.Errorf("cannot Encrypt: %w", err)
			}

			enc.addPtToCt(level, pt, ct)

			return

		default:
			return fmt.Errorf("cannot Encrypt: input ciphertext type %s is not supported", reflect.TypeOf(ct))
		}
	}
}

// EncryptNew encrypts the input plaintext using the stored public-key and returns the result on a new Ciphertext.
// The encryption procedure first samples a new encryption of zero under the public-key and
// then adds the Plaintext.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// If a Plaintext is given, then the output ciphertext MetaData will match the Plaintext MetaData.
func (enc EncryptorPublicKey) EncryptNew(pt *Plaintext) (ct *Ciphertext, err error) {
	ct = NewCiphertext(enc.params, 1, pt.Level())
	return ct, enc.Encrypt(pt, ct)
}

// EncryptZeroNew generates an encryption of zero under the stored public-key and returns it on a new Ciphertext.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc EncryptorPublicKey) EncryptZeroNew(level int) (ct *Ciphertext) {
	ct = NewCiphertext(enc.params, 1, level)
	if err := enc.EncryptZero(ct); err != nil {
		panic(err)
	}
	return
}

// EncryptZero generates an encryption of zero under the stored public-key and writes the result on ct.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc EncryptorPublicKey) EncryptZero(ct interface{}) (err error) {
	switch ct := ct.(type) {
	case *Ciphertext:
		if enc.params.PCount() > 0 {
			return enc.encryptZero(*ct.El())
		} else {
			return enc.encryptZeroNoP(ct)
		}
	case Operand[ringqp.Poly]:
		return enc.encryptZero(ct)
	default:
		return fmt.Errorf("cannot Encrypt: input ciphertext type %s is not supported", reflect.TypeOf(ct))
	}
}

func (enc EncryptorPublicKey) encryptZero(ct interface{}) (err error) {

	var ct0QP, ct1QP ringqp.Poly

	var levelQ, levelP int
	switch ct := ct.(type) {
	case Operand[ring.Poly]:

		levelQ = ct.Level()
		levelP = 0

		ct0QP = ringqp.Poly{Q: ct.Value[0], P: enc.buffP[0]}
		ct1QP = ringqp.Poly{Q: ct.Value[1], P: enc.buffP[1]}
	case Operand[ringqp.Poly]:

		levelQ = ct.LevelQ()
		levelP = ct.LevelP()

		ct0QP = ct.Value[0]
		ct1QP = ct.Value[1]
	default:
		return fmt.Errorf("invalid input: must be OperandQ or OperandQP but is %T", ct)
	}

	ringQP := enc.params.RingQP().AtLevel(levelQ, levelP)

	u := ringqp.Poly{Q: enc.buffQ[0], P: enc.buffP[2]}

	// We sample a RLWE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)
	enc.xsSampler.AtLevel(levelQ).Read(u.Q)
	ringQP.ExtendBasisSmallNormAndCenter(u.Q, levelP, u.Q, u.P)

	// (#Q + #P) NTT
	ringQP.NTT(u, u)

	// ct0 = u*pk0
	// ct1 = u*pk1
	ringQP.MulCoeffsMontgomery(u, enc.pk.Value[0], ct0QP)
	ringQP.MulCoeffsMontgomery(u, enc.pk.Value[1], ct1QP)

	// 2*(#Q + #P) NTT
	ringQP.INTT(ct0QP, ct0QP)
	ringQP.INTT(ct1QP, ct1QP)

	e := u

	enc.xeSampler.AtLevel(levelQ).Read(e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, e.Q, e.P)
	ringQP.Add(ct0QP, e, ct0QP)

	enc.xeSampler.AtLevel(levelQ).Read(e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, e.Q, e.P)
	ringQP.Add(ct1QP, e, ct1QP)

	switch ct := ct.(type) {
	case Operand[ring.Poly]:

		// ct0 = (u*pk0 + e0)/P
		enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct0QP.Q, ct0QP.P, ct.Value[0])

		// ct1 = (u*pk1 + e1)/P
		enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct1QP.Q, ct1QP.P, ct.Value[1])

		if ct.IsNTT {
			ringQP.RingQ.NTT(ct.Value[0], ct.Value[0])
			ringQP.RingQ.NTT(ct.Value[1], ct.Value[1])
		}

		if ct.IsMontgomery {
			ringQP.RingQ.MForm(ct.Value[0], ct.Value[0])
			ringQP.RingQ.MForm(ct.Value[1], ct.Value[1])
		}

	case Operand[ringqp.Poly]:
		if ct.IsNTT {
			ringQP.NTT(ct.Value[0], ct.Value[0])
			ringQP.NTT(ct.Value[1], ct.Value[1])
		}

		if ct.IsMontgomery {
			ringQP.MForm(ct.Value[0], ct.Value[0])
			ringQP.MForm(ct.Value[1], ct.Value[1])
		}
	}

	return
}

func (enc EncryptorPublicKey) encryptZeroNoP(ct *Ciphertext) (err error) {

	levelQ := ct.Level()

	ringQ := enc.params.RingQ().AtLevel(levelQ)

	buffQ0 := enc.buffQ[0]

	enc.xsSampler.AtLevel(levelQ).Read(buffQ0)
	ringQ.NTT(buffQ0, buffQ0)

	c0, c1 := ct.Value[0], ct.Value[1]

	// ct0 = NTT(u*pk0)
	ringQ.MulCoeffsMontgomery(buffQ0, enc.pk.Value[0].Q, c0)
	// ct1 = NTT(u*pk1)
	ringQ.MulCoeffsMontgomery(buffQ0, enc.pk.Value[1].Q, c1)

	// c0
	if ct.IsNTT {
		enc.xeSampler.AtLevel(levelQ).Read(buffQ0)
		ringQ.NTT(buffQ0, buffQ0)
		ringQ.Add(c0, buffQ0, c0)
	} else {
		ringQ.INTT(c0, c0)
		enc.xeSampler.AtLevel(levelQ).ReadAndAdd(c0)
	}

	// c1
	if ct.IsNTT {
		enc.xeSampler.AtLevel(levelQ).Read(buffQ0)
		ringQ.NTT(buffQ0, buffQ0)
		ringQ.Add(c1, buffQ0, c1)

	} else {
		ringQ.INTT(c1, c1)
		enc.xeSampler.AtLevel(levelQ).ReadAndAdd(c1)
	}

	return
}

// Encrypt encrypts the input plaintext using the stored secret-key and writes the result on ct.
// The method accepts only *rlwe.Ciphertext or *rgsw.Ciphertext as input and will return an error otherwise.
// If a plaintext is given, the encryptor only accepts *rlwe.Ciphertext, and the generated Ciphertext
// MetaData will match the given Plaintext MetaData.
func (enc EncryptorSecretKey) Encrypt(pt *Plaintext, ct interface{}) (err error) {
	if pt == nil {
		return enc.EncryptZero(ct)
	} else {
		switch ct := ct.(type) {
		case *Ciphertext:
			*ct.MetaData = *pt.MetaData
			level := utils.Min(pt.Level(), ct.Level())
			ct.Resize(ct.Degree(), level)
			if err = enc.EncryptZero(ct); err != nil {
				return
			}
			enc.addPtToCt(level, pt, ct)
			return
		default:
			return fmt.Errorf("cannot Encrypt: input ciphertext type %T is not supported", ct)
		}
	}
}

// EncryptNew encrypts the input plaintext using the stored secret-key and returns the result on a new Ciphertext.
// MetaData will match the given Plaintext MetaData.
func (enc EncryptorSecretKey) EncryptNew(pt *Plaintext) (ct *Ciphertext, err error) {
	ct = NewCiphertext(enc.params, 1, pt.Level())
	return ct, enc.Encrypt(pt, ct)
}

// EncryptZero generates an encryption of zero using the stored secret-key and writes the result on ct.
// The method accepts only *rlwe.Ciphertext or *rgsw.Ciphertext as input and will return an error otherwise.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc EncryptorSecretKey) EncryptZero(ct interface{}) (err error) {
	switch ct := ct.(type) {
	case *Ciphertext:

		var c1 ring.Poly
		if ct.Degree() == 1 {
			c1 = ct.Value[1]
		} else {
			c1 = enc.buffQ[1]
		}

		enc.uniformSampler.AtLevel(ct.Level(), -1).Read(ringqp.Poly{Q: c1})

		if !ct.IsNTT {
			enc.params.RingQ().AtLevel(ct.Level()).NTT(c1, c1)
		}

		return enc.encryptZero(ct.Operand, c1)

	case Operand[ringqp.Poly]:

		var c1 ringqp.Poly

		if ct.Degree() == 1 {
			c1 = ct.Value[1]
		} else {
			c1 = enc.buffQP
		}

		// ct = (e, a)
		enc.uniformSampler.AtLevel(ct.LevelQ(), ct.LevelP()).Read(c1)

		if !ct.IsNTT {
			enc.params.RingQP().AtLevel(ct.LevelQ(), ct.LevelP()).NTT(c1, c1)
		}

		return enc.encryptZeroQP(ct, c1)

	default:
		return fmt.Errorf("cannot EncryptZero: input ciphertext type %T is not supported", ct)
	}
}

// EncryptZeroNew generates an encryption of zero using the stored secret-key and writes the result on ct.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc EncryptorSecretKey) EncryptZeroNew(level int) (ct *Ciphertext) {
	ct = NewCiphertext(enc.params, 1, level)
	if err := enc.EncryptZero(ct); err != nil {
		panic(err)
	}
	return
}

func (enc EncryptorSecretKey) encryptZero(ct Operand[ring.Poly], c1 ring.Poly) (err error) {

	levelQ := ct.Level()

	ringQ := enc.params.RingQ().AtLevel(levelQ)

	c0 := ct.Value[0]

	ringQ.MulCoeffsMontgomery(c1, enc.sk.Value.Q, c0) // c0 = NTT(sc1)
	ringQ.Neg(c0, c0)                                 // c0 = NTT(-sc1)

	if ct.IsNTT {
		enc.xeSampler.AtLevel(levelQ).Read(enc.buffQ[0]) // e
		ringQ.NTT(enc.buffQ[0], enc.buffQ[0])            // NTT(e)
		ringQ.Add(c0, enc.buffQ[0], c0)                  // c0 = NTT(-sc1 + e)
	} else {
		ringQ.INTT(c0, c0) // c0 = -sc1
		if ct.Degree() == 1 {
			ringQ.INTT(c1, c1) // c1 = c1
		}

		enc.xeSampler.AtLevel(levelQ).ReadAndAdd(c0) // c0 = -sc1 + e
	}

	return
}

// EncryptZeroSeeded generates en encryption of zero under sk.
// levelQ : level of the modulus Q
// levelP : level of the modulus P
// sk     : secret key
// sampler: uniform sampler; if `sampler` is nil, then the internal sampler will be used.
// montgomery: returns the result in the Montgomery domain.
func (enc EncryptorSecretKey) encryptZeroQP(ct Operand[ringqp.Poly], c1 ringqp.Poly) (err error) {

	levelQ, levelP := ct.LevelQ(), ct.LevelP()
	ringQP := enc.params.RingQP().AtLevel(levelQ, levelP)

	c0 := ct.Value[0]

	// ct = (e, 0)
	enc.xeSampler.AtLevel(levelQ).Read(c0.Q)
	if levelP != -1 {
		ringQP.ExtendBasisSmallNormAndCenter(c0.Q, levelP, c0.Q, c0.P)
	}

	ringQP.NTT(c0, c0)
	// ct[1] is assumed to be sampled in of the Montgomery domain,
	// thus -as will also be in the Montgomery domain (s is by default), therefore 'e'
	// must be switched to the Montgomery domain.
	ringQP.MForm(c0, c0)

	// (-a*sk + e, a)
	ringQP.MulCoeffsMontgomeryThenSub(c1, enc.sk.Value, c0)

	if !ct.IsNTT {
		ringQP.INTT(c0, c0)
		ringQP.INTT(c1, c1)
	}

	return
}

// ShallowCopy creates a shallow copy of this EncryptorSecretKey in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc EncryptorPublicKey) ShallowCopy() EncryptorInterface {
	encSh, _ := NewEncryptorPublicKey(enc.params, enc.pk)
	return encSh
}

// ShallowCopy creates a shallow copy of this EncryptorSecretKey in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc EncryptorSecretKey) ShallowCopy() EncryptorInterface {
	encSh, _ := NewEncryptorSecretKey(enc.params, enc.sk)
	return encSh
}

// WithPRNG returns this encryptor with prng as its source of randomness for the uniform
// element c1.
func (enc EncryptorSecretKey) WithPRNG(prng sampling.PRNG) PRNGEncryptorInterface {
	encBase := enc.encryptorBase
	encBase.uniformSampler = ringqp.NewUniformSampler(prng, *enc.params.RingQP())
	return &EncryptorSecretKey{encBase, enc.sk}
}

func (enc encryptorBase) Encrypt(pt *Plaintext, ct interface{}) (err error) {
	return fmt.Errorf("cannot Encrypt: key hasn't been set")
}

func (enc encryptorBase) EncryptNew(pt *Plaintext) (ct *Ciphertext, err error) {
	return nil, fmt.Errorf("cannot EncryptNew: key hasn't been set")
}

func (enc encryptorBase) EncryptZero(ct interface{}) (err error) {
	return fmt.Errorf("cannot EncryptZeroNew: key hasn't been set")
}

func (enc encryptorBase) EncryptZeroNew(level int) (ct *Ciphertext) {
	panic("cannot EncryptZeroNew: key hasn't been set")
}

func (enc encryptorBase) ShallowCopy() EncryptorInterface {
	encSh, _ := NewEncryptor(enc.params, nil)
	return encSh
}

func (enc encryptorBase) WithKey(key interface{}) (EncryptorInterface, error) {
	switch key := key.(type) {
	case *SecretKey:
		if err := enc.checkSk(key); err != nil {
			return nil, fmt.Errorf("cannot WithKey: %w", err)
		}
		return &EncryptorSecretKey{enc, key}, nil
	case *PublicKey:
		if err := enc.checkPk(key); err != nil {
			return nil, fmt.Errorf("cannot WithKey: %w", err)
		}
		return &EncryptorPublicKey{enc, key}, nil
	case nil:
		return &enc, nil
	default:
		return nil, fmt.Errorf("invalid key type, want *rlwe.SecretKey, *rlwe.PublicKey or nil but have %T", key)
	}
}

// checkPk checks that a given pk is correct for the parameters.
func (enc encryptorBase) checkPk(pk *PublicKey) (err error) {
	if pk.Value[0].Q.N() != enc.params.N() || pk.Value[1].Q.N() != enc.params.N() {
		return fmt.Errorf("pk ring degree does not match params ring degree")
	}
	return
}

// checkPk checks that a given pk is correct for the parameters.
func (enc encryptorBase) checkSk(sk *SecretKey) (err error) {
	if sk.Value.Q.N() != enc.params.N() {
		return fmt.Errorf("sk ring degree does not match params ring degree")
	}
	return
}

func (enc encryptorBase) addPtToCt(level int, pt *Plaintext, ct *Ciphertext) {

	ringQ := enc.params.RingQ().AtLevel(level)
	var buff ring.Poly
	if pt.IsNTT {
		if ct.IsNTT {
			buff = pt.Value
		} else {
			buff = enc.buffQ[0]
			ringQ.NTT(pt.Value, buff)
		}
	} else {
		if ct.IsNTT {
			buff = enc.buffQ[0]
			ringQ.INTT(pt.Value, buff)
		} else {
			buff = pt.Value
		}
	}

	ringQ.Add(ct.Value[0], buff, ct.Value[0])
}
