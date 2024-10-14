package rlwe

import (
	"fmt"
	"reflect"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// EncryptionKey is an interface for encryption keys. Valid encryption
// keys are the [SecretKey] and [PublicKey] types.
type EncryptionKey interface {
	isEncryptionKey()
}

// NewEncryptor creates a new [Encryptor] from either a public key or a private key.
func NewEncryptor(params ParameterProvider, key EncryptionKey) *Encryptor {

	p := *params.GetRLWEParameters()

	enc := newEncryptor(p)
	var err error
	switch key := key.(type) {
	case *PublicKey:
		err = enc.checkPk(key)
	case *SecretKey:
		err = enc.checkSk(key)
	case nil:
		return newEncryptor(p)
	default:
		// Sanity check
		panic(fmt.Errorf("key must be either *rlwe.PublicKey, *rlwe.SecretKey or nil but have %T", key))
	}

	if err != nil {
		// Sanity check, this error should not happen.
		panic(fmt.Errorf("key is not correct: %w", err))
	}

	enc.encKey = key
	return enc
}

type Encryptor struct {
	params Parameters
	*encryptorBuffers

	encKey         EncryptionKey
	prng           sampling.PRNG
	xeSampler      ring.Sampler
	xsSampler      ring.Sampler
	basisextender  *ring.BasisExtender
	uniformSampler ringqp.UniformSampler
}

// GetRLWEParameters returns the underlying [Parameters].
func (enc Encryptor) GetRLWEParameters() *Parameters {
	return &enc.params
}

func newEncryptor(params Parameters) *Encryptor {

	prng, err := sampling.NewPRNG()
	if err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}

	var bc *ring.BasisExtender
	if params.PCount() != 0 {
		bc = ring.NewBasisExtender(params.RingQ(), params.RingP())
	}

	xeSampler, err := ring.NewSampler(prng, params.RingQ(), params.Xe(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(fmt.Errorf("newEncryptor: %w", err))
	}

	xsSampler, err := ring.NewSampler(prng, params.RingQ(), params.Xs(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(fmt.Errorf("newEncryptor: %w", err))
	}

	return &Encryptor{
		params:           params,
		prng:             prng,
		xeSampler:        xeSampler,
		xsSampler:        xsSampler,
		encryptorBuffers: newEncryptorBuffers(params),
		uniformSampler:   ringqp.NewUniformSampler(prng, *params.RingQP()),
		basisextender:    bc,
	}
}

type encryptorBuffers struct {
	buffQP [3]ringqp.Poly
}

func newEncryptorBuffers(params Parameters) *encryptorBuffers {
	return &encryptorBuffers{
		buffQP: [3]ringqp.Poly{
			params.RingQP().NewPoly(),
			params.RingQP().NewPoly(),
			params.RingQP().NewPoly(),
		},
	}
}

// Encrypt encrypts the input plaintext using the stored encryption key and writes the result on ct.
// The method currently accepts only *[Ciphertext] as ct.
// If a [Plaintext] is given, then the output [Ciphertext] [MetaData] will match the [Plaintext] [MetaData].
// The method returns an error if the ct has an unsupported type or if no encryption key is stored
// in the [Encryptor].
//
// The encryption procedure masks the plaintext by adding a fresh encryption of zero.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
func (enc Encryptor) Encrypt(pt *Plaintext, ct interface{}) (err error) {
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

// EncryptNew encrypts the input plaintext using the stored encryption key and returns a newly
// allocated [Ciphertext] containing the result.
// If a [Plaintext] is provided, then the output [Ciphertext] [MetaData] will match the [Plaintext] [MetaData].
// The method returns an error if the ct has an unsupported type or if no encryption key is stored
// in the [Encryptor].
//
// The encryption procedure masks the plaintext by adding a fresh encryption of zero.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
func (enc Encryptor) EncryptNew(pt *Plaintext) (ct *Ciphertext, err error) {
	ct = NewCiphertext(enc.params, 1, pt.Level())
	return ct, enc.Encrypt(pt, ct)
}

// EncryptZero generates an encryption of zero under the stored encryption key and writes the result on ct.
// The method accepts only *[Ciphertext] as input.
// The method returns an error if the ct has an unsupported type or if no encryption key is stored
// in the [Encryptor].
//
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The zero encryption is generated according to the given [Ciphertext] [MetaData].
func (enc Encryptor) EncryptZero(ct interface{}) (err error) {
	switch key := enc.encKey.(type) {
	case *SecretKey:
		return enc.encryptZeroSk(key, ct)
	case *PublicKey:
		if cti, isCt := ct.(*Ciphertext); isCt && enc.params.PCount() == 0 {
			return enc.encryptZeroPkNoP(key, cti.Element)
		}
		return enc.encryptZeroPk(key, ct)
	default:
		return fmt.Errorf("cannot encrypt: Encryptor has no encryption key")
	}
}

// EncryptZeroNew generates an encryption of zero under the stored encryption key and returns a newly
// allocated [Ciphertext] containing the result.
// The method returns an error if no encryption key is stored in the [Encryptor].
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
func (enc Encryptor) EncryptZeroNew(level int) (ct *Ciphertext) {
	ct = NewCiphertext(enc.params, 1, level)
	if err := enc.EncryptZero(ct); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}
	return
}

func (enc Encryptor) encryptZeroPk(pk *PublicKey, ct interface{}) (err error) {

	var ct0QP, ct1QP ringqp.Poly

	if ctCt, isCiphertext := ct.(*Ciphertext); isCiphertext {
		ct = ctCt.Element
	}

	var levelQ, levelP int
	switch ct := ct.(type) {
	case Element[ring.Poly]:

		levelQ = ct.Level()
		levelP = 0

		ct0QP = ringqp.Poly{Q: ct.Value[0], P: enc.buffQP[0].Q}
		ct1QP = ringqp.Poly{Q: ct.Value[1], P: enc.buffQP[0].P}
	case Element[ringqp.Poly]:

		levelQ = ct.LevelQ()
		levelP = ct.LevelP()

		ct0QP = ct.Value[0]
		ct1QP = ct.Value[1]
	default:
		return fmt.Errorf("invalid input: must be Element[ring.Poly] or Element[ringqp.Poly] but is %T", ct)
	}

	ringQP := enc.params.RingQP().AtLevel(levelQ, levelP)

	u := enc.buffQP[1]

	// We sample a RLWE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)
	enc.xsSampler.AtLevel(levelQ).Read(u.Q)
	ringQP.ExtendBasisSmallNormAndCenter(u.Q, levelP, u.Q, u.P)

	// (#Q + #P) NTT
	ringQP.NTT(u, u)

	// ct0 = u*pk0
	// ct1 = u*pk1
	ringQP.MulCoeffsMontgomery(u, pk.Value[0], ct0QP)
	ringQP.MulCoeffsMontgomery(u, pk.Value[1], ct1QP)

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
	case Element[ring.Poly]:

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

	case Element[ringqp.Poly]:
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

func (enc Encryptor) encryptZeroPkNoP(pk *PublicKey, ct Element[ring.Poly]) (err error) {

	levelQ := ct.Level()

	ringQ := enc.params.RingQ().AtLevel(levelQ)

	buffQ0 := enc.buffQP[0].Q

	enc.xsSampler.AtLevel(levelQ).Read(buffQ0)
	ringQ.NTT(buffQ0, buffQ0)

	c0, c1 := ct.Value[0], ct.Value[1]

	// ct0 = NTT(u*pk0)
	ringQ.MulCoeffsMontgomery(buffQ0, pk.Value[0].Q, c0)
	// ct1 = NTT(u*pk1)
	ringQ.MulCoeffsMontgomery(buffQ0, pk.Value[1].Q, c1)

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

// encryptZeroSk generates an encryption of zero using the stored secret-key and writes the result on ct.
// The method accepts only *rlwe.Ciphertext or *rgsw.Ciphertext as input and will return an error otherwise.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc Encryptor) encryptZeroSk(sk *SecretKey, ct interface{}) (err error) {

	switch ct := ct.(type) {
	case *Ciphertext:

		var c1 ring.Poly
		if ct.Degree() == 1 {
			c1 = ct.Value[1]
		} else {
			c1 = enc.buffQP[1].Q
		}

		enc.uniformSampler.AtLevel(ct.Level(), -1).Read(ringqp.Poly{Q: c1})

		if !ct.IsNTT {
			enc.params.RingQ().AtLevel(ct.Level()).NTT(c1, c1)
		}

		return enc.encryptZeroSkFromC1(sk, ct.Element, c1)

	case Element[ringqp.Poly]:

		var c1 ringqp.Poly

		if ct.Degree() == 1 {
			c1 = ct.Value[1]
		} else {
			c1 = enc.buffQP[1]
		}

		// ct = (e, a)
		enc.uniformSampler.AtLevel(ct.LevelQ(), ct.LevelP()).Read(c1)

		if !ct.IsNTT {
			enc.params.RingQP().AtLevel(ct.LevelQ(), ct.LevelP()).NTT(c1, c1)
		}

		return enc.encryptZeroSkFromC1QP(sk, ct, c1)
	default:
		return fmt.Errorf("cannot EncryptZero: input ciphertext type %T is not supported", ct)
	}
}

func (enc Encryptor) encryptZeroSkFromC1(sk *SecretKey, ct Element[ring.Poly], c1 ring.Poly) (err error) {

	levelQ := ct.Level()

	ringQ := enc.params.RingQ().AtLevel(levelQ)

	c0 := ct.Value[0]

	ringQ.MulCoeffsMontgomery(c1, sk.Value.Q, c0)
	ringQ.Neg(c0, c0)

	if ct.IsNTT {
		e := enc.buffQP[0].Q
		enc.xeSampler.AtLevel(levelQ).Read(e)
		ringQ.NTT(e, e)
		ringQ.Add(c0, e, c0)
	} else {
		ringQ.INTT(c0, c0)
		if ct.Degree() == 1 {
			ringQ.INTT(c1, c1)
		}

		enc.xeSampler.AtLevel(levelQ).ReadAndAdd(c0)
	}

	return
}

// EncryptZeroSeeded generates en encryption of zero under sk.
// levelQ : level of the modulus Q
// levelP : level of the modulus P
// sk     : secret key
// sampler: uniform sampler; if `sampler` is nil, then the internal sampler will be used.
// montgomery: returns the result in the Montgomery domain.
func (enc Encryptor) encryptZeroSkFromC1QP(sk *SecretKey, ct Element[ringqp.Poly], c1 ringqp.Poly) (err error) {

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
	ringQP.MulCoeffsMontgomeryThenSub(c1, sk.Value, c0)

	if !ct.IsNTT {
		ringQP.INTT(c0, c0)
		ringQP.INTT(c1, c1)
	}

	return
}

// WithPRNG returns this encryptor with prng as its source of randomness for the uniform
// element c1.
// The returned encryptor isn't safe to use concurrently with the original encryptor.
func (enc Encryptor) WithPRNG(prng sampling.PRNG) *Encryptor {
	enc.uniformSampler = ringqp.NewUniformSampler(prng, *enc.params.RingQP())
	return &enc
}

func (enc Encryptor) ShallowCopy() *Encryptor {
	return NewEncryptor(enc.params, enc.encKey)
}

func (enc Encryptor) WithKey(key EncryptionKey) *Encryptor {
	switch key := key.(type) {
	case *SecretKey:
		if err := enc.checkSk(key); err != nil {
			// Sanity check, this error should not happen.
			panic(fmt.Errorf("cannot WithKey: %w", err))
		}
	case *PublicKey:
		if err := enc.checkPk(key); err != nil {
			// Sanity check, this error should not happen.
			panic(fmt.Errorf("cannot WithKey: %w", err))
		}
	case nil:
		return &enc
	default:
		// Sanity check, this error should not happen.
		panic(fmt.Errorf("invalid key type, want *rlwe.SecretKey, *rlwe.PublicKey or nil but have %T", key))
	}
	enc.encKey = key
	return &enc
}

// checkPk checks that a given pk is correct for the parameters.
func (enc Encryptor) checkPk(pk *PublicKey) (err error) {
	if pk.Value[0].Q.N() != enc.params.N() || pk.Value[1].Q.N() != enc.params.N() {
		return fmt.Errorf("pk ring degree does not match params ring degree")
	}
	return
}

// checkPk checks that a given pk is correct for the parameters.
func (enc Encryptor) checkSk(sk *SecretKey) (err error) {
	if sk.Value.Q.N() != enc.params.N() {
		return fmt.Errorf("sk ring degree does not match params ring degree")
	}
	return
}

func (enc Encryptor) addPtToCt(level int, pt *Plaintext, ct *Ciphertext) {

	ringQ := enc.params.RingQ().AtLevel(level)
	var buff ring.Poly
	if pt.IsNTT {
		if ct.IsNTT {
			buff = pt.Value
		} else {
			buff = enc.buffQP[0].Q
			ringQ.NTT(pt.Value, buff)
		}
	} else {
		if ct.IsNTT {
			buff = enc.buffQP[0].Q
			ringQ.INTT(pt.Value, buff)
		} else {
			buff = pt.Value
		}
	}

	ringQ.Add(ct.Value[0], buff, ct.Value[0])
}
