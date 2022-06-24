package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/gadget"
	"github.com/tuneinsight/lattigo/v3/rlwe/rgsw"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Encryptor a generic RLWE encryption interface.
type Encryptor interface {
	Encrypt(pt *Plaintext, ct interface{})
	//EncryptZero(ct interface{})
	EncryptZeroSymetricQPNTT(levelQ, levelP int, sk ringqp.Poly, sampler *ringqp.UniformSampler, montgomery bool, ct interface{})
	ShallowCopy() Encryptor
	WithKey(key interface{}) Encryptor
}

type SymmetricEncryptor interface {
	Encryptor
	EncryptSeeded(pt *Plaintext, sampler ringqp.UniformSampler, ct interface{})
	EncryptZeroSeeded(sampler ringqp.UniformSampler, ct interface{})
}

type encryptor struct {
	*encryptorBase
	*encryptorSamplers
	*encryptorBuffers
	basisextender *ring.BasisExtender
}

type pkEncryptor struct {
	encryptor
	pk *PublicKey
}

type skEncryptor struct {
	encryptor
	sk *SecretKey
}

// NewEncryptor creates a new Encryptor
// Accepts either a secret-key or a public-key.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	enc := newEncryptor(params)
	return enc.setKey(key)
}

func NewEncryptorSymmetric(params Parameters, key *SecretKey) SymmetricEncryptor {
	enc := NewEncryptor(params, key)
	return enc.(SymmetricEncryptor)
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
	prng            utils.PRNG
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSampler  ringqp.UniformSampler
}

func newEncryptorSamplers(params Parameters) *encryptorSamplers {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &encryptorSamplers{
		prng:            prng,
		gaussianSampler: ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
		ternarySampler:  ring.NewTernarySamplerWithHammingWeight(prng, params.ringQ, params.h, false),
		uniformSampler:  ringqp.NewUniformSampler(prng, *params.RingQP()),
	}
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
// The encryption procedure first samples an new encryption of zero under the public-key and
// then adds the plaintext.
// The encryption procedures depends on the parameters. If the auxiliary modulus P is defined,
// then the encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly
// sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
func (enc *pkEncryptor) Encrypt(pt *Plaintext, ct interface{}) {
	switch el := ct.(type) {
	case *Ciphertext:
		if enc.basisextender != nil {
			enc.encryptRLWE(pt, el)
		} else {
			enc.encryptNoPRLWE(pt, el)
		}
	default:
		panic("input ciphertext type unsuported (must be *rlwe.Ciphertext or *rgsw.Ciphertext)")
	}

}

// Encrypt encrypts the input plaintext using the stored public-key and writes the result on ct.
// The encryption procedure first samples an new encryption of zero under the public-key and
// then adds the plaintext.
// The encryption procedures depends on the parameters. If the auxiliary modulus P is defined,
// then the encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly
// sampled in Q.
// The method accepts only *rlwe.Ciphertext or *rgsw.Ciphertext as input and will panic otherwise.
func (enc *skEncryptor) Encrypt(pt *Plaintext, ct interface{}) {
	enc.encrypt(pt, &enc.uniformSampler, ct)
}

// EncryptSeeded encrypts a slice of plaintext on a pre-allocated slice of cts, using a seed
// to generate a pseudo random public random polynomial.
// The pre-allocated ciphertexts can either have their degree 1 element allocated or nil.
// If the degree 1 element is allocated, then the encryptor will populate this element
// with the pseudo-random public polynomial.
// If the degree 1 element is nil, then the encryptor will only mask the degree 0 element
// using an ephemeral pseudo-random polynomial, which can later be reconstruceted using the seed.
func (enc *skEncryptor) EncryptSeeded(pt *Plaintext, sampler ringqp.UniformSampler, ct interface{}) {
	enc.encrypt(pt, &sampler, ct)
}

func (enc *skEncryptor) encrypt(pt *Plaintext, sampler *ringqp.UniformSampler, ct interface{}) {
	switch el := ct.(type) {
	case *Ciphertext:
		enc.encryptRLWE(pt, sampler, el)
	case *rgsw.Ciphertext:
		enc.encryptRGSW(pt, sampler, el)
	default:
		panic("input ciphertext type unsuported (must be *rlwe.Ciphertext or *rgsw.Ciphertext)")
	}
}

func (enc *skEncryptor) encryptZero(ct *Ciphertext, sampler *ringqp.UniformSampler) {

	ringQ := enc.params.RingQ()
	levelQ := ct.Level()
	c0, c1 := ct.Value[0], ct.Value[1]

	if sampler == nil {
		enc.uniformSampler.ReadLvl(levelQ, -1, ringqp.Poly{Q: c1})
	} else {
		sampler.ReadLvl(levelQ, -1, ringqp.Poly{Q: c1})
	}

	ringQ.MulCoeffsMontgomeryLvl(levelQ, c1, enc.sk.Value.Q, c0) // c0 = NTT(sc1)
	ringQ.NegLvl(levelQ, c0, c0)                                 // c0 = NTT(-sc1)

	if ct.Value[0].IsNTT {
		enc.gaussianSampler.ReadLvl(levelQ, enc.buffQ[0]) // e
		ringQ.NTTLvl(levelQ, enc.buffQ[0], enc.buffQ[0])  // NTT(e)
		ringQ.AddLvl(levelQ, c0, enc.buffQ[0], c0)        // c0 = NTT(-sc1 + e)
	} else {
		ringQ.InvNTTLvl(levelQ, c0, c0)               // c0 = -sc1 + e
		ringQ.InvNTTLvl(levelQ, c1, c1)               // c1 = c1
		enc.gaussianSampler.ReadAndAddLvl(levelQ, c0) // c0 = -sc1 + e
	}
}

func (enc *skEncryptor) encryptRLWE(pt *Plaintext, sampler *ringqp.UniformSampler, ct *Ciphertext) {

	ringQ := enc.params.RingQ()
	levelQ := utils.MinInt(pt.Level(), ct.Level())
	c0, c1 := ct.Value[0], ct.Value[1]
	ctNTT := ct.Value[0].IsNTT
	ptNTT := pt.Value.IsNTT

	if sampler == nil {
		enc.uniformSampler.ReadLvl(utils.MinInt(pt.Level(), ct.Level()), -1, ringqp.Poly{Q: c1})
	} else {
		sampler.ReadLvl(levelQ, -1, ringqp.Poly{Q: c1})
	}

	ringQ.MulCoeffsMontgomeryLvl(levelQ, c1, enc.sk.Value.Q, c0)
	ringQ.NegLvl(levelQ, c0, c0)

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
		enc.gaussianSampler.ReadAndAddLvl(ct.Level(), c0)
		ringQ.InvNTTLvl(levelQ, c1, c1)
	case !ctNTT && !ptNTT:
		ringQ.InvNTTLvl(levelQ, c0, c0)
		ringQ.AddLvl(levelQ, c0, pt.Value, c0)
		enc.gaussianSampler.ReadAndAddLvl(ct.Level(), c0)
		ringQ.InvNTTLvl(levelQ, c1, c1)
	}
}

func (enc *skEncryptor) encryptRGSW(pt *Plaintext, sampler *ringqp.UniformSampler, ct *rgsw.Ciphertext) {

	params := enc.params
	ringQ := params.RingQ()
	levelQ := ct.LevelQ()
	levelP := ct.LevelP()

	decompRNS := params.DecompRNS(levelQ, levelP)
	decompBIT := params.DecompBIT(levelQ, levelP)

	for j := 0; j < decompBIT; j++ {
		for i := 0; i < decompRNS; i++ {
			enc.EncryptZeroSymetricQPNTT(levelQ, levelP, enc.sk.Value, sampler, true, ct.Value[0].Value[i][j])
			enc.EncryptZeroSymetricQPNTT(levelQ, levelP, enc.sk.Value, sampler, true, ct.Value[1].Value[i][j])
		}
	}

	if pt != nil {
		ringQ.MFormLvl(levelQ, pt.Value, enc.buffQP.Q)
		if !pt.Value.IsNTT {
			ringQ.NTTLvl(levelQ, enc.buffQP.Q, enc.buffQP.Q)
		}
		gadget.AddPolyToCiphertext(
			enc.buffQP.Q,
			[]gadget.Ciphertext{ct.Value[0], ct.Value[1]},
			*params.RingQP(),
			params.LogBase2(),
			enc.buffQP.Q)
	}
}

func (enc *pkEncryptor) encryptZero(ciphertext *Ciphertext) {
	ringQP := enc.params.RingQP()
	levelQ := ciphertext.Level()
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
	ringQP.MFormLvl(levelQ, levelP, u, u)

	ct0QP := ringqp.Poly{Q: ciphertext.Value[0], P: buffP0}
	ct1QP := ringqp.Poly{Q: ciphertext.Value[1], P: buffP1}

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
	enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct0QP.Q, ct0QP.P, ciphertext.Value[0])

	// ct1 = (u*pk1 + e1)/P
	enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct1QP.Q, ct1QP.P, ciphertext.Value[1])
}

func (enc *pkEncryptor) encryptNoPRLWE(pt *Plaintext, ciphertext *Ciphertext) {

	ringQ := enc.params.RingQ()
	levelQ := utils.MinInt(pt.Level(), ciphertext.Level())

	ctNTT := ciphertext.Value[0].IsNTT
	ptNTT := pt.Value.IsNTT

	buffQ0 := enc.buffQ[0]

	enc.ternarySampler.ReadLvl(levelQ, buffQ0)
	ringQ.NTTLvl(levelQ, buffQ0, buffQ0)
	ringQ.MFormLvl(levelQ, buffQ0, buffQ0)

	c0, c1 := ciphertext.Value[0], ciphertext.Value[1]

	// ct0 = NTT(u*pk0)
	ringQ.MulCoeffsMontgomeryLvl(levelQ, buffQ0, enc.pk.Value[0].Q, c0)
	// ct1 = NTT(u*pk1)
	ringQ.MulCoeffsMontgomeryLvl(levelQ, buffQ0, enc.pk.Value[1].Q, c1)

	// c0
	switch { // TODO, same as skEncryptor
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
		enc.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), c0)
	case !ctNTT && !ptNTT:
		ringQ.InvNTTLvl(levelQ, c0, c0)
		ringQ.AddLvl(levelQ, c0, pt.Value, c0)
		enc.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), c0)
	}

	// c1
	if ctNTT {
		enc.gaussianSampler.ReadLvl(levelQ, buffQ0)
		ringQ.NTTLvl(levelQ, buffQ0, buffQ0)
		ringQ.AddLvl(levelQ, c1, buffQ0, c1) // ct1 = u*pk1 + e1

	} else {
		ringQ.InvNTTLvl(levelQ, c1, c1)
		// ct[1] = pk[1]*u + e1
		enc.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), c1)
	}

	c1.IsNTT = ciphertext.Value[0].IsNTT
	ciphertext.Resize(ciphertext.Degree(), levelQ)
}

func (enc *pkEncryptor) encryptRLWE(plaintext *Plaintext, ciphertext *Ciphertext) {
	ringQ := enc.params.RingQ()
	levelQ := utils.MinInt(plaintext.Level(), ciphertext.Level())
	c0, c1 := ciphertext.Value[0], ciphertext.Value[1]
	ctNTT := ciphertext.Value[0].IsNTT
	ptNTT := plaintext.Value.IsNTT

	buffQ0 := enc.buffQ[0]

	enc.encryptZero(ciphertext)

	switch {
	case ctNTT && ptNTT:
		ringQ.NTTLvl(levelQ, c0, c0)
		ringQ.NTTLvl(levelQ, c1, c1)
		ringQ.AddLvl(levelQ, c0, plaintext.Value, c0)
	case ctNTT && !ptNTT:
		ringQ.AddLvl(levelQ, c0, plaintext.Value, c0)
		ringQ.NTTLvl(levelQ, c0, c0)
		ringQ.NTTLvl(levelQ, c1, c1)
	case !ctNTT && ptNTT:
		ringQ.InvNTTLvl(levelQ, plaintext.Value, buffQ0)
		ringQ.AddLvl(levelQ, c0, buffQ0, c0)
	case !ctNTT && !ptNTT:
		ringQ.AddLvl(levelQ, c0, plaintext.Value, c0)
	}

	c1.IsNTT = c0.IsNTT
	ciphertext.Resize(ciphertext.Degree(), levelQ)
}

// EncryptZeroSymetricQPNTT generates en encryption of zero under sk.
// levelQ : level of the modulus Q
// levelP : level of the modulus P
// sk     : secret key
// sampler: uniform sampler, if `sampler` is nil, then will sample using the internal sampler.
// montgomery: returns the result in the Montgomery domain.
func (enc *encryptor) EncryptZeroSymetricQPNTT(levelQ, levelP int, sk ringqp.Poly, sampler *ringqp.UniformSampler, montgomery bool, ct interface{}) {

	var c0, c1 ringqp.Poly
	switch ct := ct.(type) {
	case [2]ringqp.Poly:
		c0, c1 = ct[0], ct[1]
	case ringqp.Poly:

		if sampler == nil {
			panic("cannot EncryptZeroSymetricQPNTT: sampler cannot be nil if ct.(type) is ringqp.Poly")
		}

		c0, c1 = ct, enc.buffQP
	default:
		panic("cannot EncryptZeroSymetricQPNTT: ct must be either [2]3ingqp.Poly or ringqp.Poly")
	}

	ringQP := enc.params.RingQP()

	// ct = (e, 0)
	enc.gaussianSampler.ReadLvl(levelQ, c0.Q)
	if levelP != -1 {
		ringQP.ExtendBasisSmallNormAndCenter(c0.Q, levelP, nil, c0.P)
	}

	ringQP.NTTLvl(levelQ, levelP, c0, c0)

	// If montgomery, then ct[1] is assumed to be sampled in of the Montgomery domain,
	// thus -as will also be in the montgomery domain (s is by default), therefore 'e'
	// must be switched to the montgomery domain.
	if montgomery {
		ringQP.MFormLvl(levelQ, levelP, c0, c0)
	}

	// ct = (e, a)
	if sampler == nil {
		enc.uniformSampler.ReadLvl(levelQ, levelP, c1)
	} else {
		sampler.ReadLvl(levelQ, levelP, c1)
	}

	// (-a*sk + e, a)
	ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, c1, sk, c0)
}

// ShallowCopy creates a shallow copy of this pkEncryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *pkEncryptor) ShallowCopy() Encryptor {
	return &pkEncryptor{*enc.encryptor.ShallowCopy(), enc.pk}
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

// WithKey creates a shallow copy of this encryptor with a new key in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *encryptor) WithKey(key interface{}) Encryptor {
	return enc.ShallowCopy().setKey(key)
}

func (enc *encryptor) setKey(key interface{}) Encryptor {
	switch key := key.(type) {
	case *PublicKey:
		if key.Value[0].Q.N() != enc.params.N() || key.Value[1].Q.N() != enc.params.N() {
			panic("cannot setKey: pk ring degree does not match params ring degree")
		}
		return &pkEncryptor{*enc, key}
	case *SecretKey:
		if key.Value.Q.N() != enc.params.N() {
			panic("cannot setKey: sk ring degree does not match params ring degree")
		}
		return &skEncryptor{*enc, key}
	default:
		panic("cannot setKey: key must be either *rlwe.PublicKey or *rlwe.SecretKey")
	}
}
