package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"math/big"
)

// Encryptor a generic RLWE encryption interface.
type Encryptor interface {
	// Encrypt encrypts the input plaintext and write the result on ct.
	// The encryption algorithm depends on the implementor.
	Encrypt(pt *Plaintext, ct *Ciphertext)

	// EncryptFromCRP encrypts the input plaintext and writes the result on ct.
	// The encryption algorithm depends on the implementor.
	EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ctOut *Ciphertext)

	// EncryptFromCRP encrypts the input plaintext and write the result on ct.
	EncryptRGSW(pt *Plaintext, ct *RGSWCiphertext)
}

// encryptorBase is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptorBase struct {
	params Parameters

	ringQ *ring.Ring
	ringP *ring.Ring

	poolQ [2]*ring.Poly
	poolP [3]*ring.Poly

	baseconverter   *ring.FastBasisExtender
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSamplerQ *ring.UniformSampler
	uniformSamplerP *ring.UniformSampler
}

type pkEncryptor struct {
	encryptorBase
	pk *PublicKey
}

type skEncryptor struct {
	encryptorBase
	sk *SecretKey
}

func (enc *pkEncryptor) Encrypt(pt *Plaintext, ct *Ciphertext) {
	enc.encrypt(pt, enc.pk, ct)
}

func (enc *pkEncryptor) EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ct *Ciphertext) {
	enc.encrypt(pt, enc.pk, ct)
}

func (enc *skEncryptor) Encrypt(pt *Plaintext, ct *Ciphertext) {
	enc.encrypt(pt, enc.sk, ct)
}

func (enc *skEncryptor) EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ct *Ciphertext) {
	enc.encrypt(pt, enc.sk, ct)
}

// NewEncryptor creates a new Encryptor
// Accepts either a secret-key or a public-key
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	switch key := key.(type) {
	case *PublicKey:
		if key.Value[0].Q.Degree() != params.N() || key.Value[1].Q.Degree() != params.N() {
			panic("cannot newEncryptor: pk ring degree does not match params ring degree")
		}
		return &pkEncryptor{newEncryptorBase(params), key}
	case *SecretKey:
		if key.Value.Q.Degree() != params.N() {
			panic("cannot newEncryptor: sk ring degree does not match params ring degree")
		}
		return &skEncryptor{newEncryptorBase(params), key}
	default:
		panic("key must be either rlwe.PublicKey or rlwe.SecretKey")
	}
}

func newEncryptorBase(params Parameters) encryptorBase {

	ringQ := params.RingQ()
	ringP := params.RingP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var poolP [3]*ring.Poly
	var bc *ring.FastBasisExtender
	var uniformSamplerP *ring.UniformSampler
	if params.PCount() != 0 {
		poolP = [3]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
		bc = ring.NewFastBasisExtender(ringQ, ringP)
		uniformSamplerP = ring.NewUniformSampler(prng, ringP)
	}

	return encryptorBase{
		params:          params,
		ringQ:           ringQ,
		ringP:           ringP,
		poolQ:           [2]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
		poolP:           poolP,
		baseconverter:   bc,
		gaussianSampler: ring.NewGaussianSampler(prng, ringQ, params.Sigma(), int(6*params.Sigma())),
		ternarySampler:  ring.NewTernarySampler(prng, ringQ, 0.5, false),
		uniformSamplerQ: ring.NewUniformSampler(prng, ringQ),
		uniformSamplerP: uniformSamplerP,
	}
}

// Encrypt encrypts the input Plaintext and write the result on the output Ciphertext.
func (enc *encryptorBase) encrypt(plaintext *Plaintext, key interface{}, ciphertext *Ciphertext) {
	switch key := key.(type) {
	case *PublicKey:

		if key.Value[0].Q.Degree() != enc.params.N() || key.Value[1].Q.Degree() != enc.params.N() {
			panic("cannot newEncryptor: pk ring degree does not match params ring degree")
		}

		enc.uniformSamplerQ.ReadLvl(utils.MinInt(plaintext.Level(), ciphertext.Level()), ciphertext.Value[1])

		if enc.baseconverter != nil {
			enc.encryptPk(plaintext, key, ciphertext)
		} else {
			enc.encryptPkNoP(plaintext, key, ciphertext)
		}

	case *SecretKey:

		if key.Value.Q.Degree() != enc.params.N() {
			panic("cannot newEncryptor: sk ring degree does not match params ring degree")
		}

		enc.uniformSamplerQ.ReadLvl(utils.MinInt(plaintext.Level(), ciphertext.Level()), ciphertext.Value[1])

		enc.encryptSk(plaintext, key, ciphertext)

	default:
		panic("key must be either rlwe.PublicKey or rlwe.SecretKey")
	}
}

func (enc *encryptorBase) encryptFromCRP(plaintext *Plaintext, key interface{}, crp *ring.Poly, ciphertext *Ciphertext) {
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

func (enc *encryptorBase) encryptPk(plaintext *Plaintext, pk *PublicKey, ciphertext *Ciphertext) {
	ringQ := enc.ringQ
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
	enc.baseconverter.ModDownQPtoQ(levelQ, levelP, ct0QP.Q, ct0QP.P, ct0QP.Q)

	// ct1 = (u*pk1 + e1)/P
	enc.baseconverter.ModDownQPtoQ(levelQ, levelP, ct1QP.Q, ct1QP.P, ct1QP.Q)

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

// EncryptRGSW encrypts the input Plaintext and writes the result on the output RGSW ciphertext.
func (enc *pkEncryptor) EncryptRGSW(plaintext *Plaintext, ciphertext *RGSWCiphertext) {
	panic("method not implemented")
}

func (enc *encryptorBase) encryptPkNoP(plaintext *Plaintext, pk *PublicKey, ciphertext *Ciphertext) {
	levelQ := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := enc.poolQ[0]

	ringQ := enc.ringQ

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

func (enc *encryptorBase) encryptSk(plaintext *Plaintext, sk *SecretKey, ciphertext *Ciphertext) {

	ringQ := enc.ringQ

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

func (enc *skEncryptor) EncryptRGSW(plaintext *Plaintext, ciphertext *RGSWCiphertext) {

	params := enc.params
	ringQ := params.RingQ()
	ringP := params.RingP()
	ringQP := params.RingQP()
	isNTT := ciphertext.Value[0][0][0].Q.IsNTT
	levelQ := ciphertext.LevelQ()
	levelP := ciphertext.LevelP()

	var pBigInt *big.Int
	if levelP == params.PCount()-1 {
		pBigInt = ringP.ModulusBigint
	} else {
		P := ringP.Modulus
		pBigInt = new(big.Int).SetUint64(P[0])
		for i := 1; i < levelP+1; i++ {
			pBigInt.Mul(pBigInt, ring.NewUint(P[i]))
		}
	}

	ptTimesP := enc.poolQ[1]

	if plaintext != nil {
		ringQ.MulScalarBigintLvl(levelQ, plaintext.Value, pBigInt, ptTimesP)
		if !plaintext.Value.IsNTT {
			ringQ.NTTLvl(levelQ, ptTimesP, ptTimesP)
		}
	}

	alpha := levelP + 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	var index int
	for i := 0; i < beta; i++ {

		enc.encryptZeroSymetricQP(levelQ, levelP, enc.sk.Value, true, isNTT, ciphertext.Value[i][0][0], ciphertext.Value[i][0][1])
		enc.encryptZeroSymetricQP(levelQ, levelP, enc.sk.Value, true, isNTT, ciphertext.Value[i][1][0], ciphertext.Value[i][1][1])

		if plaintext != nil {
			for j := 0; j < alpha; j++ {

				index = i*alpha + j

				// It handles the case where nb pj does not divide nb qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.Modulus[index]
				p0tmp := ptTimesP.Coeffs[index]
				p1tmp := ciphertext.Value[i][0][0].Q.Coeffs[index]
				p2tmp := ciphertext.Value[i][1][1].Q.Coeffs[index]

				for w := 0; w < ringQ.N; w++ {
					p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
					p2tmp[w] = ring.CRed(p2tmp[w]+p0tmp[w], qi)
				}
			}
		}

		ringQP.MFormLvl(levelQ, levelP, ciphertext.Value[i][0][0], ciphertext.Value[i][0][0])
		ringQP.MFormLvl(levelQ, levelP, ciphertext.Value[i][0][1], ciphertext.Value[i][0][1])
		ringQP.MFormLvl(levelQ, levelP, ciphertext.Value[i][1][0], ciphertext.Value[i][1][0])
		ringQP.MFormLvl(levelQ, levelP, ciphertext.Value[i][1][1], ciphertext.Value[i][1][1])
	}
}

func (enc *encryptorBase) encryptZeroSymetricQP(levelQ, levelP int, sk PolyQP, sample, ntt bool, a, b PolyQP) {

	params := enc.params
	ringQP := params.RingQP()

	if ntt {
		enc.gaussianSampler.ReadLvl(levelQ, a.Q)
		ringQP.ExtendBasisSmallNormAndCenter(a.Q, levelP, nil, a.P)
		ringQP.NTTLvl(levelQ, levelP, a, a)
	}

	if sample {
		enc.uniformSamplerQ.ReadLvl(levelQ, b.Q)
		enc.uniformSamplerP.ReadLvl(levelP, b.P)
	}

	ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, b, sk, a)

	if !ntt {
		ringQP.InvNTTLvl(levelQ, levelP, a, a)
		ringQP.InvNTTLvl(levelQ, levelP, b, b)

		e := PolyQP{Q: enc.poolQ[0], P: enc.poolP[0]}
		enc.gaussianSampler.ReadLvl(levelQ, e.Q)
		ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
		ringQP.AddLvl(levelQ, levelP, a, e, a)
	}
}

func (enc *encryptorBase) encryptZeroSymetricQ(levelQ int, sk *ring.Poly, sample, ntt bool, a, b *ring.Poly) {

	params := enc.params
	ringQ := params.RingQ()

	if ntt {
		enc.gaussianSampler.ReadLvl(levelQ, a)
		ringQ.NTTLazyLvl(levelQ, a, a)
	}

	if sample {
		enc.uniformSamplerQ.ReadLvl(levelQ, b)
	}

	ringQ.MulCoeffsMontgomeryAndSubLvl(levelQ, b, sk, a)

	if !ntt {
		ringQ.InvNTTLvl(levelQ, a, a)
		ringQ.InvNTTLvl(levelQ, b, b)
		enc.gaussianSampler.ReadAndAddLvl(levelQ, a)
	}
}
