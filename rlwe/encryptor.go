package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"math/big"
)

// Encryptor a generic RLWE encryption interface.
type Encryptor interface {
	// Encrypt encrypts the input plaintext and write the result on ctOut.
	// The encryption algorithm depends on the implementor.
	Encrypt(pt *Plaintext, ctOut *Ciphertext)

	// EncryptFromCRP encrypts the input plaintext and writes the result in ctOut.
	// The encryption algorithm depends on the implementor.
	EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ctOut *Ciphertext)

	EncryptRGSW(plaintext *Plaintext, ciphertext *RGSWCiphertext)
}

// encryptorBase is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptorBase struct {
	params Parameters

	ringQ *ring.Ring
	ringP *ring.Ring

	poolQ [2]*ring.Poly
	poolP [3]*ring.Poly

	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSamplerQ *ring.UniformSampler
	uniformSamplerP *ring.UniformSampler
}

type pkEncryptor struct {
	encryptorBase
	pk            *PublicKey
	baseconverter *ring.FastBasisExtender
}

type pkFastEncryptor struct {
	encryptorBase
	pk *PublicKey
}

type skEncryptor struct {
	encryptorBase
	sk *SecretKey
}

// NewEncryptor instatiates a new generic RLWE Encryptor. The key argument can
// be either a *rlwe.PublicKey or a *rlwe.SecretKey.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	switch key := key.(type) {
	case *PublicKey:
		if key.Value[0].Q.Degree() != params.N() || key.Value[1].Q.Degree() != params.N() {
			panic("cannot newEncryptor: pk ring degree does not match params ring degree")
		}
		encryptorBase := newEncryptorBase(params)
		if params.PCount() > 0 {
			baseconverter := ring.NewFastBasisExtender(params.ringQ, params.ringP)
			return &pkEncryptor{encryptorBase, key, baseconverter}
		}
		return &pkFastEncryptor{encryptorBase, key}
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
	if params.PCount() != 0 {
		poolP = [3]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return encryptorBase{
		params:          params,
		ringQ:           ringQ,
		ringP:           ringP,
		poolQ:           [2]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
		poolP:           poolP,
		gaussianSampler: ring.NewGaussianSampler(prng, ringQ, params.Sigma(), int(6*params.Sigma())),
		ternarySampler:  ring.NewTernarySampler(prng, ringQ, 0.5, false),
		uniformSamplerQ: ring.NewUniformSampler(prng, ringQ),
		uniformSamplerP: ring.NewUniformSampler(prng, ringP),
	}
}

// NewFastEncryptor instantiates a new  generic RLWE Encryptor.
// This encryptor's Encrypt method first encrypts zero in Q and then adds the plaintext.
// This method is faster than the normal encryptor but result in a noisier ciphertext.
func NewFastEncryptor(params Parameters, key *PublicKey) Encryptor {
	return &pkFastEncryptor{newEncryptorBase(params), key}
}

// Encrypt encrypts the input Plaintext and write the result in ctOut.
func (encryptor *pkEncryptor) Encrypt(plaintext *Plaintext, ctOut *Ciphertext) {
	ringQ := encryptor.ringQ
	ringQP := encryptor.params.RingQP()

	levelQ := utils.MinInt(plaintext.Level(), ctOut.Level())
	levelP := 0

	poolQ0 := encryptor.poolQ[0]
	poolP0 := encryptor.poolP[0]
	poolP1 := encryptor.poolP[1]
	poolP2 := encryptor.poolP[2]

	// We sample a R-WLE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)

	ciphertextNTT := ctOut.Value[0].IsNTT

	u := PolyQP{Q: poolQ0, P: poolP2}

	encryptor.ternarySampler.ReadLvl(levelQ, u.Q)
	ringQP.ExtendBasisSmallNormAndCenter(u.Q, levelP, nil, u.P)

	// (#Q + #P) NTT
	ringQP.NTTLvl(levelQ, levelP, u, u)
	ringQP.MFormLvl(levelQ, levelP, u, u)

	ct0QP := PolyQP{Q: ctOut.Value[0], P: poolP0}
	ct1QP := PolyQP{Q: ctOut.Value[1], P: poolP1}

	// ct0 = u*pk0
	// ct1 = u*pk1
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, u, encryptor.pk.Value[0], ct0QP)
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, u, encryptor.pk.Value[1], ct1QP)

	// 2*(#Q + #P) NTT
	ringQP.InvNTTLvl(levelQ, levelP, ct0QP, ct0QP)
	ringQP.InvNTTLvl(levelQ, levelP, ct1QP, ct1QP)

	e := PolyQP{Q: poolQ0, P: poolP2}

	encryptor.gaussianSampler.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.AddLvl(levelQ, levelP, ct0QP, e, ct0QP)

	encryptor.gaussianSampler.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.AddLvl(levelQ, levelP, ct1QP, e, ct1QP)

	// ct0 = (u*pk0 + e0)/P
	encryptor.baseconverter.ModDownQPtoQ(levelQ, levelP, ct0QP.Q, ct0QP.P, ct0QP.Q)

	// ct1 = (u*pk1 + e1)/P
	encryptor.baseconverter.ModDownQPtoQ(levelQ, levelP, ct1QP.Q, ct1QP.P, ct1QP.Q)

	if ciphertextNTT {

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		}

		// 2*#Q NTT
		ringQ.NTTLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.NTTLvl(levelQ, ctOut.Value[1], ctOut.Value[1])

		if plaintext.Value.IsNTT {
			// ct0 = (u*pk0 + e0)/P + m
			ringQ.AddLvl(levelQ, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		}

	} else {

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		} else {
			ringQ.InvNTTLvl(levelQ, plaintext.Value, poolQ0)
			ringQ.AddLvl(levelQ, ctOut.Value[0], poolQ0, ctOut.Value[0])
		}
	}

	ctOut.Value[1].IsNTT = ctOut.Value[0].IsNTT
	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:levelQ+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:levelQ+1]
}

func (encryptor *pkEncryptor) EncryptRGSW(plaintext *Plaintext, ciphertext *RGSWCiphertext) {
	panic("method not implemented")
}

func (encryptor *pkFastEncryptor) EncryptRGSW(plaintext *Plaintext, ciphertext *RGSWCiphertext) {
	panic("method not implemented")
}

// Encrypt encrypts the input Plaintext and write the result in ctOut.
// It first encrypts zero in Q and then adds the plaintext.
// This method is faster than the normal encryptor but result in a noisier ciphertext.
func (encryptor *pkFastEncryptor) Encrypt(plaintext *Plaintext, ctOut *Ciphertext) {
	levelQ := utils.MinInt(plaintext.Level(), ctOut.Level())

	poolQ0 := encryptor.poolQ[0]

	ringQ := encryptor.ringQ

	ciphertextNTT := ctOut.Value[0].IsNTT

	encryptor.ternarySampler.ReadLvl(levelQ, poolQ0)
	ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
	ringQ.MFormLvl(levelQ, poolQ0, poolQ0)

	// ct0 = u*pk0
	ringQ.MulCoeffsMontgomeryLvl(levelQ, poolQ0, encryptor.pk.Value[0].Q, ctOut.Value[0])
	// ct1 = u*pk1
	ringQ.MulCoeffsMontgomeryLvl(levelQ, poolQ0, encryptor.pk.Value[1].Q, ctOut.Value[1])

	if ciphertextNTT {

		// ct1 = u*pk1 + e1
		encryptor.gaussianSampler.ReadLvl(levelQ, poolQ0)
		ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
		ringQ.AddLvl(levelQ, ctOut.Value[1], poolQ0, ctOut.Value[1])

		// ct0 = u*pk0 + e0
		encryptor.gaussianSampler.ReadLvl(levelQ, poolQ0)

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, poolQ0, plaintext.Value, poolQ0)
			ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
			ringQ.AddLvl(levelQ, ctOut.Value[0], poolQ0, ctOut.Value[0])
		} else {
			ringQ.NTTLvl(levelQ, poolQ0, poolQ0)
			ringQ.AddLvl(levelQ, ctOut.Value[0], poolQ0, ctOut.Value[0])
			ringQ.AddLvl(levelQ, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		}

	} else {

		ringQ.InvNTTLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.InvNTTLvl(levelQ, ctOut.Value[1], ctOut.Value[1])

		// ct[0] = pk[0]*u + e0
		encryptor.gaussianSampler.ReadAndAddLvl(ctOut.Level(), ctOut.Value[0])

		// ct[1] = pk[1]*u + e1
		encryptor.gaussianSampler.ReadAndAddLvl(ctOut.Level(), ctOut.Value[1])

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		} else {
			ringQ.InvNTTLvl(levelQ, plaintext.Value, poolQ0)
			ringQ.AddLvl(levelQ, ctOut.Value[0], poolQ0, ctOut.Value[0])
		}
	}

	ctOut.Value[1].IsNTT = ctOut.Value[0].IsNTT

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:levelQ+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:levelQ+1]
}

// Encrypt encrypts the input Plaintext and write the result in ctOut.
func (encryptor *skEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.uniformSamplerQ.ReadLvl(utils.MinInt(plaintext.Level(), ciphertext.Level()), ciphertext.Value[1])
	encryptor.encrypt(plaintext, ciphertext)
}

// EncryptFromCRP encrypts the input Plaintext given the uniformly random element c1 and write the result in ctOut.
func (encryptor *skEncryptor) EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ctOut *Ciphertext) {
	ring.CopyValues(crp, ctOut.Value[1])
	encryptor.encrypt(plaintext, ctOut)
}

func (encryptor *skEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	ringQ := encryptor.ringQ

	levelQ := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := encryptor.poolQ[0]

	ciphertextNTT := ciphertext.Value[0].IsNTT

	ringQ.MulCoeffsMontgomeryLvl(levelQ, ciphertext.Value[1], encryptor.sk.Value.Q, ciphertext.Value[0])
	ringQ.NegLvl(levelQ, ciphertext.Value[0], ciphertext.Value[0])

	if ciphertextNTT {

		encryptor.gaussianSampler.ReadLvl(levelQ, poolQ0)

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

		encryptor.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), ciphertext.Value[0])

		ringQ.InvNTTLvl(levelQ, ciphertext.Value[1], ciphertext.Value[1])

		ciphertext.Value[0].IsNTT = false
		ciphertext.Value[1].IsNTT = false

	}

	ciphertext.Value[0].Coeffs = ciphertext.Value[0].Coeffs[:levelQ+1]
	ciphertext.Value[1].Coeffs = ciphertext.Value[1].Coeffs[:levelQ+1]
}

func (encryptor *encryptorBase) EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ctOut *Ciphertext) {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *skEncryptor) EncryptRGSW(plaintext *Plaintext, ciphertext *RGSWCiphertext) {

	params := encryptor.params
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

	ptTimesP := encryptor.poolQ[1]

	ringQ.MulScalarBigintLvl(levelQ, plaintext.Value, pBigInt, ptTimesP)
	if !plaintext.Value.IsNTT {
		ringQ.NTTLvl(levelQ, ptTimesP, ptTimesP)
	}

	alpha := levelP + 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	var index int
	for i := 0; i < beta; i++ {
		encryptor.encryptZeroSymetricQP(levelQ, levelP, encryptor.sk.Value, true, isNTT, ciphertext.Value[i][0][0], ciphertext.Value[i][0][1])
		encryptor.encryptZeroSymetricQP(levelQ, levelP, encryptor.sk.Value, true, isNTT, ciphertext.Value[i][1][0], ciphertext.Value[i][1][1])

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

		ringQP.MFormLvl(levelQ, levelP, ciphertext.Value[i][0][0], ciphertext.Value[i][0][0])
		ringQP.MFormLvl(levelQ, levelP, ciphertext.Value[i][0][1], ciphertext.Value[i][0][1])
		ringQP.MFormLvl(levelQ, levelP, ciphertext.Value[i][1][0], ciphertext.Value[i][1][0])
		ringQP.MFormLvl(levelQ, levelP, ciphertext.Value[i][1][1], ciphertext.Value[i][1][1])
	}
}

func (encryptor *encryptorBase) encryptZeroSymetricQP(levelQ, levelP int, sk PolyQP, sample, ntt bool, a, b PolyQP) {

	params := encryptor.params
	ringQP := params.RingQP()

	if ntt {
		encryptor.gaussianSampler.ReadLvl(levelQ, a.Q)
		ringQP.ExtendBasisSmallNormAndCenter(a.Q, levelP, nil, a.P)
		ringQP.NTTLvl(levelQ, levelP, a, a)
	}

	if sample {
		encryptor.uniformSamplerQ.ReadLvl(levelQ, b.Q)
		encryptor.uniformSamplerP.ReadLvl(levelP, b.P)
	}

	ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, b, sk, a)

	if !ntt {
		ringQP.InvNTTLvl(levelQ, levelP, a, a)
		ringQP.InvNTTLvl(levelQ, levelP, b, b)

		e := PolyQP{Q: encryptor.poolQ[0], P: encryptor.poolP[0]}
		encryptor.gaussianSampler.ReadLvl(levelQ, e.Q)
		ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
		ringQP.AddLvl(levelQ, levelP, a, e, a)
	}
}

func (encryptor *encryptorBase) encryptZeroSymetricQ(levelQ int, sk *ring.Poly, sample, ntt bool, a, b *ring.Poly) {

	params := encryptor.params
	ringQ := params.RingQ()

	if ntt {
		encryptor.gaussianSampler.ReadLvl(levelQ, a)
		ringQ.NTTLazyLvl(levelQ, a, a)
	}

	if sample {
		encryptor.uniformSamplerQ.ReadLvl(levelQ, b)
	}

	ringQ.MulCoeffsMontgomeryAndSubLvl(levelQ, b, sk, a)

	if !ntt {
		ringQ.InvNTTLvl(levelQ, a, a)
		ringQ.InvNTTLvl(levelQ, b, b)
		encryptor.gaussianSampler.ReadAndAddLvl(levelQ, a)
	}
}
