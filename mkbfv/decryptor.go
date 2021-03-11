package mkbfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKDecryptor is a type for bfv decryptor in a multi key context
type MKDecryptor interface {
	PartDec(ct *ring.Poly, sk *MKSecretKey, out *ring.Poly)
	MergeDec(c0 *ring.Poly, partialKeys []*ring.Poly, out *bfv.Plaintext)
}

type mkDecryptor struct {
	params          *bfv.Parameters
	ringQ           *ring.Ring
	ringQMul        *ring.Ring
	pHalf           *big.Int
	samplerGaussian *ring.GaussianSampler
	polyPoolQ1      []*ring.Poly
	polyPoolQ2      []*ring.Poly
	convertor       *ring.FastBasisExtender
}

// NewMKDecryptor returns a decryptor for bfv in a multi key context
func NewMKDecryptor(params *bfv.Parameters) MKDecryptor {

	ringQ := GetRingQ(params)
	ringQMul := GetRingQMul(params)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	sampler := GetGaussianSampler(params, ringQ, prng)
	convertor := ring.NewFastBasisExtender(ringQ, ringQMul)

	pHalf := new(big.Int).Rsh(ringQMul.ModulusBigint, 1)

	pool1 := make([]*ring.Poly, 4)
	pool2 := make([]*ring.Poly, 4)

	for i := 0; i < len(pool1); i++ {
		pool1[i] = ringQ.NewPoly()
		pool2[i] = ringQ.NewPoly()
	}

	return &mkDecryptor{
		params:          params,
		ringQ:           ringQ,
		ringQMul:        ringQMul,
		pHalf:           pHalf,
		samplerGaussian: sampler,
		polyPoolQ1:      pool1,
		polyPoolQ2:      pool2,
		convertor:       convertor,
	}

}

// PartDec computes a partial decription key for the ciphertext component of a given participant
// for participant i, ski and cti must be used
func (dec *mkDecryptor) PartDec(ct *ring.Poly, sk *MKSecretKey, out *ring.Poly) {

	// mu_i = c_i * sk_i + e_i mod q

	out = dec.samplerGaussian.ReadNew()

	dec.ringQ.MulCoeffsAndAdd(ct, sk.key.Value, out)

}

// MergeDec merges the partial decription parts and returns the plaintext. The first component of the ciphertext vector must be provided (c0)
func (dec *mkDecryptor) MergeDec(c0 *ring.Poly, partialKeys []*ring.Poly, out *bfv.Plaintext) {

	res := out.Value()[0]

	dec.ringQ.Copy(c0, res)

	for _, k := range partialKeys {
		dec.ringQ.MulCoeffsAndAdd(res, k, res)
	}

	dec.quantize(out.El())
}

// function taken from bfv evaluator to scale by t/q
func (dec *mkDecryptor) quantize(ctOut *bfv.Element) {

	levelQ := uint64(len(dec.ringQ.Modulus) - 1)
	levelQMul := uint64(len(dec.ringQMul.Modulus) - 1)

	c2Q1 := dec.polyPoolQ1
	c2Q2 := dec.polyPoolQ2

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.Value() {
		dec.ringQ.InvNTTLazy(c2Q1[i], c2Q1[i])
		dec.ringQMul.InvNTTLazy(c2Q2[i], c2Q2[i])

		// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
		dec.convertor.ModDownSplitQP(levelQ, levelQMul, c2Q1[i], c2Q2[i], c2Q2[i])

		// Centers (ct(x)Q -> P)/Q by (P-1)/2 and extends ((ct(x)Q -> P)/Q) to the basis Q
		dec.ringQMul.AddScalarBigint(c2Q2[i], dec.pHalf, c2Q2[i])
		dec.convertor.ModUpSplitPQ(levelQMul, c2Q2[i], ctOut.Value()[i])
		dec.ringQ.SubScalarBigint(ctOut.Value()[i], dec.pHalf, ctOut.Value()[i])

		// Option (2) (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
		dec.ringQ.MulScalar(ctOut.Value()[i], dec.params.T(), ctOut.Value()[i])
	}
}
