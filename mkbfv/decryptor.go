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
	MergeDec(c0 *ring.Poly, partialKeys []*ring.Poly) *bfv.Plaintext
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

	return &mkDecryptor{
		params:          params,
		ringQ:           ringQ,
		ringQMul:        ringQMul,
		pHalf:           pHalf,
		samplerGaussian: sampler,
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
func (dec *mkDecryptor) MergeDec(c0 *ring.Poly, partialKeys []*ring.Poly) *bfv.Plaintext {

	plaintext := new(bfv.Element)

	res := dec.ringQ.NewPoly()

	dec.ringQ.Copy(c0, res)

	for _, k := range partialKeys {
		dec.ringQ.Add(res, k, res)
	}

	plaintext.SetValue([]*ring.Poly{res})

	return plaintext.Plaintext()
}
