package mkbfv

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKDecryptor is a type for bfv decryptor in a multi key context
type MKDecryptor interface {
	PartDec(ct *ring.Poly, sk *MKSecretKey) *ring.Poly
	MergeDec(c0 *ring.Poly, partialKeys []*ring.Poly) *ckks.Plaintext
}

type mkDecryptor struct {
	params          *ckks.Parameters
	ringQ           *ring.Ring
	samplerGaussian *ring.GaussianSampler
}

// NewMKDecryptor returns a decryptor for bfv in a multi key context
func NewMKDecryptor(params *ckks.Parameters) MKDecryptor {

	ringQ := GetRingQ(params)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	sampler := GetGaussianSampler(params, ringQ, prng)

	return &mkDecryptor{
		params:          params,
		ringQ:           ringQ,
		samplerGaussian: sampler,
	}

}

// PartDec computes a partial decription key for the ciphertext component of a given participant
// for participant i, ski and cti must be used
func (dec *mkDecryptor) PartDec(ct *ring.Poly, sk *MKSecretKey) *ring.Poly {

	// mu_i = c_i * sk_i + e_i mod q

	out := dec.samplerGaussian.ReadNew() // TODO: in paper they want sigma > 3.2 for this error... but they don't tell how much...
	dec.ringQ.NTT(out, out)

	tmp := dec.ringQ.NewPoly()
	dec.ringQ.NTTLazy(ct, tmp)

	dec.ringQ.MulCoeffsMontgomeryAndAdd(tmp, sk.key.Value, out)

	return out
}

// MergeDec merges the partial decription parts and returns the plaintext. The first component of the ciphertext vector must be provided (c0)
func (dec *mkDecryptor) MergeDec(c0 *ring.Poly, partialKeys []*ring.Poly) *ckks.Plaintext {

	plaintext := new(ckks.Element)

	res := dec.ringQ.NewPoly()
	dec.ringQ.NTTLazy(c0, res)

	for _, k := range partialKeys {
		dec.ringQ.Add(res, k, res)
	}

	dec.ringQ.Reduce(res, res)
	dec.ringQ.InvNTT(res, res)

	plaintext.SetValue([]*ring.Poly{res})

	return plaintext.Plaintext()
}
