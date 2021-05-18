package mkrlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKDecryptor is a type for mkrlwe decryptor in a multi key context
type MKDecryptor interface {
	PartDec(ct *ring.Poly, level uint64, sk *MKSecretKey) *ring.Poly
	MergeDec(c0 *ring.Poly, level uint64, partialKeys []*ring.Poly) *ring.Poly
}

type mkDecryptor struct {
	params          *rlwe.Parameters
	ringQ           *ring.Ring
	samplerGaussian *ring.GaussianSampler
}

// NewMKDecryptor returns a decryptor for rlwe in a multi key context
// the standard deviation for the partial decryption must be provided
func NewMKDecryptor(params *rlwe.Parameters, sigmaSmudging float64) MKDecryptor {

	ringQ := GetRingQ(params)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	sampler := ring.NewGaussianSampler(prng, ringQ, sigmaSmudging, uint64(6*sigmaSmudging))

	return &mkDecryptor{
		params:          params,
		ringQ:           ringQ,
		samplerGaussian: sampler,
	}

}

// PartDec computes a partial decription key for the ciphertext component of a given participant
// for participant i, ski and cti must be used
func (dec *mkDecryptor) PartDec(ct *ring.Poly, level uint64, sk *MKSecretKey) *ring.Poly {

	// mu_i = c_i * sk_i + e_i mod q
	out := dec.samplerGaussian.ReadLvlNew(level)
	dec.ringQ.NTTLvl(level, out, out)

	dec.ringQ.MulCoeffsMontgomeryAndAddLvl(level, ct, sk.Key.Value, out)

	out.Coeffs = out.Coeffs[:level+1]

	return out
}

// MergeDec merges the partial decription parts and returns the plaintext. The first component of the ciphertext vector must be provided (c0)
func (dec *mkDecryptor) MergeDec(c0 *ring.Poly, level uint64, partialKeys []*ring.Poly) *ring.Poly {

	res := dec.ringQ.NewPoly()
	dec.ringQ.CopyLvl(level, c0, res)

	for _, k := range partialKeys {
		dec.ringQ.AddLvl(level, res, k, res)
	}

	dec.ringQ.ReduceLvl(level, res, res)
	res.Coeffs = res.Coeffs[:level+1]

	return res
}
