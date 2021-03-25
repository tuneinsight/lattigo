package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKDecryptor is a type for ckks decryptor in a multi key context
type MKDecryptor interface {
	PartDec(ct *ring.Poly, level uint64, sk *MKSecretKey) *ring.Poly
	MergeDec(c0 *ring.Poly, scale float64, level uint64, partialKeys []*ring.Poly) *ckks.Plaintext
}

type mkDecryptor struct {
	params          *ckks.Parameters
	ringQ           *ring.Ring
	samplerGaussian *ring.GaussianSampler
}

// NewMKDecryptor returns a decryptor for ckks in a multi key context
// the standard deviation for the partial decryption must be provided
func NewMKDecryptor(params *ckks.Parameters, sigmaSmudging float64) MKDecryptor {

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

	out := dec.samplerGaussian.ReadLvlNew(level) // TODO: in paper they want sigma > 3.2 for this error... but they don't tell how much...
	dec.ringQ.NTTLvl(level, out, out)

	tmp := dec.ringQ.NewPoly()
	dec.ringQ.CopyLvl(level, ct, tmp)

	dec.ringQ.MulCoeffsMontgomeryAndAddLvl(level, tmp, sk.key.Value, out)

	out.Coeffs = out.Coeffs[:level+1]

	return out
}

// MergeDec merges the partial decription parts and returns the plaintext. The first component of the ciphertext vector must be provided (c0)
func (dec *mkDecryptor) MergeDec(c0 *ring.Poly, scale float64, level uint64, partialKeys []*ring.Poly) *ckks.Plaintext {

	plaintext := ckks.NewPlaintext(dec.params, level, scale)

	res := dec.ringQ.NewPoly()
	dec.ringQ.CopyLvl(level, c0, res)

	for _, k := range partialKeys {
		dec.ringQ.AddLvl(level, res, k, res)
	}

	dec.ringQ.ReduceLvl(level, res, res)
	res.Coeffs = res.Coeffs[:level+1]

	plaintext.SetValue([]*ring.Poly{res})

	return plaintext.Plaintext()
}
