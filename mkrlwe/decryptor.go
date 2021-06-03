package mkrlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// MKDecryptor is a type for mkrlwe decryptor in a multi key context
type MKDecryptor interface {
	PartDec(ct *rlwe.Element, level uint64, sk *MKSecretKey, sigmaSmudging float64) *ring.Poly
	MergeDec(ct *rlwe.Element, level uint64, partialKeys []*ring.Poly) *ring.Poly
}

type mkDecryptor struct {
	params *rlwe.Parameters
	ringQ  *ring.Ring
}

// NewMKDecryptor returns a decryptor for rlwe in a multi key context
// the standard deviation for the partial decryption must be provided
func NewMKDecryptor(params *rlwe.Parameters) MKDecryptor {

	//sampler := ring.NewGaussianSampler(prng, ringQ, sigmaSmudging, uint64(6*sigmaSmudging))

	return &mkDecryptor{
		params: params,
		ringQ:  GetRingQ(params),
	}

}

// PartDec computes a partial decription key for a given ciphertext.
func (dec *mkDecryptor) PartDec(ct *rlwe.Element, level uint64, sk *MKSecretKey, sigmaSmudging float64) *ring.Poly {

	if ct == nil || ct.Degree() < 1 {
		panic("Uninitialized Ciphertext")
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	sampler := ring.NewGaussianSampler(prng, dec.ringQ, sigmaSmudging, uint64(6*sigmaSmudging))

	if !ct.IsNTT {
		dec.ringQ.NTT(ct.Value[1], ct.Value[1])
	}

	// mu_i = c_i * sk_i + e_i mod q
	out := sampler.ReadLvlNew(level)
	dec.ringQ.NTTLvl(level, out, out)

	dec.ringQ.MulCoeffsMontgomeryAndAddLvl(level, ct.Value[1], sk.Key.Value, out)

	out.Coeffs = out.Coeffs[:level+1]

	return out
}

// MergeDec merges the partial decription parts and returns the plaintext.
// the same ciphertext that was used for PartDec must be provided
func (dec *mkDecryptor) MergeDec(ct *rlwe.Element, level uint64, partialKeys []*ring.Poly) *ring.Poly {

	if ct == nil || ct.Degree() < 1 {
		panic("Uninitilaized polynomial")
	}

	res := dec.ringQ.NewPoly()

	if !ct.IsNTT {
		dec.ringQ.NTT(ct.Value[0], res)
	} else {
		dec.ringQ.CopyLvl(level, ct.Value[0], res)
	}

	for _, k := range partialKeys {
		dec.ringQ.AddLvl(level, res, k, res)
	}

	dec.ringQ.ReduceLvl(level, res, res)
	res.Coeffs = res.Coeffs[:level+1]

	if !ct.IsNTT {
		dec.ringQ.InvNTT(res, res)
	}

	return res
}
