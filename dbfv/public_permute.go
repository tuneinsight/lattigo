package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// PermuteProtocol is a struct storing the parameters for the PermuteProtocol protocol.
type PermuteProtocol struct {
	context         *dbfvContext
	indexMatrix     []uint64
	tmp1            *ring.Poly
	tmp2            *ring.Poly
	hP              *ring.Poly
	baseconverter   *ring.FastBasisExtender
	scaler          ring.Scaler
	gaussianSampler *ring.GaussianSampler
	sigma           float64
	uniformSampler  *ring.UniformSampler
}

// NewPermuteProtocol creates a new instance of the PermuteProtocol.
func NewPermuteProtocol(params *bfv.Parameters) (refreshProtocol *PermuteProtocol) {

	context := newDbfvContext(params)

	refreshProtocol = new(PermuteProtocol)
	refreshProtocol.context = context
	refreshProtocol.tmp1 = context.ringQP.NewPoly()
	refreshProtocol.tmp2 = context.ringQP.NewPoly()
	refreshProtocol.hP = context.ringP.NewPoly()

	refreshProtocol.baseconverter = ring.NewFastBasisExtender(context.ringQ, context.ringP)

	var m, pos, index1, index2 int

	indexMatrix := make([]uint64, params.N())

	logN := params.LogN()

	rowSize := params.N() >> 1
	m = (params.N() << 1)
	pos = 1

	for i := 0; i < rowSize; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		indexMatrix[i] = utils.BitReverse64(uint64(index1), uint64(logN))
		indexMatrix[i|rowSize] = utils.BitReverse64(uint64(index2), uint64(logN))

		pos *= bfv.GaloisGen
		pos &= (m - 1)
	}

	refreshProtocol.indexMatrix = indexMatrix
	refreshProtocol.scaler = ring.NewRNSScaler(params.T(), context.ringQ)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	refreshProtocol.gaussianSampler = ring.NewGaussianSampler(prng)
	refreshProtocol.sigma = params.Sigma()
	refreshProtocol.uniformSampler = ring.NewUniformSampler(prng, context.ringT)

	return
}

// AllocateShares allocates the shares of the PermuteProtocol
func (pp *PermuteProtocol) AllocateShares() RefreshShare {
	return RefreshShare{pp.context.ringQ.NewPoly(),
		pp.context.ringQ.NewPoly()}
}

// GenShares generates the shares of the PermuteProtocol
func (pp *PermuteProtocol) GenShares(sk *ring.Poly, ciphertext *bfv.Ciphertext, crs *ring.Poly, permutation []uint64, share RefreshShare) {

	level := len(ciphertext.Value()[1].Coeffs) - 1

	ringQ := pp.context.ringQ
	ringT := pp.context.ringT
	ringQP := pp.context.ringQP

	// h0 = s*ct[1]
	ringQ.NTTLazy(ciphertext.Value()[1], pp.tmp1)
	ringQ.MulCoeffsMontgomeryConstant(sk, pp.tmp1, share.RefreshShareDecrypt)
	ringQ.InvNTTLazy(share.RefreshShareDecrypt, share.RefreshShareDecrypt)

	// h0 = s*ct[1]*P
	ringQ.MulScalarBigint(share.RefreshShareDecrypt, pp.context.ringP.ModulusBigint, share.RefreshShareDecrypt)

	// h0 = s*ct[1]*P + e
	pp.gaussianSampler.ReadLvl(len(ringQP.Modulus)-1, pp.tmp1, ringQP, pp.sigma, int(6*pp.sigma))
	ringQ.Add(share.RefreshShareDecrypt, pp.tmp1, share.RefreshShareDecrypt)

	for x, i := 0, len(ringQ.Modulus); i < len(pp.context.ringQP.Modulus); x, i = x+1, i+1 {
		tmphP := pp.hP.Coeffs[x]
		tmp1 := pp.tmp1.Coeffs[i]
		for j := 0; j < ringQ.N; j++ {
			tmphP[j] += tmp1[j]
		}
	}

	// h0 = (s*ct[1]*P + e)/P
	pp.baseconverter.ModDownSplitPQ(level, share.RefreshShareDecrypt, pp.hP, share.RefreshShareDecrypt)

	// h1 = -s*a
	ringQP.Neg(crs, pp.tmp1)
	ringQP.NTTLazy(pp.tmp1, pp.tmp1)
	ringQP.MulCoeffsMontgomeryConstant(sk, pp.tmp1, pp.tmp2)
	ringQP.InvNTTLazy(pp.tmp2, pp.tmp2)

	// h1 = s*a + e'
	pp.gaussianSampler.ReadAndAdd(pp.tmp2, ringQP, pp.sigma, int(6*pp.sigma))

	// h1 = (-s*a + e')/P
	pp.baseconverter.ModDownPQ(level, pp.tmp2, share.RefreshShareRecrypt)

	// mask = (uniform plaintext in [0, T-1]) * floor(Q/T)

	// Mask in the time domain
	coeffs := pp.uniformSampler.ReadNew()

	// Multiply by Q/t
	lift(coeffs, pp.tmp1, pp.context)

	// h0 = (s*ct[1]*P + e)/P + mask
	ringQ.Add(share.RefreshShareDecrypt, pp.tmp1, share.RefreshShareDecrypt)

	// Mask in the spectral domain
	ringT.NTT(coeffs, coeffs)

	// Permutation over the mask
	pp.permuteWithIndex(coeffs, permutation, pp.tmp1)

	// Switch back the mask in the time domain
	ringT.InvNTTLazy(pp.tmp1, coeffs)

	// Multiply by Q/t
	lift(coeffs, pp.tmp1, pp.context)

	// h1 = (-s*a + e')/P - permute(mask)
	ringQ.Sub(share.RefreshShareRecrypt, pp.tmp1, share.RefreshShareRecrypt)
}

// Aggregate sums share1 and share2 on shareOut.
func (pp *PermuteProtocol) Aggregate(share1, share2, shareOut RefreshShare) {
	pp.context.ringQ.Add(share1.RefreshShareDecrypt, share2.RefreshShareDecrypt, shareOut.RefreshShareDecrypt)
	pp.context.ringQ.Add(share1.RefreshShareRecrypt, share2.RefreshShareRecrypt, shareOut.RefreshShareRecrypt)
}

// Decrypt operates a masked decryption on the input ciphertext using the provided decryption shares.
func (pp *PermuteProtocol) Decrypt(ciphertext *bfv.Ciphertext, shareDecrypt RefreshShareDecrypt, sharePlaintext *ring.Poly) {
	pp.context.ringQ.Add(ciphertext.Value()[0], shareDecrypt, sharePlaintext)
}

// Permute decodes and re-encode (removing the error) the masked decrypted ciphertext with a permutation of the plaintext slots.
func (pp *PermuteProtocol) Permute(sharePlaintext *ring.Poly, permutation []uint64, sharePlaintextOut *ring.Poly) {

	ringT := pp.context.ringT

	pp.scaler.DivByQOverTRounded(sharePlaintext, sharePlaintextOut)

	ringT.NTT(sharePlaintextOut, sharePlaintextOut)

	pp.permuteWithIndex(sharePlaintextOut, permutation, pp.tmp1)

	ringT.InvNTTLazy(pp.tmp1, sharePlaintextOut)

	lift(sharePlaintextOut, sharePlaintextOut, pp.context)
}

// Recrypt recrypts the input masked decrypted ciphertext with the recryption shares.
func (pp *PermuteProtocol) Recrypt(sharePlaintext *ring.Poly, crs *ring.Poly, shareRecrypt RefreshShareRecrypt, ciphertextOut *bfv.Ciphertext) {

	// ciphertext[0] = (-crs*s + e')/P + permute(m)
	pp.context.ringQ.Add(sharePlaintext, shareRecrypt, ciphertextOut.Value()[0])

	// ciphertext[1] = crs/P
	pp.baseconverter.ModDownPQ(len(ciphertextOut.Value()[1].Coeffs)-1, crs, ciphertextOut.Value()[1])

}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (pp *PermuteProtocol) Finalize(ciphertext *bfv.Ciphertext, permutation []uint64, crs *ring.Poly, share RefreshShare, ciphertextOut *bfv.Ciphertext) {
	pp.Decrypt(ciphertext, share.RefreshShareDecrypt, pp.tmp1)
	pp.Permute(pp.tmp1, permutation, pp.tmp1)
	pp.Recrypt(pp.tmp1, crs, share.RefreshShareRecrypt, ciphertextOut)
}

func (pp *PermuteProtocol) permuteWithIndex(polIn *ring.Poly, index []uint64, polOut *ring.Poly) {
	for j := 0; j < len(polIn.Coeffs[0]); j++ {
		polOut.Coeffs[0][pp.indexMatrix[j]] = polIn.Coeffs[0][pp.indexMatrix[index[j]]]
	}
}
