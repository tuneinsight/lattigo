package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKSProtocol is the structure storing the parameters and and precomputations for the collective key-switching protocol.
type CKSProtocol struct {
	ringQ           *ring.Ring
	ringP           *ring.Ring
	ringQP          *ring.Ring
	gaussianSampler *ring.GaussianSampler
	baseconverter   *ring.FastBasisExtender

	tmpP     *ring.Poly
	tmpQ     *ring.Poly
	tmpDelta *ring.Poly
	tmpNtt   *ring.Poly
	hP       *ring.Poly
}

// CKSShare is a type for the CKS protocol shares.
type CKSShare struct {
	Value *ring.Poly
}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) *CKSProtocol {
	cks := new(CKSProtocol)
	var err error
	cks.ringQ, err = ring.NewRing(params.N(), params.Q())
	if err != nil {
		panic(err)
	}
	cks.ringP, err = ring.NewRing(params.N(), params.P())
	if err != nil {
		panic(err)
	}
	cks.ringQP, err = ring.NewRing(params.N(), params.QP())
	if err != nil {
		panic(err)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	cks.gaussianSampler = ring.NewGaussianSampler(prng, cks.ringQP, sigmaSmudging, int(6*sigmaSmudging))

	cks.baseconverter = ring.NewFastBasisExtender(cks.ringQ, cks.ringP)

	cks.tmpQ = cks.ringQ.NewPoly()
	cks.tmpP = cks.ringP.NewPoly()
	cks.tmpDelta = cks.ringQ.NewPoly()
	cks.tmpNtt = cks.ringQP.NewPoly()
	cks.hP = cks.ringP.NewPoly()

	return cks
}

// AllocateShare allocates the shares of the CKSProtocol
func (cks *CKSProtocol) AllocateShare() CKSShare {
	return CKSShare{cks.ringQ.NewPoly()}
}

// GenShare generates a party's share in the CKSProtocol
func (cks *CKSProtocol) GenShare(skInput, skOutput *ring.Poly, ct *rlwe.Element, shareOut CKSShare) {

	cks.ringQ.Sub(skInput, skOutput, cks.tmpDelta)

	cks.genShareDeltaBFV(cks.tmpDelta, ct, shareOut)
}

func (cks *CKSProtocol) genShareDeltaBFV(skDelta *ring.Poly, ct *rlwe.Element, shareOut CKSShare) { // BFV

	level := len(ct.Value[1].Coeffs) - 1
	cks.ringQ.NTTLazy(ct.Value[1], cks.tmpNtt)

	cks.ringQ.MulCoeffsMontgomeryConstant(cks.tmpNtt, skDelta, shareOut.Value)
	cks.ringQ.MulScalarBigint(shareOut.Value, cks.ringP.ModulusBigint, shareOut.Value)

	cks.ringQ.InvNTTLazy(shareOut.Value, shareOut.Value)

	cks.gaussianSampler.ReadLvl(len(cks.ringQP.Modulus)-1, cks.tmpNtt)
	cks.ringQ.AddNoMod(shareOut.Value, cks.tmpNtt, shareOut.Value)

	for x, i := 0, uint64(len(cks.ringQ.Modulus)); i < uint64(len(cks.ringQP.Modulus)); x, i = x+1, i+1 {
		tmphP := cks.hP.Coeffs[x]
		tmpNTT := cks.tmpNtt.Coeffs[i]
		for j := 0; j < cks.ringQ.N; j++ {
			tmphP[j] += tmpNTT[j]
		}
	}

	cks.baseconverter.ModDownSplitPQ(level, shareOut.Value, cks.hP, shareOut.Value)

	cks.tmpNtt.Zero()
	cks.hP.Zero()
}

func (cks *CKSProtocol) genShareDeltaCKKS(skDelta *ring.Poly, ct *rlwe.Element, shareOut CKSShare) { // CKKS

	cks.ringQ.MulCoeffsMontgomeryConstantLvl(ct.Level(), ct.Value[1], skDelta, shareOut.Value)
	cks.ringQ.MulScalarBigintLvl(ct.Level(), shareOut.Value, cks.ringP.ModulusBigint, shareOut.Value)

	cks.gaussianSampler.ReadLvl(ct.Level(), cks.tmpQ)
	extendBasisSmallNormAndCenter(cks.ringQ, cks.ringP, cks.tmpQ, cks.tmpP)

	cks.ringQ.NTTLvl(ct.Level(), cks.tmpQ, cks.tmpQ)
	cks.ringP.NTT(cks.tmpP, cks.tmpP)

	cks.ringQ.AddLvl(ct.Level(), shareOut.Value, cks.tmpQ, shareOut.Value)

	cks.baseconverter.ModDownSplitNTTPQ(ct.Level(), shareOut.Value, cks.tmpP, shareOut.Value)

	cks.tmpQ.Zero()
	cks.tmpP.Zero()
}

func extendBasisSmallNormAndCenter(ringQ, ringP *ring.Ring, polQ, polP *ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = ringQ.Modulus[0]
	QHalf = Q >> 1

	for j := 0; j < ringQ.N; j++ {

		coeff = polQ.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i, pi := range ringP.Modulus {
			polP.Coeffs[i][j] = (coeff * sign) | (pi-coeff)*(sign^1)
		}
	}
}
