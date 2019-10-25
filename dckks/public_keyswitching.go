package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	ckksContext *ckks.CkksContext

	sigmaSmudging         float64
	gaussianSamplerSmudge *ring.KYSampler

	tmp *ring.Poly

	baseconverter *ring.FastBasisExtender
}

type PCKSShare [2]*ring.Poly

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(ckksContext *ckks.CkksContext, sigmaSmudging float64) *PCKSProtocol {

	pcks := new(PCKSProtocol)

	pcks.ckksContext = ckksContext

	pcks.gaussianSamplerSmudge = pcks.ckksContext.ContextKeys().NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	pcks.tmp = pcks.ckksContext.ContextKeys().NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(ckksContext.ContextQ().Modulus, ckksContext.KeySwitchPrimes())

	return pcks
}

func (pcks *PCKSProtocol) AllocateShares(level uint64) (s PCKSShare) {
	s[0] = pcks.ckksContext.ContextKeys().NewPolyLvl(level + uint64(len(pcks.ckksContext.KeySwitchPrimes())))
	s[1] = pcks.ckksContext.ContextKeys().NewPolyLvl(level + uint64(len(pcks.ckksContext.KeySwitchPrimes())))
	return
}

// GenShareRoundThree is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + u_i * pk[0] + e_0i, u_i * pk[1] + e_1i]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *ring.Poly, pk *ckks.PublicKey, ct *ckks.Ciphertext, shareOut PCKSShare) {

	contextQ := pcks.ckksContext.ContextQ()
	contextP := pcks.ckksContext.ContextP()

	//u_i
	_ = pcks.ckksContext.TernarySampler().SampleMontgomeryNTT(0.5, pcks.tmp)

	// h_0 = u_i * pk_0
	contextQ.MulCoeffsMontgomeryLvl(ct.Level(), pcks.tmp, pk.Get()[0], shareOut[0])

	// h_1 = u_i * pk_1
	contextQ.MulCoeffsMontgomeryLvl(ct.Level(), pcks.tmp, pk.Get()[1], shareOut[1])

	// h0 = u_i * pk_0 + s_i*c_1
	contextQ.MulCoeffsMontgomeryAndAddLvl(ct.Level(), sk, ct.Value()[1], shareOut[0])

	// TODO : improve by pre-computing prd(pj) for each qi
	for _, pj := range pcks.ckksContext.KeySwitchPrimes() {
		contextQ.MulScalarLvl(ct.Level(), shareOut[0], pj, shareOut[0])
		contextQ.MulScalarLvl(ct.Level(), shareOut[1], pj, shareOut[1])
	}

	share0P := contextP.NewPoly()

	// h_0 = s_i*c_1 + u_i * pk_0 + e0
	pcks.gaussianSamplerSmudge.SampleNTT(pcks.tmp)
	contextQ.Add(shareOut[0], pcks.tmp, shareOut[0])

	for x, i := 0, uint64(len(contextQ.Modulus)); i < uint64(len(pcks.ckksContext.ContextKeys().Modulus)); x, i = x+1, i+1 {
		for j := uint64(0); j < contextP.N; j++ {
			share0P.Coeffs[x][j] += pcks.tmp.Coeffs[i][j]
		}
	}

	share1P := contextP.NewPoly()

	// h_1 = u_i * pk_1 + e1
	pcks.ckksContext.GaussianSampler().SampleNTT(pcks.tmp)
	contextQ.Add(shareOut[1], pcks.tmp, shareOut[1])

	for x, i := 0, uint64(len(contextQ.Modulus)); i < uint64(len(pcks.ckksContext.ContextKeys().Modulus)); x, i = x+1, i+1 {
		for j := uint64(0); j < contextP.N; j++ {
			share1P.Coeffs[x][j] += pcks.tmp.Coeffs[i][j]
		}
	}

	pcks.baseconverter.ModDownSplitedNTT(contextQ, contextP, pcks.ckksContext.RescaleParamsKeys(), ct.Level(), shareOut[0], share0P, shareOut[0], pcks.tmp)
	pcks.baseconverter.ModDownSplitedNTT(contextQ, contextP, pcks.ckksContext.RescaleParamsKeys(), ct.Level(), shareOut[1], share1P, shareOut[1], pcks.tmp)

	shareOut[0].Coeffs = shareOut[0].Coeffs[:ct.Level()+1]
	shareOut[1].Coeffs = shareOut[1].Coeffs[:ct.Level()+1]

	pcks.tmp.Zero()
}

// GenShareRoundTwo is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {

	level := uint64(len(share1[0].Coeffs)) - 1
	pcks.ckksContext.ContextQ().AddLvl(level, share1[0], share2[0], shareOut[0])
	pcks.ckksContext.ContextQ().AddLvl(level, share1[1], share2[1], shareOut[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *ckks.Ciphertext) {

	pcks.ckksContext.ContextQ().AddLvl(ct.Level(), ct.Value()[0], combined[0], ctOut.Value()[0])
	pcks.ckksContext.ContextQ().CopyLvl(ct.Level(), combined[1], ctOut.Value()[1])
}
