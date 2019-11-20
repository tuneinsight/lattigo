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

	share0tmp *ring.Poly
	share1tmp *ring.Poly

	baseconverter *ring.FastBasisExtender
}

type PCKSShare [2]*ring.Poly

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(ckksContext *ckks.CkksContext, sigmaSmudging float64) *PCKSProtocol {

	pcks := new(PCKSProtocol)

	pcks.ckksContext = ckksContext

	pcks.gaussianSamplerSmudge = ckksContext.ContextKeys().NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	pcks.tmp = ckksContext.ContextKeys().NewPoly()
	pcks.share0tmp = ckksContext.ContextKeys().NewPoly()
	pcks.share1tmp = ckksContext.ContextKeys().NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(ckksContext.ContextQ().Modulus, ckksContext.KeySwitchPrimes())

	return pcks
}

func (pcks *PCKSProtocol) AllocateShares(level uint64) (s PCKSShare) {
	s[0] = pcks.ckksContext.ContextQ().NewPolyLvl(level)
	s[1] = pcks.ckksContext.ContextQ().NewPolyLvl(level)
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
	contextKeys := pcks.ckksContext.ContextKeys()

	contextKeys.SampleTernaryMontgomeryNTT(pcks.tmp, 0.5)

	// h_0 = u_i * pk_0
	contextKeys.MulCoeffsMontgomery(pcks.tmp, pk.Get()[0], pcks.share0tmp)
	// h_1 = u_i * pk_1
	contextKeys.MulCoeffsMontgomery(pcks.tmp, pk.Get()[1], pcks.share1tmp)

	// h_0 = u_i * pk_0 + e0
	pcks.gaussianSamplerSmudge.SampleNTT(pcks.tmp)
	contextKeys.Add(pcks.share0tmp, pcks.tmp, pcks.share0tmp)
	// h_1 = u_i * pk_1 + e1
	pcks.ckksContext.GaussianSampler().SampleNTT(pcks.tmp)
	contextKeys.Add(pcks.share1tmp, pcks.tmp, pcks.share1tmp)

	// h_0 = (u_i * pk_0 + e0)/P
	pcks.baseconverter.ModDownNTT(contextQ, contextP, pcks.ckksContext.RescaleParamsKeys(), ct.Level(), pcks.share0tmp, shareOut[0], pcks.tmp)

	// h_1 = (u_i * pk_1 + e1)/P
	// Cound be moved to the keyswitch part of the protocol, but the second element of the shares will be larger.
	pcks.baseconverter.ModDownNTT(contextQ, contextP, pcks.ckksContext.RescaleParamsKeys(), ct.Level(), pcks.share1tmp, shareOut[1], pcks.tmp)

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	contextQ.MulCoeffsMontgomeryAndAddLvl(ct.Level(), ct.Value()[1], sk, shareOut[0])

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
