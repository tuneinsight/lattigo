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

	polypool *ring.Poly

	baseconverter *ring.FastBasisExtender
}

type PCKSShare [2]*ring.Poly

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(ckksContext *ckks.CkksContext, sigmaSmudging float64) *PCKSProtocol {

	pcks := new(PCKSProtocol)

	pcks.ckksContext = ckksContext

	pcks.gaussianSamplerSmudge = pcks.ckksContext.ContextKey(ckksContext.Levels()-1).NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	pcks.tmp = pcks.ckksContext.ContextKey(ckksContext.Levels() - 1).NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(ckksContext.ContextQ().Modulus, ckksContext.KeySwitchPrimes())

	return pcks
}

func (pcks *PCKSProtocol) AllocateShares(level uint64) (s PCKSShare) {
	s[0] = pcks.ckksContext.ContextKey(level).NewPoly()
	s[1] = pcks.ckksContext.ContextKey(level).NewPoly()
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

	//u_i
	_ = pcks.ckksContext.TernarySampler().SampleMontgomeryNTT(0.5, pcks.tmp)

	// h_0 = u_i * pk_0
	contextQ.MulCoeffsMontgomery(pcks.tmp, pk.Get()[0], shareOut[0])

	// h_1 = u_i * pk_1
	contextQ.MulCoeffsMontgomery(pcks.tmp, pk.Get()[1], shareOut[1])

	contextQ.Copy(ct.Value()[1], pcks.tmp)

	contextQ.InvNTT(pcks.tmp, pcks.tmp)
	pcks.baseconverter.ModUp(ct.Level(), pcks.tmp, pcks.tmp)
	contextKeys.NTT(pcks.tmp, pcks.tmp)

	// h0 = u_i * pk_0 + s_i*c_1
	contextKeys.MulCoeffsMontgomeryAndAdd(sk, pcks.tmp, shareOut[0])

	for _, pj := range pcks.ckksContext.KeySwitchPrimes() {
		contextKeys.MulScalar(shareOut[0], pj, shareOut[0])
		contextKeys.MulScalar(shareOut[1], pj, shareOut[1])
	}
	// h_0 = s_i*c_1 + u_i * pk_0 + e0
	pcks.gaussianSamplerSmudge.SampleNTT(pcks.tmp)
	contextKeys.Add(shareOut[0], pcks.tmp, shareOut[0])

	// h_1 = u_i * pk_1 + e1
	pcks.ckksContext.GaussianSampler().SampleNTT(pcks.tmp)
	contextKeys.Add(shareOut[1], pcks.tmp, shareOut[1])

	pcks.baseconverter.ModDownNTT(contextQ, contextP, pcks.ckksContext.RescaleParamsKeys(), ct.Level(), shareOut[0], shareOut[0], pcks.tmp)
	pcks.baseconverter.ModDownNTT(contextQ, contextP, pcks.ckksContext.RescaleParamsKeys(), ct.Level(), shareOut[1], shareOut[1], pcks.tmp)

	shareOut[0].Coeffs = shareOut[0].Coeffs[:ct.Level()+1]
	shareOut[1].Coeffs = shareOut[1].Coeffs[:ct.Level()+1]

	pcks.tmp.Zero()
}

// GenShareRoundTwo is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {
	pcks.ckksContext.Context(uint64(len(share1[0].Coeffs))-1).Add(share1[0], share2[0], shareOut[0])
	pcks.ckksContext.Context(uint64(len(share1[0].Coeffs))-1).Add(share1[1], share2[1], shareOut[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *ckks.Ciphertext) {

	pcks.ckksContext.Context(ct.Level()).Add(ct.Value()[0], combined[0], ctOut.Value()[0])
	ctOut.Value()[1].Copy(combined[1])
}
