package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	bfvContext *bfv.BfvContext

	sigmaSmudging         float64
	gaussianSamplerSmudge *ring.KYSampler

	tmp  *ring.Poly
	tmpP *ring.Poly

	baseconverter *ring.FastBasisExtender
}

type PCKSShare [2]*ring.Poly

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(bfvContext *bfv.BfvContext, sigmaSmudging float64) *PCKSProtocol {

	pcks := new(PCKSProtocol)

	pcks.bfvContext = bfvContext

	pcks.gaussianSamplerSmudge = bfvContext.ContextKeys().NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	pcks.tmp = bfvContext.ContextKeys().NewPoly()
	pcks.tmpP = bfvContext.ContextPKeys().NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(bfvContext.ContextQ().Modulus, bfvContext.KeySwitchPrimes())

	return pcks
}

func (pcks *PCKSProtocol) AllocateShares() (s PCKSShare) {
	s[0] = pcks.bfvContext.ContextQ().NewPoly()
	s[1] = pcks.bfvContext.ContextQ().NewPoly()
	return
}

// GenShareRoundThree is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + u_i * pk[0] + e_0i, u_i * pk[1] + e_1i]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *ring.Poly, pk *bfv.PublicKey, ct *bfv.Ciphertext, shareOut PCKSShare) {

	contextQ := pcks.bfvContext.ContextQ()
	contextP := pcks.bfvContext.ContextPKeys()

	level := uint64(len(ct.Value()[1].Coeffs) - 1)

	//u_i
	_ = pcks.bfvContext.TernarySampler().SampleMontgomeryNTT(0.5, pcks.tmp)

	// h_0 = u_i * pk_0 (NTT)
	contextQ.MulCoeffsMontgomery(pcks.tmp, pk.Get()[0], shareOut[0])
	// h_1 = u_i * pk_1 (NTT)
	contextQ.MulCoeffsMontgomery(pcks.tmp, pk.Get()[1], shareOut[1])

	// h0 = u_i * pk_0 + s_i*c_1 (NTT)

	contextQ.NTT(ct.Value()[1], pcks.tmp)
	contextQ.MulCoeffsMontgomeryAndAdd(sk, pcks.tmp, shareOut[0])

	for _, pj := range pcks.bfvContext.KeySwitchPrimes() {
		contextQ.MulScalar(shareOut[0], pj, shareOut[0])
		contextQ.MulScalar(shareOut[1], pj, shareOut[1])
	}

	contextQ.InvNTT(shareOut[0], shareOut[0])
	contextQ.InvNTT(shareOut[1], shareOut[1])

	// h_0 = InvNTT(s_i*c_1 + u_i * pk_0) + e0
	pcks.gaussianSamplerSmudge.Sample(pcks.tmp)
	contextQ.Add(shareOut[0], pcks.tmp, shareOut[0])

	for x, i := 0, uint64(len(contextQ.Modulus)); i < uint64(len(pcks.bfvContext.ContextKeys().Modulus)); x, i = x+1, i+1 {
		for j := uint64(0); j < contextP.N; j++ {
			pcks.tmpP.Coeffs[x][j] = pcks.tmp.Coeffs[i][j]
		}
	}

	pcks.baseconverter.ModDownSplited(contextQ, contextP, pcks.bfvContext.RescaleParamsKeys(), level, shareOut[0], pcks.tmpP, shareOut[0], pcks.tmp)

	// h_1 = u_i * pk_1 + e1
	pcks.bfvContext.GaussianSampler().Sample(pcks.tmp)
	contextQ.Add(shareOut[1], pcks.tmp, shareOut[1])

	for x, i := 0, uint64(len(contextQ.Modulus)); i < uint64(len(pcks.bfvContext.ContextKeys().Modulus)); x, i = x+1, i+1 {
		for j := uint64(0); j < contextP.N; j++ {
			pcks.tmpP.Coeffs[x][j] = pcks.tmp.Coeffs[i][j]
		}
	}

	pcks.baseconverter.ModDownSplited(contextQ, contextP, pcks.bfvContext.RescaleParamsKeys(), level, shareOut[1], pcks.tmpP, shareOut[1], pcks.tmp)

	pcks.tmpP.Zero()
	pcks.tmp.Zero()
}

// GenShareRoundTwo is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {

	pcks.bfvContext.ContextQ().Add(share1[0], share2[0], shareOut[0])
	pcks.bfvContext.ContextQ().Add(share1[1], share2[1], shareOut[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *bfv.Ciphertext) {

	pcks.bfvContext.ContextQ().Add(ct.Value()[0], combined[0], ctOut.Value()[0])
	pcks.bfvContext.ContextQ().Copy(combined[1], ctOut.Value()[1])
}
