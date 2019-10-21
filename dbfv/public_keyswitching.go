package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	contextCiphertexts *ring.Context
	contextKeys        *ring.Context

	sigmaSmudging         float64
	gaussianSamplerSmudge *ring.KYSampler
	gaussianSampler       *ring.KYSampler
	ternarySampler        *ring.TernarySampler

	tmp *ring.Poly

	polypool *ring.Poly

	baseconverter *bfv.FastBasisExtender

	rescaleParamsKeys []uint64
	keyswitchprimes   []uint64
	rescalepool       []uint64
}

type PCKSShare [2]*ring.Poly

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(bfvContext *bfv.BfvContext, sigmaSmudging float64) *PCKSProtocol {

	pcks := new(PCKSProtocol)

	pcks.contextCiphertexts = bfvContext.ContextQ()
	pcks.contextKeys = bfvContext.ContextKeys()

	pcks.gaussianSamplerSmudge = pcks.contextKeys.NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))
	pcks.gaussianSampler = bfvContext.GaussianSampler()
	pcks.ternarySampler = bfvContext.TernarySampler()

	pcks.tmp = pcks.contextKeys.NewPoly()

	pcks.keyswitchprimes = make([]uint64, len(bfvContext.KeySwitchPrimes()))
	for i, pi := range bfvContext.KeySwitchPrimes() {
		pcks.keyswitchprimes[i] = pi
	}

	pcks.baseconverter = bfv.NewFastBasisExtender(pcks.contextCiphertexts.Modulus, pcks.keyswitchprimes)

	pcks.rescaleParamsKeys = make([]uint64, len(pcks.contextCiphertexts.Modulus))

	PBig := ring.NewUint(1)
	for _, pj := range pcks.keyswitchprimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := ring.NewUint(0)
	bredParams := pcks.contextCiphertexts.GetBredParams()
	for i, Qi := range pcks.contextCiphertexts.Modulus {
		tmp.Mod(PBig, ring.NewUint(Qi))
		pcks.rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	return pcks
}

func (pcks *PCKSProtocol) AllocateShares() (s PCKSShare) {
	s[0] = pcks.contextKeys.NewPoly()
	s[1] = pcks.contextKeys.NewPoly()
	return
}

// GenShareRoundThree is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + u_i * pk[0] + e_0i, u_i * pk[1] + e_1i]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *ring.Poly, pk *bfv.PublicKey, ct *bfv.Ciphertext, shareOut PCKSShare) {

	level := uint64(len(ct.Value()[1].Coeffs) - 1)

	//u_i
	_ = pcks.ternarySampler.SampleMontgomeryNTT(0.5, pcks.tmp)

	// h_0 = u_i * pk_0 (NTT)
	pcks.contextKeys.MulCoeffsMontgomery(pcks.tmp, pk.Get()[0], shareOut[0])
	// h_1 = u_i * pk_1 (NTT)
	pcks.contextKeys.MulCoeffsMontgomery(pcks.tmp, pk.Get()[1], shareOut[1])

	// h0 = u_i * pk_0 + s_i*c_1 (NTT)

	pcks.contextCiphertexts.Copy(ct.Value()[1], pcks.tmp)
	pcks.baseconverter.ModUp(level, pcks.tmp, pcks.tmp)

	pcks.contextKeys.NTT(pcks.tmp, pcks.tmp)
	pcks.contextKeys.MulCoeffsMontgomeryAndAdd(sk, pcks.tmp, shareOut[0])

	for _, pj := range pcks.keyswitchprimes {
		pcks.contextKeys.MulScalar(shareOut[0], pj, shareOut[0])
		pcks.contextKeys.MulScalar(shareOut[1], pj, shareOut[1])
	}

	// h_0 = InvNTT(s_i*c_1 + u_i * pk_0) + e0
	pcks.gaussianSamplerSmudge.SampleNTT(pcks.tmp)
	pcks.contextKeys.Add(shareOut[0], pcks.tmp, shareOut[0])

	// h_1 = InvNTT(u_i * pk_1) + e1
	pcks.gaussianSampler.SampleNTT(pcks.tmp)
	pcks.contextKeys.Add(shareOut[1], pcks.tmp, shareOut[1])

	pcks.baseconverter.ModDown(pcks.contextKeys, pcks.rescaleParamsKeys, level, shareOut[0], shareOut[0], pcks.tmp)
	pcks.baseconverter.ModDown(pcks.contextKeys, pcks.rescaleParamsKeys, level, shareOut[1], shareOut[1], pcks.tmp)

	pcks.contextCiphertexts.InvNTT(shareOut[0], shareOut[0])
	pcks.contextCiphertexts.InvNTT(shareOut[1], shareOut[1])

	shareOut[0].Coeffs = shareOut[0].Coeffs[:level+1]
	shareOut[1].Coeffs = shareOut[1].Coeffs[:level+1]

	pcks.tmp.Zero()
}

// GenShareRoundTwo is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {
	pcks.contextCiphertexts.Add(share1[0], share2[0], shareOut[0])
	pcks.contextCiphertexts.Add(share1[1], share2[1], shareOut[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *bfv.Ciphertext) {

	pcks.contextCiphertexts.Add(ct.Value()[0], combined[0], ctOut.Value()[0])
	ctOut.Value()[1].Copy(combined[1])
}
