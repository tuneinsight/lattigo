package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	dckksContext *dckksContext

	sigmaSmudging float64

	tmp *ring.Poly

	share0tmp *ring.Poly
	share1tmp *ring.Poly

	baseconverter *ring.FastBasisExtender
}

// PCKSShare is a struct storing the share of the PCKS protocol.
type PCKSShare [2]*ring.Poly

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(params *ckks.Parameters, sigmaSmudging float64) *PCKSProtocol {

	if !params.IsValid() {
		panic("cannot NewPCKSProtocol : params not valid (check if they where generated properly)")
	}

	pcks := new(PCKSProtocol)

	dckksContext := newDckksContext(params)

	pcks.dckksContext = dckksContext

	pcks.tmp = dckksContext.contextQP.NewPoly()
	pcks.share0tmp = dckksContext.contextQP.NewPoly()
	pcks.share1tmp = dckksContext.contextQP.NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(dckksContext.contextQ, dckksContext.contextP)

	return pcks
}

// AllocateShares allocates the share of the PCKS protocol.
func (pcks *PCKSProtocol) AllocateShares(level uint64) (s PCKSShare) {
	s[0] = pcks.dckksContext.contextQ.NewPolyLvl(level)
	s[1] = pcks.dckksContext.contextQ.NewPolyLvl(level)
	return
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + u_i * pk[0] + e_0i, u_i * pk[1] + e_1i]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *ring.Poly, pk *ckks.PublicKey, ct *ckks.Ciphertext, shareOut PCKSShare) {

	contextQ := pcks.dckksContext.contextQ
	contextKeys := pcks.dckksContext.contextQP
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	gaussianSampler := ring.NewGaussianSampler(prng, contextKeys)
	ternarySampler := ring.NewTernarySampler(prng, contextKeys)

	ternarySampler.SampleTernaryMontgomeryNTT(pcks.tmp, 0.5)

	// h_0 = u_i * pk_0
	contextKeys.MulCoeffsMontgomery(pcks.tmp, pk.Get()[0], pcks.share0tmp)
	// h_1 = u_i * pk_1
	contextKeys.MulCoeffsMontgomery(pcks.tmp, pk.Get()[1], pcks.share1tmp)

	// h_0 = u_i * pk_0 + e0
	gaussianSampler.SampleGaussianNTTLvl(uint64(len(contextKeys.Modulus)-1), pcks.tmp, pcks.sigmaSmudging, uint64(6*pcks.sigmaSmudging))
	contextKeys.Add(pcks.share0tmp, pcks.tmp, pcks.share0tmp)
	// h_1 = u_i * pk_1 + e1
	gaussianSampler.SampleGaussianNTTLvl(uint64(len(contextKeys.Modulus)-1), pcks.tmp, pcks.sigmaSmudging, uint64(6*pcks.sigmaSmudging))
	contextKeys.Add(pcks.share1tmp, pcks.tmp, pcks.share1tmp)

	// h_0 = (u_i * pk_0 + e0)/P
	pcks.baseconverter.ModDownNTTPQ(ct.Level(), pcks.share0tmp, shareOut[0])

	// h_1 = (u_i * pk_1 + e1)/P
	// Cound be moved to the keyswitch part of the protocol, but the second element of the shares will be larger.
	pcks.baseconverter.ModDownNTTPQ(ct.Level(), pcks.share1tmp, shareOut[1])

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	contextQ.MulCoeffsMontgomeryAndAddLvl(ct.Level(), ct.Value()[1], sk, shareOut[0])

	pcks.tmp.Zero()
}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {

	level := uint64(len(share1[0].Coeffs)) - 1
	pcks.dckksContext.contextQ.AddLvl(level, share1[0], share2[0], shareOut[0])
	pcks.dckksContext.contextQ.AddLvl(level, share1[1], share2[1], shareOut[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *ckks.Ciphertext) {

	ctOut.SetScale(ct.Scale())

	pcks.dckksContext.contextQ.AddLvl(ct.Level(), ct.Value()[0], combined[0], ctOut.Value()[0])
	pcks.dckksContext.contextQ.CopyLvl(ct.Level(), combined[1], ctOut.Value()[1])
}
