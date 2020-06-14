package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	context *dbfvContext

	sigmaSmudging float64

	tmp       *ring.Poly
	share0tmp *ring.Poly
	share1tmp *ring.Poly

	baseconverter   *ring.FastBasisExtender
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
}

// PCKSShare is a type for the PCKS protocol shares.
type PCKSShare [2]*ring.Poly

//type PCKSShare struct {
//	share [2]*ring.Poly
//}

// MarshalBinary encodes a PCKS share on a slice of bytes.
func (share *PCKSShare) MarshalBinary() ([]byte, error) {
	//TODO discuss choice here.
	//Maybe not worth it to have the metadata separated. we "lose" two bytes but complexity of the code would be higher in Unmarshalling.
	lenR1 := share[0].GetDataLen(true)
	lenR2 := share[1].GetDataLen(true)

	data := make([]byte, lenR1+lenR2)
	_, err := share[0].WriteTo(data[0:lenR1])
	if err != nil {
		return []byte{}, err
	}

	_, err = share[1].WriteTo(data[lenR1 : lenR1+lenR2])
	if err != nil {
		return []byte{}, err
	}

	return data, nil
}

// UnmarshalBinary decodes marshaled PCKS share on the target PCKS share.
func (share *PCKSShare) UnmarshalBinary(data []byte) error {

	if share[0] == nil {
		share[0] = new(ring.Poly)
	}

	if share[1] == nil {
		share[1] = new(ring.Poly)
	}

	err := share[0].UnmarshalBinary(data[0 : len(data)/2])
	if err != nil {
		return err
	}

	err = share[1].UnmarshalBinary(data[len(data)/2:])
	if err != nil {
		return err
	}

	return nil
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params *bfv.Parameters, sigmaSmudging float64) *PCKSProtocol {

	if !params.IsValid() {
		panic("cannot NewPCKSProtocol : params not valid (check if they where generated properly)")
	}

	context := newDbfvContext(params)

	pcks := new(PCKSProtocol)

	pcks.context = context

	pcks.sigmaSmudging = sigmaSmudging

	pcks.tmp = context.contextQP.NewPoly()
	pcks.share0tmp = context.contextQP.NewPoly()
	pcks.share1tmp = context.contextQP.NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(context.contextQ, context.contextP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pcks.gaussianSampler = ring.NewGaussianSamplerLvl(prng, context.contextQP)
	pcks.ternarySampler = ring.NewTernarySampler(prng, context.contextQP, 0.5, true)

	return pcks
}

// AllocateShares allocates the shares of the PCKS protocol
func (pcks *PCKSProtocol) AllocateShares() (s PCKSShare) {
	s[0] = pcks.context.contextQ.NewPoly()
	s[1] = pcks.context.contextQ.NewPoly()
	return
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + (u_i * pk[0] + e_0i)/P, (u_i * pk[1] + e_1i)/P]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *ring.Poly, pk *bfv.PublicKey, ct *bfv.Ciphertext, shareOut PCKSShare) {

	contextQ := pcks.context.contextQ
	contextKeys := pcks.context.contextQP

	pcks.ternarySampler.ReadNTT(pcks.tmp)

	// h_0 = u_i * pk_0
	contextKeys.MulCoeffsMontgomery(pcks.tmp, pk.Get()[0], pcks.share0tmp)
	// h_1 = u_i * pk_1
	contextKeys.MulCoeffsMontgomery(pcks.tmp, pk.Get()[1], pcks.share1tmp)

	contextKeys.InvNTT(pcks.share0tmp, pcks.share0tmp)
	contextKeys.InvNTT(pcks.share1tmp, pcks.share1tmp)

	// h_0 = u_i * pk_0 + e0
	pcks.gaussianSampler.Read(uint64(len(contextKeys.Modulus)-1), pcks.share0tmp, pcks.sigmaSmudging, uint64(6*pcks.sigmaSmudging))

	// h_1 = u_i * pk_1 + e1
	pcks.gaussianSampler.Read(uint64(len(contextKeys.Modulus)-1), pcks.share1tmp, pcks.sigmaSmudging, uint64(6*pcks.sigmaSmudging))

	// h_0 = (u_i * pk_0 + e0)/P
	pcks.baseconverter.ModDownPQ(uint64(len(contextQ.Modulus))-1, pcks.share0tmp, shareOut[0])

	// h_0 = (u_i * pk_0 + e0)/P
	// Could be moved to the keyswitch phase, but the second element of the shares will be larger
	pcks.baseconverter.ModDownPQ(uint64(len(contextQ.Modulus))-1, pcks.share1tmp, shareOut[1])

	// tmp = s_i*c_1
	contextQ.NTT(ct.Value()[1], pcks.tmp)
	contextQ.MulCoeffsMontgomery(pcks.tmp, sk, pcks.tmp)
	contextQ.InvNTT(pcks.tmp, pcks.tmp)

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	contextQ.Add(shareOut[0], pcks.tmp, shareOut[0])

	pcks.tmp.Zero()

}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {

	pcks.context.contextQ.Add(share1[0], share2[0], shareOut[0])
	pcks.context.contextQ.Add(share1[1], share2[1], shareOut[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *bfv.Ciphertext) {

	pcks.context.contextQ.Add(ct.Value()[0], combined[0], ctOut.Value()[0])
	pcks.context.contextQ.Copy(combined[1], ctOut.Value()[1])
}
