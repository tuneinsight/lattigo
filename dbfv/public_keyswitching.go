package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	ringContext *ring.Context

	sigmaSmudging         float64
	gaussianSamplerSmudge *ring.KYSampler
	gaussianSampler       *ring.KYSampler
	ternarySampler        *ring.TernarySampler

	tmp *ring.Poly
}

type PCKSShare struct {
	share [2]*ring.Poly
}
//TODO write marshalling and unmarshalling
func (share *PCKSShare) MarshalBinary() ([]byte,error){
	lenR1 := share.share[0].GetDataLen()
	lenR2 := share.share[1].GetDataLen()


	data := make([]byte,lenR1+lenR2)
	_,err := share.share[0].WriteTo(data[0:lenR1])
	if err != nil{
		return []byte{},err
	}

	_, err =share.share[1].WriteTo(data[lenR1:lenR1+lenR2])
	if err != nil{
		return []byte{},err
	}

	return data,nil
}

func (share *PCKSShare) UnmarshalBinary(data []byte )(error){
	if share == nil{
		share = new(PCKSShare)
	}

	share.share[0] = new(ring.Poly)
	share.share[1] = new(ring.Poly)

	err := share.share[0].UnmarshalBinary(data[0:len(data)/2])
	if err != nil{
		return err
	}

	err = share.share[1].UnmarshalBinary(data[len(data)/2:])
	if err != nil{
		return err
	}

	return nil
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(bfvContext *bfv.BfvContext, sigmaSmudging float64) *PCKSProtocol {

	pcks := new(PCKSProtocol)
	pcks.ringContext = bfvContext.ContextQ()
	pcks.gaussianSamplerSmudge = pcks.ringContext.NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))
	pcks.gaussianSampler = bfvContext.GaussianSampler()
	pcks.ternarySampler = bfvContext.TernarySampler()

	pcks.tmp = pcks.ringContext.NewPoly()
	return pcks
}

func (pcks *PCKSProtocol) AllocateShares() (s PCKSShare) {
	s.share[0] = pcks.ringContext.NewPoly()
	s.share[1] = pcks.ringContext.NewPoly()
	return
}

// GenShareRoundThree is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + u_i * pk[0] + e_0i, u_i * pk[1] + e_1i]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *ring.Poly, pk *bfv.PublicKey, ct *bfv.Ciphertext, shareOut PCKSShare) {

	//u_i
	_ = pcks.ternarySampler.SampleMontgomeryNTT(0.5, pcks.tmp)

	// h_0 = u_i * pk_0 (NTT)
	pcks.ringContext.MulCoeffsMontgomery(pcks.tmp, pk.Get()[0], shareOut.share[0])
	// h_1 = u_i * pk_1 (NTT)
	pcks.ringContext.MulCoeffsMontgomery(pcks.tmp, pk.Get()[1], shareOut.share[1])

	// h0 = u_i * pk_0 + s_i*c_1 (NTT)
	pcks.ringContext.NTT(ct.Value()[1], pcks.tmp)
	pcks.ringContext.MulCoeffsMontgomeryAndAdd(sk, pcks.tmp, shareOut.share[0])

	pcks.ringContext.InvNTT(shareOut.share[0], shareOut.share[0])
	pcks.ringContext.InvNTT(shareOut.share[1], shareOut.share[1])

	// h_0 = InvNTT(s_i*c_1 + u_i * pk_0) + e0
	pcks.gaussianSamplerSmudge.Sample(pcks.tmp)
	pcks.ringContext.Add(shareOut.share[0], pcks.tmp, shareOut.share[0])

	// h_1 = InvNTT(u_i * pk_1) + e1
	pcks.gaussianSampler.Sample(pcks.tmp)
	pcks.ringContext.Add(shareOut.share[1], pcks.tmp, shareOut.share[1])

	pcks.tmp.Zero()
}

// GenShareRoundTwo is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {
	pcks.ringContext.Add(share1.share[0], share2.share[0], shareOut.share[0])
	pcks.ringContext.Add(share1.share[1], share2.share[1], shareOut.share[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *bfv.Ciphertext) {

	pcks.ringContext.Add(ct.Value()[0], combined.share[0], ctOut.Value()[0])
	ctOut.Value()[1].Copy(combined.share[1])
}
