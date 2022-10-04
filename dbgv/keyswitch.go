package dbgv

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v3/bgv"
	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	drlwe.CKSProtocol
	params   bgv.Parameters
	tInvModQ []*big.Int
	buffQ    *ring.Poly
}

// NewCKSProtocol creates a new CKSProtocol that will be used to perform a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params bgv.Parameters, sigmaSmudging float64) (cks *CKSProtocol) {

	ringQ := params.RingQ()
	T := params.T()
	tInvModQ := make([]*big.Int, len(ringQ.Modulus))
	for i := range ringQ.Modulus {
		tInvModQ[i] = ring.NewUint(T)
		tInvModQ[i].ModInverse(tInvModQ[i], ringQ.ModulusAtLevel[i])
	}

	buffQ := params.RingQ().NewPoly()
	buffQ.IsNTT = true

	return &CKSProtocol{
		CKSProtocol: *drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging),
		params:      params,
		tInvModQ:    tInvModQ,
		buffQ:       buffQ,
	}
}

// GenShare computes a party's share in the CKS protocol.
// ct1 is the degree 1 element of the rlwe.Ciphertext to keyswitch, i.e. ct1 = rlwe.Ciphertext.Value[1].
func (cks *CKSProtocol) GenShare(skInput, skOutput *rlwe.SecretKey, ct1 *ring.Poly, shareOut *drlwe.CKSShare) {
	level := utils.MinInt(ct1.Level(), shareOut.Value.Level())
	cks.params.RingQ().MulScalarBigintLvl(level, ct1, cks.tInvModQ[level], cks.buffQ)
	cks.CKSProtocol.GenShare(skInput, skOutput, cks.buffQ, shareOut)
	cks.params.RingQ().MulScalarLvl(level, shareOut.Value, cks.params.T(), shareOut.Value)
}

// KeySwitch performs the actual keyswitching operation on a Ciphertext ct and stores the result in ctOut
func (cks *CKSProtocol) KeySwitch(ctIn *bgv.Ciphertext, combined *drlwe.CKSShare, ctOut *bgv.Ciphertext) {
	cks.CKSProtocol.KeySwitch(ctIn.Ciphertext, combined, ctOut.Ciphertext)
	ctOut.SetScale(ctIn.Scale())
}

// ShallowCopy creates a shallow copy of CKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// CKSProtocol can be used concurrently.
func (cks *CKSProtocol) ShallowCopy() *CKSProtocol {
	buffQ := cks.params.RingQ().NewPoly()
	buffQ.IsNTT = true
	return &CKSProtocol{
		CKSProtocol: *cks.CKSProtocol.ShallowCopy(),
		params:      cks.params,
		tInvModQ:    cks.tInvModQ,
		buffQ:       buffQ,
	}
}

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	drlwe.PCKSProtocol
	params   bgv.Parameters
	buffQ    *ring.Poly
	tInvModQ []*big.Int
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a Ciphertext ctx encrypted under a key secret-shared among j parties under a new
// collective public-key.
func NewPCKSProtocol(params bgv.Parameters, sigmaSmudging float64) *PCKSProtocol {
	ringQ := params.RingQ()
	T := params.T()
	tInvModQ := make([]*big.Int, len(ringQ.Modulus))
	for i := range ringQ.Modulus {
		tInvModQ[i] = ring.NewUint(T)
		tInvModQ[i].ModInverse(tInvModQ[i], ringQ.ModulusAtLevel[i])
	}

	buffQ := params.RingQ().NewPoly()
	buffQ.IsNTT = true

	return &PCKSProtocol{
		params:       params,
		PCKSProtocol: *drlwe.NewPCKSProtocol(params.Parameters, sigmaSmudging),
		buffQ:        buffQ,
		tInvModQ:     tInvModQ,
	}
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ct[1] + (u_i * pk[0] + e_0i)/P, (u_i * pk[1] + e_1i)/P]
//
// and broadcasts the result to the other j-1 parties.
// ct1 is the degree 1 element of the rlwe.Ciphertext to keyswitch, i.e. ct1 = rlwe.Ciphertext.Value[1].
func (pcks *PCKSProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct1 *ring.Poly, shareOut *drlwe.PCKSShare) {
	level := utils.MinInt(ct1.Level(), shareOut.Value[0].Level())
	pcks.params.RingQ().MulScalarBigintLvl(level, ct1, pcks.tInvModQ[level], pcks.buffQ)
	pcks.PCKSProtocol.GenShare(sk, pk, pcks.buffQ, shareOut)
	pcks.params.RingQ().MulScalarLvl(level, shareOut.Value[0], pcks.params.T(), shareOut.Value[0])
	pcks.params.RingQ().MulScalarLvl(level, shareOut.Value[1], pcks.params.T(), shareOut.Value[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut.
func (pcks *PCKSProtocol) KeySwitch(ctIn *bgv.Ciphertext, combined *drlwe.PCKSShare, ctOut *bgv.Ciphertext) {
	pcks.PCKSProtocol.KeySwitch(ctIn.Ciphertext, combined, ctOut.Ciphertext)
	ctOut.SetScale(ctIn.Scale())
}

// ShallowCopy creates a shallow copy of PCKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// PCKSProtocol can be used concurrently.
func (pcks *PCKSProtocol) ShallowCopy() *PCKSProtocol {

	buffQ := pcks.params.RingQ().NewPoly()
	buffQ.IsNTT = true

	return &PCKSProtocol{
		PCKSProtocol: *pcks.PCKSProtocol.ShallowCopy(),
		params:       pcks.params,
		buffQ:        buffQ,
		tInvModQ:     pcks.tInvModQ,
	}
}
