package drlwe

import (
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// CKSProtocol is the structure storing the parameters and and precomputations for the collective key-switching protocol.
type CKSProtocol struct {
	params          rlwe.Parameters
	sigmaSmudging   float64
	gaussianSampler *ring.GaussianSampler
	basisExtender   *ring.BasisExtender
	buff            *ring.Poly
	buffDelta       *ring.Poly
}

// ShallowCopy creates a shallow copy of CKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// CKSProtocol can be used concurrently.
func (cks *CKSProtocol) ShallowCopy() *CKSProtocol {
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := cks.params

	return &CKSProtocol{
		params:          params,
		gaussianSampler: ring.NewGaussianSampler(prng, params.RingQ(), cks.sigmaSmudging, int(6*cks.sigmaSmudging)),
		basisExtender:   cks.basisExtender.ShallowCopy(),
		buff:            params.RingQ().NewPoly(),
		buffDelta:       params.RingQ().NewPoly(),
	}
}

// CKSCRP is a type for common reference polynomials in the CKS protocol.
type CKSCRP ring.Poly

// NewCKSProtocol creates a new CKSProtocol that will be used to perform a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) *CKSProtocol {
	cks := new(CKSProtocol)
	cks.params = params
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	// EncFreshSK + sigmaSmudging
	cks.sigmaSmudging = math.Sqrt(params.Sigma()*params.Sigma() + sigmaSmudging*sigmaSmudging)

	cks.gaussianSampler = ring.NewGaussianSampler(prng, params.RingQ(), cks.sigmaSmudging, int(6*cks.sigmaSmudging))

	if cks.params.RingP() != nil {
		cks.basisExtender = ring.NewBasisExtender(params.RingQ(), params.RingP())
	}
	cks.buff = params.RingQ().NewPoly()
	cks.buffDelta = params.RingQ().NewPoly()
	return cks
}

// AllocateShare allocates the shares of the CKSProtocol
func (cks *CKSProtocol) AllocateShare(level int) *CKSShare {
	return &CKSShare{cks.params.RingQ().AtLevel(level).NewPoly()}
}

// SampleCRP samples a common random polynomial to be used in the CKS protocol from the provided
// common reference string.
func (cks *CKSProtocol) SampleCRP(level int, crs CRS) CKSCRP {
	ringQ := cks.params.RingQ().AtLevel(level)
	crp := ringQ.NewPoly()
	ring.NewUniformSampler(crs, ringQ).Read(crp)
	return CKSCRP(*crp)
}

// GenShare computes a party's share in the CKS protocol from secret-key skInput to secret-key skOutput.
// ct is the rlwe.Ciphertext to keyswitch. Note that ct.Value[0] is not used by the function and can be nil/zero.
//
// Expected noise: ctNoise + encFreshSk + smudging
func (cks *CKSProtocol) GenShare(skInput, skOutput *rlwe.SecretKey, ct *rlwe.Ciphertext, shareOut *CKSShare) {

	levelQ := utils.MinInt(shareOut.Value.Level(), ct.Value[1].Level())

	shareOut.Value.Resize(levelQ)

	ringQ := cks.params.RingQ().AtLevel(levelQ)

	ringQ.Sub(skInput.Value.Q, skOutput.Value.Q, cks.buffDelta)

	var c1NTT *ring.Poly
	if !ct.IsNTT {
		ringQ.NTTLazy(ct.Value[1], cks.buff)
		c1NTT = cks.buff
	} else {
		c1NTT = ct.Value[1]
	}

	// c1NTT * (skIn - skOut)
	ringQ.MulCoeffsMontgomeryLazy(c1NTT, cks.buffDelta, shareOut.Value)

	if !ct.IsNTT {
		// InvNTT(c1NTT * (skIn - skOut)) + e
		ringQ.INTTLazy(shareOut.Value, shareOut.Value)
		cks.gaussianSampler.AtLevel(levelQ).ReadAndAdd(shareOut.Value)
	} else {
		// c1NTT * (skIn - skOut) + e
		cks.gaussianSampler.AtLevel(levelQ).Read(cks.buff)
		ringQ.NTT(cks.buff, cks.buff)
		ringQ.Add(shareOut.Value, cks.buff, shareOut.Value)
	}
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut *CKSShare) {
	if share1.Level() != share2.Level() || share1.Level() != shareOut.Level() {
		panic("shares levels do not match")
	}

	cks.params.RingQ().AtLevel(share1.Level()).Add(share1.Value, share2.Value, shareOut.Value)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(ctIn *rlwe.Ciphertext, combined *CKSShare, ctOut *rlwe.Ciphertext) {

	level := ctIn.Level()

	if ctIn != ctOut {

		ctOut.Resize(ctIn.Degree(), level)

		ring.CopyLvl(level, ctIn.Value[1], ctOut.Value[1])

		ctOut.MetaData = ctIn.MetaData
	}

	cks.params.RingQ().AtLevel(level).Add(ctIn.Value[0], combined.Value, ctOut.Value[0])
}

// CKSShare is a type for the CKS protocol shares.
type CKSShare struct {
	Value *ring.Poly
}

// Level returns the level of the target share.
func (ckss *CKSShare) Level() int {
	return ckss.Value.Level()
}

// BinarySize returns the size in bytes that the object once marshalled into a binary form.
func (ckss *CKSShare) BinarySize() int {
	return ckss.Value.BinarySize()
}

// MarshalBinary encodes a CKS share on a slice of bytes.
func (ckss *CKSShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, ckss.BinarySize())
	_, err = ckss.MarshalBinaryInPlace(data)
	return
}

// MarshalBinaryInPlace encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (ckss *CKSShare) MarshalBinaryInPlace(data []byte) (ptr int, err error) {
	return ckss.Value.Read(data)
}

// UnmarshalBinary decodes a slice of bytes generated by MarshalBinary
// or MarshalBinaryInPlace on the object.
func (ckss *CKSShare) UnmarshalBinary(data []byte) (err error) {
	_, err = ckss.UnmarshalBinaryInPlace(data)
	return
}

// UnmarshalBinaryInPlace decodes a slice of bytes generated by MarshalBinary or
// MarshalBinaryInPlace on the object and returns the number of bytes read.
func (ckss *CKSShare) UnmarshalBinaryInPlace(data []byte) (ptr int, err error) {
	if ckss.Value == nil {
		ckss.Value = new(ring.Poly)
	}

	return ckss.Value.Write(data)
}
