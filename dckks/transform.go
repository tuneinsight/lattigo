package dckks

import (
	"fmt"
	"math/big"

	"encoding/binary"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	sigmaSmudging float64

	defaultScale *big.Int
	precision    int

	tmpMask []*big.Int
	encoder ckks.EncoderBigComplex
}

// ShallowCopy creates a shallow copy of MaskedTransformProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// MaskedTransformProtocol can be used concurrently.
func (rfp *MaskedTransformProtocol) ShallowCopy() *MaskedTransformProtocol {

	params := rfp.e2s.params
	precision := rfp.precision

	tmpMask := make([]*big.Int, params.N())
	for i := range rfp.tmpMask {
		tmpMask[i] = new(big.Int)
	}

	return &MaskedTransformProtocol{
		e2s:          *rfp.e2s.ShallowCopy(),
		s2e:          *rfp.s2e.ShallowCopy(),
		precision:    precision,
		defaultScale: rfp.defaultScale,
		tmpMask:      tmpMask,
		encoder:      rfp.encoder.ShallowCopy(),
	}
}

// WithParams creates a shallow copy of the target MaskedTransformProtocol but with new output parameters.
// The expected input parameters remain unchnaged.
func (rfp *MaskedTransformProtocol) WithParams(paramsOut ckks.Parameters) *MaskedTransformProtocol {

	precision := rfp.precision

	tmpMask := make([]*big.Int, rfp.e2s.params.N())
	for i := range rfp.tmpMask {
		tmpMask[i] = new(big.Int)
	}

	return &MaskedTransformProtocol{
		e2s:          *rfp.e2s.ShallowCopy(),
		s2e:          *NewS2EProtocol(paramsOut, rfp.sigmaSmudging),
		precision:    precision,
		defaultScale: rfp.defaultScale,
		tmpMask:      tmpMask,
		encoder:      rfp.encoder.ShallowCopy(),
	}
}

// MaskedTransformFunc represents a user-defined in-place function that can be evaluated on masked CKKS plaintexts, as a part of the
// Masked Transform Protocol.
// The function is called with a vector of *ring.Complex modulo ckks.Parameters.Slots() as input, and must write
// its output on the same buffer.
type MaskedTransformFunc func(coeffs []*ring.Complex)

// MaskedTransformShare is a struct storing the decryption and recryption shares.
type MaskedTransformShare struct {
	e2sShare drlwe.CKSShare
	s2eShare drlwe.CKSShare
}

// MarshalBinary encodes a RefreshShare on a slice of bytes.
func (share *MaskedTransformShare) MarshalBinary() (data []byte, err error) {
	var e2sData, s2eData []byte
	if e2sData, err = share.e2sShare.MarshalBinary(); err != nil {
		return nil, err
	}
	if s2eData, err = share.s2eShare.MarshalBinary(); err != nil {
		return nil, err
	}
	data = make([]byte, 8)
	binary.LittleEndian.PutUint64(data, uint64(len(e2sData)))
	data = append(data, e2sData...)
	data = append(data, s2eData...)
	return data, nil
}

// UnmarshalBinary decodes a marshaled RefreshShare on the target RefreshShare.
func (share *MaskedTransformShare) UnmarshalBinary(data []byte) error {

	e2sDataLen := binary.LittleEndian.Uint64(data[:8])

	if err := share.e2sShare.UnmarshalBinary(data[8 : e2sDataLen+8]); err != nil {
		return err
	}
	if err := share.s2eShare.UnmarshalBinary(data[8+e2sDataLen:]); err != nil {
		return err
	}
	return nil
}

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
// paramsIn: the ckks.Parameters of the ciphertext before the protocol.
// paramsOut: the ckks.Parameters of the ciphertext after the protocol.
// precision : the log2 of decimal precision of the internal encoder.
// The method will return an error if the maximum number of slots of the output parameters is smaller than the number of slots of the input ciphertext.
func NewMaskedTransformProtocol(paramsIn, paramsOut ckks.Parameters, precision int, sigmaSmudging float64) (rfp *MaskedTransformProtocol, err error) {

	if paramsIn.Slots() > paramsOut.MaxSlots() {
		return nil, fmt.Errorf("newMaskedTransformProtocol: paramsOut.N()/2 < paramsIn.Slots()")
	}

	rfp = new(MaskedTransformProtocol)

	rfp.sigmaSmudging = sigmaSmudging

	rfp.e2s = *NewE2SProtocol(paramsIn, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(paramsOut, sigmaSmudging)

	rfp.precision = precision

	rfp.defaultScale = new(big.Int)
	ring.NewFloat(paramsOut.DefaultScale(), precision).Int(rfp.defaultScale)

	rfp.tmpMask = make([]*big.Int, paramsIn.N())
	for i := range rfp.tmpMask {
		rfp.tmpMask[i] = new(big.Int)
	}
	rfp.encoder = ckks.NewEncoderBigComplex(paramsIn, precision)
	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShare(levelDecrypt, levelRecrypt int) *MaskedTransformShare {
	return &MaskedTransformShare{*rfp.e2s.AllocateShare(levelDecrypt), *rfp.s2e.AllocateShare(levelRecrypt)}
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string. The CRP is considered to be in the NTT domain.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs utils.PRNG) drlwe.CKSCRP {
	crp := rfp.s2e.SampleCRP(level, crs)
	crp.IsNTT = true
	return crp
}

// GenShare generates the shares of the PermuteProtocol
// This protocol requires additional inputs which are :
// skIn     : the secret-key if the input ciphertext.
// skOut    : the secret-key of the output ciphertext.
// logBound : the bit length of the masks.
// logSlots : the bit length of the number of slots.
// ct1      : the degree 1 element the ciphertext to refresh, i.e. ct1 = ckk.Ciphetext.Value[1].
// scale    : the scale of the ciphertext when entering the refresh.
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which the masked transform can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (rfp *MaskedTransformProtocol) GenShare(skIn, skOut *rlwe.SecretKey, logBound, logSlots int, ct1 *ring.Poly, scale float64, crs drlwe.CKSCRP, transform MaskedTransformFunc, shareOut *MaskedTransformShare) {

	ringQ := rfp.s2e.params.RingQ()

	if ct1.Level() < shareOut.e2sShare.Value.Level() {
		panic("ct[1] level must be at least equal to e2sShare level")
	}

	if (*ring.Poly)(&crs).Level() != shareOut.s2eShare.Value.Level() {
		panic("crs level must be equal to s2eShare")
	}

	slots := 1 << logSlots

	dslots := slots
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	// Generates the decryption share
	// Returns [M_i] on rfp.tmpMask and [a*s_i -M_i + e] on e2sShare
	rfp.e2s.GenShare(skIn, logBound, logSlots, ct1, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask}, &shareOut.e2sShare)

	// Applies LT(M_i)
	if transform != nil {

		bigComplex := make([]*ring.Complex, slots)

		for i := range bigComplex {
			bigComplex[i] = ring.NewComplex(ring.NewFloat(0, rfp.precision), ring.NewFloat(0, rfp.precision))
		}

		// Extracts sparse coefficients
		for i := 0; i < slots; i++ {
			bigComplex[i][0].SetInt(rfp.tmpMask[i])
		}

		switch rfp.e2s.params.RingType() {
		case ring.Standard:
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				bigComplex[i][1].SetInt(rfp.tmpMask[j])
			}
		case ring.ConjugateInvariant:
			for i := 1; i < slots; i++ {
				bigComplex[i][1].Neg(bigComplex[slots-i][0])
			}
		default:
			panic("invalid ring type")
		}

		// Decodes
		rfp.encoder.FFT(bigComplex, 1<<logSlots)

		// Applies the linear transform
		transform(bigComplex)

		// Recodes
		rfp.encoder.InvFFT(bigComplex, 1<<logSlots)

		// Puts the coefficient back
		for i := 0; i < slots; i++ {
			bigComplex[i].Real().Int(rfp.tmpMask[i])
		}

		if rfp.e2s.params.RingType() == ring.Standard {
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				bigComplex[i].Imag().Int(rfp.tmpMask[j])
			}
		}
	}

	// Applies LT(M_i) * diffscale
	inputScaleInt := new(big.Int)
	ring.NewFloat(scale, 256).Int(inputScaleInt)

	// Scales the mask by the ratio between the two scales
	for i := 0; i < dslots; i++ {
		rfp.tmpMask[i].Mul(rfp.tmpMask[i], rfp.defaultScale)
		rfp.tmpMask[i].Quo(rfp.tmpMask[i], inputScaleInt)
	}

	// Returns [-a*s_i + LT(M_i) * diffscale + e] on s2eShare
	rfp.s2e.GenShare(skOut, crs, logSlots, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask}, &shareOut.s2eShare)
}

// AggregateShare sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) AggregateShare(share1, share2, shareOut *MaskedTransformShare) {

	if share1.e2sShare.Value.Level() != share2.e2sShare.Value.Level() || share1.e2sShare.Value.Level() != shareOut.e2sShare.Value.Level() {
		panic("all e2s shares must be at the same level")
	}

	if share1.s2eShare.Value.Level() != share2.s2eShare.Value.Level() || share1.s2eShare.Value.Level() != shareOut.s2eShare.Value.Level() {
		panic("all s2e shares must be at the same level")
	}

	rfp.e2s.params.RingQ().AddLvl(share1.e2sShare.Value.Level(), share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.s2e.params.RingQ().AddLvl(share1.s2eShare.Value.Level(), share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
// The ciphertext scale is reset to the default scale.
func (rfp *MaskedTransformProtocol) Transform(ct *ckks.Ciphertext, logSlots int, transform MaskedTransformFunc, crs drlwe.CKSCRP, share *MaskedTransformShare, ciphertextOut *ckks.Ciphertext) {

	if ct.Level() < share.e2sShare.Value.Level() {
		panic("input ciphertext level must be at least equal to e2s level")
	}

	maxLevel := (*ring.Poly)(&crs).Level()

	if maxLevel != share.s2eShare.Value.Level() {
		panic("crs level and s2e level must be the same")
	}

	ringQ := rfp.s2e.params.RingQ()

	slots := 1 << logSlots

	dslots := slots
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	// Returns -sum(M_i) + x (outside of the NTT domain)

	rfp.e2s.GetShare(nil, &share.e2sShare, logSlots, ct, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask[:dslots]})

	// Returns LT(-sum(M_i) + x)
	if transform != nil {

		bigComplex := make([]*ring.Complex, slots)

		for i := range bigComplex {
			bigComplex[i] = ring.NewComplex(ring.NewFloat(0, rfp.precision), ring.NewFloat(0, rfp.precision))
		}

		// Extracts sparse coefficients
		for i := 0; i < slots; i++ {
			bigComplex[i][0].SetInt(rfp.tmpMask[i])
		}

		switch rfp.e2s.params.RingType() {
		case ring.Standard:
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				bigComplex[i][1].SetInt(rfp.tmpMask[j])
			}
		case ring.ConjugateInvariant:
			for i := 1; i < slots; i++ {
				bigComplex[i][1].Neg(bigComplex[slots-i][0])
			}
		default:
			panic("invalid ring type")
		}

		// Decodes
		rfp.encoder.FFT(bigComplex, 1<<logSlots)

		// Applies the linear transform
		transform(bigComplex)

		// Recodes
		rfp.encoder.InvFFT(bigComplex, 1<<logSlots)

		// Puts the coefficient back
		for i := 0; i < slots; i++ {
			bigComplex[i].Real().Int(rfp.tmpMask[i])
		}

		if rfp.e2s.params.RingType() == ring.Standard {
			for i := 0; i < slots; i++ {
				bigComplex[i].Imag().Int(rfp.tmpMask[i+slots])
			}
		}
	}

	// Returns LT(-sum(M_i) + x) * diffscale
	inputScaleInt := new(big.Int)
	ring.NewFloat(ct.Scale, 256).Int(inputScaleInt)

	// Scales the mask by the ratio between the two scales
	for i := 0; i < dslots; i++ {
		rfp.tmpMask[i].Mul(rfp.tmpMask[i], rfp.defaultScale)
		rfp.tmpMask[i].Quo(rfp.tmpMask[i], inputScaleInt)
	}

	// Extend the levels of the ciphertext for future allocation
	if ciphertextOut.Value[0].N() != ringQ.N {
		for i := range ciphertextOut.Value {
			isNTT := ciphertextOut.Value[i].IsNTT
			ciphertextOut.Value[i] = ringQ.NewPolyLvl(maxLevel)
			ciphertextOut.Value[i].IsNTT = isNTT
		}
	} else {
		ciphertextOut.Resize(ciphertextOut.Degree(), maxLevel)
	}

	// Sets LT(-sum(M_i) + x) * diffscale in the RNS domain
	ringQ.SetCoefficientsBigintLvl(maxLevel, rfp.tmpMask[:dslots], ciphertextOut.Value[0])

	ckks.NttAndMontgomeryLvl(maxLevel, logSlots, ringQ, false, ciphertextOut.Value[0])

	// LT(-sum(M_i) + x) * diffscale + [-a*s + LT(M_i) * diffscale + e] = [-a*s + LT(x) * diffscale + e]
	ringQ.AddLvl(maxLevel, ciphertextOut.Value[0], share.s2eShare.Value, ciphertextOut.Value[0])

	// Copies the result on the out ciphertext
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)

	ciphertextOut.Scale = rfp.s2e.params.DefaultScale()
}
