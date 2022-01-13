package dckks

import (
	"math/big"

	"encoding/binary"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	defaultScale *big.Int
	precision    int

	tmpMask       []*big.Int
	tmpBigComplex []*ring.Complex
	encoder       ckks.EncoderBigComplex
}

// ShallowCopy creates a shallow copy of MaskedTransformProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// MaskedTransformProtocol can be used concurrently.
func (rfp *MaskedTransformProtocol) ShallowCopy() *MaskedTransformProtocol {

	params := rfp.e2s.params
	precision := rfp.precision

	tmpMask := make([]*big.Int, params.N())
	tmpBigComplex := make([]*ring.Complex, params.MaxSlots())
	for i := range rfp.tmpMask {
		tmpMask[i] = new(big.Int)
	}
	for i := range rfp.tmpBigComplex {
		tmpBigComplex[i] = ring.NewComplex(ring.NewFloat(0, precision), ring.NewFloat(0, precision))
	}

	return &MaskedTransformProtocol{
		e2s:           *rfp.e2s.ShallowCopy(),
		s2e:           *rfp.s2e.ShallowCopy(),
		precision:     precision,
		defaultScale:  rfp.defaultScale,
		tmpMask:       tmpMask,
		tmpBigComplex: tmpBigComplex,
		encoder:       rfp.encoder.ShallowCopy(),
	}
}

// MaskedTransformFunc is a method template for linear transforms that can be
// evaluated on a ciphertext during its collective refresh
type MaskedTransformFunc func(coeffsIn, coeffsOut []*ring.Complex)

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
func NewMaskedTransformProtocol(params ckks.Parameters, precision int, sigmaSmudging float64) (rfp *MaskedTransformProtocol) {

	rfp = new(MaskedTransformProtocol)
	rfp.e2s = *NewE2SProtocol(params, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(params, sigmaSmudging)
	rfp.precision = precision

	rfp.defaultScale = new(big.Int)
	ring.NewFloat(params.DefaultScale(), precision).Int(rfp.defaultScale)

	rfp.tmpMask = make([]*big.Int, params.N())
	rfp.tmpBigComplex = make([]*ring.Complex, params.MaxSlots())
	for i := range rfp.tmpMask {
		rfp.tmpMask[i] = new(big.Int)
	}
	for i := range rfp.tmpBigComplex {
		rfp.tmpBigComplex[i] = ring.NewComplex(ring.NewFloat(0, precision), ring.NewFloat(0, precision))
	}
	rfp.encoder = ckks.NewEncoderBigComplex(params, precision)
	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShare(levelDecrypt, levelRecrypt int) *MaskedTransformShare {
	return &MaskedTransformShare{*rfp.e2s.AllocateShare(levelDecrypt), *rfp.s2e.AllocateShare(levelRecrypt)}
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs utils.PRNG) drlwe.CKSCRP {
	return rfp.s2e.SampleCRP(level, crs)
}

// GenShares generates the shares of the PermuteProtocol
// This protocol requires additional inputs which are :
// logBound : the bit length of the masks
// logSlots : the bit length of the number of slots
//
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which the masked transform can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (rfp *MaskedTransformProtocol) GenShares(sk *rlwe.SecretKey, logBound, logSlots int, c1 *ring.Poly, scale float64, crs drlwe.CKSCRP, transform MaskedTransformFunc, shareOut *MaskedTransformShare) {

	ringQ := rfp.s2e.params.RingQ()

	if c1.Level() < shareOut.e2sShare.Value.Level() {
		panic("ct[1] level must be at least equal to e2sShare level")
	}

	if (*ring.Poly)(&crs).Level() != shareOut.s2eShare.Value.Level() {
		panic("crs level must be equal to s2eShare")
	}

	// Generates the decryption share
	// Returns [M_i] on rfp.tmpMask and [a*s_i -M_i + e] on e2sShare
	rfp.e2s.GenShare(sk, logBound, logSlots, c1, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask}, &shareOut.e2sShare)

	slots := 1 << logSlots
	gap := rfp.e2s.params.MaxSlots() / slots

	// Applies LT(M_i)
	if transform != nil {

		// Extracts sparse coefficients
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			rfp.tmpBigComplex[i][0].SetInt(rfp.tmpMask[idx])
		}

		switch rfp.e2s.params.RingType() {
		case ring.Standard:
			for i, idx := 0, ringQ.N>>1; i < slots; i, idx = i+1, idx+gap {
				rfp.tmpBigComplex[i][1].SetInt(rfp.tmpMask[idx])
			}
		case ring.ConjugateInvariant:
			for i := 1; i < slots; i++ {
				rfp.tmpBigComplex[i][1].Neg(rfp.tmpBigComplex[slots-i][0])
			}
		default:
			panic("invalid ring type")
		}

		// Decodes
		rfp.encoder.FFT(rfp.tmpBigComplex, 1<<logSlots)

		// Applies the linear transform
		transform(rfp.tmpBigComplex, rfp.tmpBigComplex)

		// Recodes
		rfp.encoder.InvFFT(rfp.tmpBigComplex, 1<<logSlots)

		// Puts the coefficient back
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			rfp.tmpBigComplex[i].Real().Int(rfp.tmpMask[idx])
		}

		if rfp.e2s.params.RingType() == ring.Standard {
			for i, jdx := 0, ringQ.N>>1; i < slots; i, jdx = i+1, jdx+gap {
				rfp.tmpBigComplex[i].Imag().Int(rfp.tmpMask[jdx])
			}
		}
	}

	// Applies LT(M_i) * diffscale
	inputScaleInt := new(big.Int)
	ring.NewFloat(scale, 256).Int(inputScaleInt)

	// Scales the mask by the ratio between the two scales
	for i := 0; i < ringQ.N; i++ {
		if i%gap == 0 {
			rfp.tmpMask[i].Mul(rfp.tmpMask[i], rfp.defaultScale)
			rfp.tmpMask[i].Quo(rfp.tmpMask[i], inputScaleInt)
		}
	}

	// Returns [-a*s_i + LT(M_i) * diffscale + e] on s2eShare
	rfp.s2e.GenShare(sk, crs, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask}, &shareOut.s2eShare)
}

// Aggregate sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) Aggregate(share1, share2, shareOut *MaskedTransformShare) {

	if share1.e2sShare.Value.Level() != share2.e2sShare.Value.Level() || share1.e2sShare.Value.Level() != shareOut.e2sShare.Value.Level() {
		panic("all e2s shares must be at the same level")
	}

	if share1.s2eShare.Value.Level() != share2.s2eShare.Value.Level() || share1.s2eShare.Value.Level() != shareOut.s2eShare.Value.Level() {
		panic("all s2e shares must be at the same level")
	}

	ringQ := rfp.s2e.params.RingQ()

	ringQ.AddLvl(share1.e2sShare.Value.Level(), share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	ringQ.AddLvl(share1.s2eShare.Value.Level(), share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ct *ckks.Ciphertext, logSlots int, transform MaskedTransformFunc, crs drlwe.CKSCRP, share *MaskedTransformShare, ciphertextOut *ckks.Ciphertext) {

	if ct.Level() < share.e2sShare.Value.Level() {
		panic("input ciphertext level must be at least equal to e2s level")
	}

	maxLevel := (*ring.Poly)(&crs).Level()

	if maxLevel != share.s2eShare.Value.Level() {
		panic("crs level and s2e level must be the same")
	}

	ringQ := rfp.s2e.params.RingQ()

	// Returns -sum(M_i) + x (outside of the NTT domain)
	rfp.e2s.GetShare(nil, &share.e2sShare, ct, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask})

	slots := 1 << logSlots
	gap := rfp.e2s.params.MaxSlots() / slots

	// Returns LT(-sum(M_i) + x)
	if transform != nil {
		// Extracts sparse coefficients
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			rfp.tmpBigComplex[i][0].SetInt(rfp.tmpMask[idx])
		}

		switch rfp.e2s.params.RingType() {
		case ring.Standard:
			for i, idx := 0, ringQ.N>>1; i < slots; i, idx = i+1, idx+gap {
				rfp.tmpBigComplex[i][1].SetInt(rfp.tmpMask[idx])
			}
		case ring.ConjugateInvariant:
			for i := 1; i < slots; i++ {
				rfp.tmpBigComplex[i][1].Neg(rfp.tmpBigComplex[slots-i][0])
			}
		default:
			panic("invalid ring type")
		}

		// Decodes
		rfp.encoder.FFT(rfp.tmpBigComplex, 1<<logSlots)

		// Applies the linear transform
		transform(rfp.tmpBigComplex, rfp.tmpBigComplex)

		// Recodes
		rfp.encoder.InvFFT(rfp.tmpBigComplex, 1<<logSlots)

		// Puts the coefficient back
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			rfp.tmpBigComplex[i].Real().Int(rfp.tmpMask[idx])
		}

		if rfp.e2s.params.RingType() == ring.Standard {
			for i, jdx := 0, ringQ.N>>1; i < slots; i, jdx = i+1, jdx+gap {
				rfp.tmpBigComplex[i].Imag().Int(rfp.tmpMask[jdx])
			}
		}
	}

	// Returns LT(-sum(M_i) + x) * diffscale
	inputScaleInt := new(big.Int)
	ring.NewFloat(ct.Scale, 256).Int(inputScaleInt)

	// Scales the mask by the ratio between the two scales
	for i := range rfp.tmpMask {
		rfp.tmpMask[i].Mul(rfp.tmpMask[i], rfp.defaultScale)
		rfp.tmpMask[i].Quo(rfp.tmpMask[i], inputScaleInt)
	}

	// Extend the levels of the ciphertext for future allocation
	for ciphertextOut.Level() != maxLevel {
		level := ciphertextOut.Level() + 1

		ciphertextOut.Value[0].Coeffs = append(ciphertextOut.Value[0].Coeffs, make([][]uint64, 1)...)
		ciphertextOut.Value[0].Coeffs[level] = make([]uint64, ringQ.N)

		ciphertextOut.Value[1].Coeffs = append(ciphertextOut.Value[1].Coeffs, make([][]uint64, 1)...)
		ciphertextOut.Value[1].Coeffs[level] = make([]uint64, ringQ.N)
	}

	// Sets LT(-sum(M_i) + x) * diffscale in the RNS domain
	ringQ.SetCoefficientsBigintLvl(maxLevel, rfp.tmpMask, ciphertextOut.Value[0])

	// Sets LT(-sum(M_i) + x) * diffscale in the NTT domain
	ringQ.NTTLvl(maxLevel, ciphertextOut.Value[0], ciphertextOut.Value[0])

	// LT(-sum(M_i) + x) * diffscale + [-a*s + LT(M_i) * diffscale + e] = [-a*s + LT(x) * diffscale + e]
	ringQ.AddLvl(maxLevel, ciphertextOut.Value[0], share.s2eShare.Value, ciphertextOut.Value[0])

	// Copies the result on the out ciphertext
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
