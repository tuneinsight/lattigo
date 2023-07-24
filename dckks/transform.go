package dckks

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s EncToShareProtocol
	s2e ShareToEncProtocol

	noise ring.DistributionParameters

	defaultScale *big.Int
	prec         uint

	tmpMask []*big.Int
	encoder *ckks.Encoder
}

// ShallowCopy creates a shallow copy of MaskedTransformProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// MaskedTransformProtocol can be used concurrently.
func (rfp MaskedTransformProtocol) ShallowCopy() MaskedTransformProtocol {

	params := rfp.e2s.params

	tmpMask := make([]*big.Int, params.N())
	for i := range rfp.tmpMask {
		tmpMask[i] = new(big.Int)
	}

	return MaskedTransformProtocol{
		e2s:          rfp.e2s.ShallowCopy(),
		s2e:          rfp.s2e.ShallowCopy(),
		prec:         rfp.prec,
		defaultScale: rfp.defaultScale,
		tmpMask:      tmpMask,
		encoder:      rfp.encoder.ShallowCopy(),
	}
}

// WithParams creates a shallow copy of the target MaskedTransformProtocol but with new output parameters.
// The expected input parameters remain unchanged.
func (rfp MaskedTransformProtocol) WithParams(paramsOut ckks.Parameters) MaskedTransformProtocol {

	tmpMask := make([]*big.Int, rfp.e2s.params.N())
	for i := range rfp.tmpMask {
		tmpMask[i] = new(big.Int)
	}

	s2e, err := NewShareToEncProtocol(paramsOut, rfp.noise)

	if err != nil {
		panic(err)
	}

	return MaskedTransformProtocol{
		e2s:          rfp.e2s.ShallowCopy(),
		s2e:          s2e,
		prec:         rfp.prec,
		defaultScale: rfp.defaultScale,
		tmpMask:      tmpMask,
		encoder:      rfp.encoder.ShallowCopy(),
	}
}

// MaskedTransformFunc represents a user-defined in-place function that can be evaluated on masked CKKS plaintexts, as a part of the
// Masked Transform Protocol.
// The function is called with a vector of *Complex modulo ckks.Parameters.Slots() as input, and must write
// its output on the same buffer.
// Transform can be the identity.
// Decode: if true, then the masked CKKS plaintext will be decoded before applying Transform.
// Recode: if true, then the masked CKKS plaintext will be recoded after applying Transform.
// i.e. : Decode (true/false) -> Transform -> Recode (true/false).
type MaskedTransformFunc struct {
	Decode bool
	Func   func(coeffs []*bignum.Complex)
	Encode bool
}

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
// paramsIn: the ckks.Parameters of the ciphertext before the protocol.
// paramsOut: the ckks.Parameters of the ciphertext after the protocol.
// prec : the log2 of decimal precision of the internal encoder.
// The method will return an error if the maximum number of slots of the output parameters is smaller than the number of slots of the input ciphertext.
func NewMaskedTransformProtocol(paramsIn, paramsOut ckks.Parameters, prec uint, noise ring.DistributionParameters) (rfp MaskedTransformProtocol, err error) {

	rfp = MaskedTransformProtocol{}

	rfp.noise = noise

	if rfp.e2s, err = NewEncToShareProtocol(paramsIn, noise); err != nil {
		return
	}

	if rfp.s2e, err = NewShareToEncProtocol(paramsOut, noise); err != nil {
		return
	}

	rfp.prec = prec

	scale := paramsOut.DefaultScale().Value

	rfp.defaultScale, _ = new(big.Float).SetPrec(prec).Set(&scale).Int(nil)

	rfp.tmpMask = make([]*big.Int, paramsIn.N())
	for i := range rfp.tmpMask {
		rfp.tmpMask[i] = new(big.Int)
	}

	rfp.encoder = ckks.NewEncoder(paramsIn, prec)

	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp MaskedTransformProtocol) AllocateShare(levelDecrypt, levelRecrypt int) drlwe.RefreshShare {
	return drlwe.RefreshShare{EncToShareShare: rfp.e2s.AllocateShare(levelDecrypt), ShareToEncShare: rfp.s2e.AllocateShare(levelRecrypt)}
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string. The CRP is considered to be in the NTT domain.
func (rfp MaskedTransformProtocol) SampleCRP(level int, crs sampling.PRNG) drlwe.KeySwitchCRP {
	return rfp.s2e.SampleCRP(level, crs)
}

// GenShare generates the shares of the PermuteProtocol
// This protocol requires additional inputs which are :
// skIn     : the secret-key if the input ciphertext.
// skOut    : the secret-key of the output ciphertext.
// logBound : the bit length of the masks.
// ct1      : the degree 1 element the ciphertext to refresh, i.e. ct1 = ckk.Ciphetext.Value[1].
// scale    : the scale of the ciphertext when entering the refresh.
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which the masked transform can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (rfp MaskedTransformProtocol) GenShare(skIn, skOut *rlwe.SecretKey, logBound uint, ct *rlwe.Ciphertext, crs drlwe.KeySwitchCRP, transform *MaskedTransformFunc, shareOut *drlwe.RefreshShare) (err error) {

	ringQ := rfp.s2e.params.RingQ()

	ct1 := ct.Value[1]

	if ct1.Level() < shareOut.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot GenShare: ct[1] level must be at least equal to EncToShareShare level")
	}

	if crs.Value.Level() != shareOut.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot GenShare: crs level must be equal to ShareToEncShare")
	}

	slots := ct.Slots()

	dslots := slots
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	// Generates the decryption share
	// Returns [M_i] on rfp.tmpMask and [a*s_i -M_i + e] on EncToShareShare
	if err = rfp.e2s.GenShare(skIn, logBound, ct, &drlwe.AdditiveShareBigint{Value: rfp.tmpMask}, &shareOut.EncToShareShare); err != nil {
		return
	}

	// Applies LT(M_i)
	if transform != nil {

		bigComplex := make([]*bignum.Complex, slots)

		for i := range bigComplex {
			bigComplex[i] = bignum.NewComplex()
			bigComplex[i][0].SetPrec(rfp.prec)
			bigComplex[i][1].SetPrec(rfp.prec)
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
			return fmt.Errorf("cannot GenShare: invalid ring type")
		}

		// Decodes if asked to
		if transform.Decode {
			if err := rfp.encoder.FFT(bigComplex[:slots], ct.LogSlots()); err != nil {
				return err
			}
		}

		// Applies the linear transform
		transform.Func(bigComplex)

		// Recodes if asked to
		if transform.Encode {
			if err := rfp.encoder.IFFT(bigComplex[:slots], ct.LogSlots()); err != nil {
				return err
			}
		}

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
	inputScaleInt, _ := new(big.Float).SetPrec(256).Set(&ct.Scale.Value).Int(nil)

	// Scales the mask by the ratio between the two scales
	for i := 0; i < dslots; i++ {
		rfp.tmpMask[i].Mul(rfp.tmpMask[i], rfp.defaultScale)
		rfp.tmpMask[i].Quo(rfp.tmpMask[i], inputScaleInt)
	}

	// Returns [-a*s_i + LT(M_i) * diffscale + e] on ShareToEncShare
	return rfp.s2e.GenShare(skOut, crs, ct.MetaData, drlwe.AdditiveShareBigint{Value: rfp.tmpMask}, &shareOut.ShareToEncShare)
}

// AggregateShares sums share1 and share2 on shareOut.
func (rfp MaskedTransformProtocol) AggregateShares(share1, share2, shareOut *drlwe.RefreshShare) (err error) {

	if share1.EncToShareShare.Value.Level() != share2.EncToShareShare.Value.Level() || share1.EncToShareShare.Value.Level() != shareOut.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot AggregateShares: all e2s shares must be at the same level")
	}

	if share1.ShareToEncShare.Value.Level() != share2.ShareToEncShare.Value.Level() || share1.ShareToEncShare.Value.Level() != shareOut.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot AggregateShares: all s2e shares must be at the same level")
	}

	rfp.e2s.params.RingQ().AtLevel(share1.EncToShareShare.Value.Level()).Add(share1.EncToShareShare.Value, share2.EncToShareShare.Value, shareOut.EncToShareShare.Value)
	rfp.s2e.params.RingQ().AtLevel(share1.ShareToEncShare.Value.Level()).Add(share1.ShareToEncShare.Value, share2.ShareToEncShare.Value, shareOut.ShareToEncShare.Value)

	return
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
// The ciphertext scale is reset to the default scale.
func (rfp MaskedTransformProtocol) Transform(ct *rlwe.Ciphertext, transform *MaskedTransformFunc, crs drlwe.KeySwitchCRP, share drlwe.RefreshShare, ciphertextOut *rlwe.Ciphertext) (err error) {

	if ct.Level() < share.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot Transform: input ciphertext level must be at least equal to e2s level")
	}

	maxLevel := crs.Value.Level()

	if maxLevel != share.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot Transform: crs level and s2e level must be the same")
	}

	ringQ := rfp.s2e.params.RingQ().AtLevel(maxLevel)

	slots := ct.Slots()

	dslots := slots
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	// Returns -sum(M_i) + x (outside of the NTT domain)

	rfp.e2s.GetShare(nil, share.EncToShareShare, ct, &drlwe.AdditiveShareBigint{Value: rfp.tmpMask[:dslots]})

	// Returns LT(-sum(M_i) + x)
	if transform != nil {

		bigComplex := make([]*bignum.Complex, slots)

		for i := range bigComplex {
			bigComplex[i] = bignum.NewComplex()
			bigComplex[i][0].SetPrec(rfp.prec)
			bigComplex[i][1].SetPrec(rfp.prec)
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
			return fmt.Errorf("cannot Transform: invalid ring type")
		}

		// Decodes if asked to
		if transform.Decode {
			if err := rfp.encoder.FFT(bigComplex[:slots], ct.LogSlots()); err != nil {
				return err
			}
		}

		// Applies the linear transform
		transform.Func(bigComplex)

		// Recodes if asked to
		if transform.Encode {
			if err := rfp.encoder.IFFT(bigComplex[:slots], ct.LogSlots()); err != nil {
				return err
			}
		}

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

	scale := ct.Scale.Value

	// Returns LT(-sum(M_i) + x) * diffscale
	inputScaleInt, _ := new(big.Float).Set(&scale).Int(nil)

	// Scales the mask by the ratio between the two scales
	for i := 0; i < dslots; i++ {
		rfp.tmpMask[i].Mul(rfp.tmpMask[i], rfp.defaultScale)
		rfp.tmpMask[i].Quo(rfp.tmpMask[i], inputScaleInt)
	}

	// Extend the levels of the ciphertext for future allocation
	if ciphertextOut.Value[0].N() != ringQ.N() {
		for i := range ciphertextOut.Value {
			ciphertextOut.Value[i] = ringQ.NewPoly()
		}
	} else {
		ciphertextOut.Resize(ciphertextOut.Degree(), maxLevel)
	}

	// Sets LT(-sum(M_i) + x) * diffscale in the RNS domain
	ringQ.SetCoefficientsBigint(rfp.tmpMask[:dslots], ciphertextOut.Value[0])

	rlwe.NTTSparseAndMontgomery(ringQ, ct.MetaData, ciphertextOut.Value[0])

	// LT(-sum(M_i) + x) * diffscale + [-a*s + LT(M_i) * diffscale + e] = [-a*s + LT(x) * diffscale + e]
	ringQ.Add(ciphertextOut.Value[0], share.ShareToEncShare.Value, ciphertextOut.Value[0])

	// Copies the result on the out ciphertext
	if err = rfp.s2e.GetEncryption(drlwe.KeySwitchShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut); err != nil {
		return
	}

	*ciphertextOut.MetaData = *ct.MetaData
	ciphertextOut.Scale = rfp.s2e.params.DefaultScale()

	return
}
