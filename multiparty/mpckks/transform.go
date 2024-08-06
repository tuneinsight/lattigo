package mpckks

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// MaskedLinearTransformationProtocol is a struct storing the parameters for the [MaskedLinearTransformationProtocol] protocol.
type MaskedLinearTransformationProtocol struct {
	e2s EncToShareProtocol
	s2e ShareToEncProtocol

	noise ring.DistributionParameters

	defaultScale *big.Int
	prec         uint

	mask    []*big.Int
	encoder *ckks.Encoder
}

// ShallowCopy creates a shallow copy of [MaskedLinearTransformationProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [MaskedLinearTransformationProtocol] can be used concurrently.
func (mltp MaskedLinearTransformationProtocol) ShallowCopy() MaskedLinearTransformationProtocol {

	mask := make([]*big.Int, mltp.e2s.params.N())
	for i := range mask {
		mask[i] = new(big.Int)
	}

	return MaskedLinearTransformationProtocol{
		e2s:          mltp.e2s.ShallowCopy(),
		s2e:          mltp.s2e.ShallowCopy(),
		prec:         mltp.prec,
		defaultScale: mltp.defaultScale,
		mask:         mask,
		encoder:      mltp.encoder.ShallowCopy(),
	}
}

// WithParams creates a shallow copy of the target [MaskedLinearTransformationProtocol] but with new output parameters.
// The expected input parameters remain unchanged.
func (mltp MaskedLinearTransformationProtocol) WithParams(paramsOut ckks.Parameters) MaskedLinearTransformationProtocol {

	s2e, err := NewShareToEncProtocol(paramsOut, mltp.noise)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	mask := make([]*big.Int, mltp.e2s.params.N())
	for i := range mask {
		mask[i] = new(big.Int)
	}

	scale := paramsOut.DefaultScale().Value

	defaultScale, _ := new(big.Float).SetPrec(mltp.prec).Set(&scale).Int(nil)

	return MaskedLinearTransformationProtocol{
		e2s:          mltp.e2s.ShallowCopy(),
		s2e:          s2e,
		prec:         mltp.prec,
		defaultScale: defaultScale,
		mask:         mask,
		encoder:      ckks.NewEncoder(paramsOut, mltp.prec),
	}
}

// MaskedLinearTransformationFunc represents a user-defined in-place function that can be evaluated on masked float plaintexts, as a part of the
// Masked Transform Protocol.
// The function is called with a vector of *Complex modulo ckks.Parameters.Slots() as input, and must write
// its output on the same buffer.
// Transform can be the identity.
//
//   - Decode: if true, then the masked float plaintext will be decoded before applying Transform.
//   - Recode: if true, then the masked float plaintext will be recoded after applying Transform.
//
// Decode (true/false) -> Transform -> Recode (true/false).
type MaskedLinearTransformationFunc struct {
	Decode bool
	Func   func(coeffs []*bignum.Complex)
	Encode bool
}

// NewMaskedLinearTransformationProtocol creates a new instance of the PermuteProtocol.
// paramsIn: the ckks.Parameters of the ciphertext before the protocol.
// paramsOut: the ckks.Parameters of the ciphertext after the protocol.
// prec : the log2 of decimal precision of the internal encoder.
// The method will return an error if the maximum number of slots of the output parameters is smaller than the number of slots of the input ciphertext.
func NewMaskedLinearTransformationProtocol(paramsIn, paramsOut ckks.Parameters, prec uint, noise ring.DistributionParameters) (mltp MaskedLinearTransformationProtocol, err error) {

	mltp = MaskedLinearTransformationProtocol{}

	mltp.noise = noise

	if mltp.e2s, err = NewEncToShareProtocol(paramsIn, noise); err != nil {
		return
	}

	if mltp.s2e, err = NewShareToEncProtocol(paramsOut, noise); err != nil {
		return
	}

	mltp.prec = prec

	scale := paramsOut.DefaultScale().Value

	mltp.defaultScale, _ = new(big.Float).SetPrec(prec).Set(&scale).Int(nil)

	mltp.mask = make([]*big.Int, paramsIn.N())
	for i := range mltp.mask {
		mltp.mask[i] = new(big.Int)
	}

	mltp.encoder = ckks.NewEncoder(paramsOut, prec)

	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (mltp MaskedLinearTransformationProtocol) AllocateShare(levelDecrypt, levelRecrypt int) multiparty.RefreshShare {
	return multiparty.RefreshShare{EncToShareShare: mltp.e2s.AllocateShare(levelDecrypt), ShareToEncShare: mltp.s2e.AllocateShare(levelRecrypt)}
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string. The CRP is considered to be in the NTT domain.
func (mltp MaskedLinearTransformationProtocol) SampleCRP(level int, crs sampling.PRNG) multiparty.KeySwitchCRP {
	return mltp.s2e.SampleCRP(level, crs)
}

// GenShare generates the shares of the PermuteProtocol
// This protocol requires additional inputs which are:
//
//   - skIn     : the secret-key if the input ciphertext.
//   - skOut    : the secret-key of the output ciphertext.
//   - logBound : the bit length of the masks.
//   - ct1      : the degree 1 element the ciphertext to refresh, i.e. ct1 = ckk.Ciphetext.Value[1].
//   - scale    : the scale of the ciphertext when entering the refresh.
//
// The method [GetMinimumLevelForRefresh] should be used to get the minimum level at which the masked transform can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (mltp MaskedLinearTransformationProtocol) GenShare(skIn, skOut *rlwe.SecretKey, logBound uint, ct *rlwe.Ciphertext, crs multiparty.KeySwitchCRP, transform *MaskedLinearTransformationFunc, shareOut *multiparty.RefreshShare) (err error) {

	ct1 := ct.Value[1]

	if ct1.Level() < shareOut.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot GenShare: ct[1] level must be at least equal to EncToShareShare level")
	}

	if crs.Value.Level() != shareOut.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot GenShare: crs level must be equal to ShareToEncShare")
	}

	if transform != nil {

		if transform.Decode && !ct.IsBatched {
			return fmt.Errorf("cannot GenShare: trying to decode a non-batched ciphertext (transform.Decode = true but ciphertext.IsBatched = false)")
		}

		if transform.Encode && !transform.Decode && ct.IsBatched {
			return fmt.Errorf("cannot GenShare: trying to encode a batched ciphertext (transform.Decode = false, transform.Encode = true but ciphertext.IsBatched = true")
		}
	}

	dslots := ct.Slots()
	if mltp.e2s.params.RingType() == ring.Standard {
		dslots *= 2
	}

	mask := mltp.mask[:dslots]

	// Generates the decryption share
	// Returns [M_i] on mltp.tmpMask and [a*s_i -M_i + e] on EncToShareShare
	if err = mltp.e2s.GenShare(skIn, logBound, ct, &multiparty.AdditiveShareBigint{Value: mask}, &shareOut.EncToShareShare); err != nil {
		return
	}

	// Applies LT(M_i)
	if err = mltp.applyTransformAndScale(transform, *ct.MetaData, mask); err != nil {
		return
	}

	// Stores the metadata of the ciphertext
	shareOut.MetaData = *ct.MetaData

	// Returns [-a*s_i + LT(M_i) * diffscale + e] on ShareToEncShare
	return mltp.s2e.GenShare(skOut, crs, ct.MetaData, multiparty.AdditiveShareBigint{Value: mask}, &shareOut.ShareToEncShare)
}

// AggregateShares sums share1 and share2 on shareOut.
func (mltp MaskedLinearTransformationProtocol) AggregateShares(share1, share2, shareOut *multiparty.RefreshShare) (err error) {

	if share1.EncToShareShare.Value.Level() != share2.EncToShareShare.Value.Level() || share1.EncToShareShare.Value.Level() != shareOut.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot AggregateShares: all e2s shares must be at the same level")
	}

	if share1.ShareToEncShare.Value.Level() != share2.ShareToEncShare.Value.Level() || share1.ShareToEncShare.Value.Level() != shareOut.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot AggregateShares: all s2e shares must be at the same level")
	}

	mltp.e2s.params.RingQ().AtLevel(share1.EncToShareShare.Value.Level()).Add(share1.EncToShareShare.Value, share2.EncToShareShare.Value, shareOut.EncToShareShare.Value)
	mltp.s2e.params.RingQ().AtLevel(share1.ShareToEncShare.Value.Level()).Add(share1.ShareToEncShare.Value, share2.ShareToEncShare.Value, shareOut.ShareToEncShare.Value)

	return
}

// Transform decrypts the ciphertext to LSSS-shares, applies the linear transformation on the LSSS-shares and re-encrypts the LSSS-shares to an RLWE ciphertext.
// The re-encrypted ciphertext's scale is set to the default scaling factor of the output parameters.
func (mltp MaskedLinearTransformationProtocol) Transform(ct *rlwe.Ciphertext, transform *MaskedLinearTransformationFunc, crs multiparty.KeySwitchCRP, share multiparty.RefreshShare, ciphertextOut *rlwe.Ciphertext) (err error) {

	if ct.Level() < share.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot Transform: input ciphertext level must be at least equal to e2s level")
	}

	if !ct.MetaData.Equal(&share.MetaData) {
		return fmt.Errorf("cannot Transform: input ciphertext MetaData is not equal to share.MetaData")
	}

	maxLevel := crs.Value.Level()

	if maxLevel != share.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot Transform: crs level and s2e level must be the same")
	}

	if transform != nil {

		if transform.Decode && !ct.IsBatched {
			return fmt.Errorf("cannot Transform: trying to decode a non-batched ciphertext (transform.Decode = true but ciphertext.IsBatched = false)")
		}

		if transform.Encode && !transform.Decode && ct.IsBatched {
			return fmt.Errorf("cannot Transform: trying to encode a batched ciphertext (transform.Decode = false, transform.Encode = true but ciphertext.IsBatched = true")
		}
	}

	ringQ := mltp.s2e.params.RingQ().AtLevel(maxLevel)

	dslots := ct.Slots()
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	mask := mltp.mask[:dslots]

	// Returns -sum(M_i) + x (outside of the NTT domain)
	mltp.e2s.GetShare(nil, share.EncToShareShare, ct, &multiparty.AdditiveShareBigint{Value: mask})

	// Returns LT(-sum(M_i) + x)
	if err = mltp.applyTransformAndScale(transform, *ct.MetaData, mask); err != nil {
		return
	}

	// Extend the levels of the ciphertext for future allocation
	if ciphertextOut.Value[0].N() != ringQ.N() {
		for i := range ciphertextOut.Value {
			ciphertextOut.Value[i] = ringQ.NewPoly()
		}
	} else {
		ciphertextOut.Resize(ciphertextOut.Degree(), maxLevel)
	}

	// Updates the ciphertext metadata if the output dimensions is smaller
	if logSlots := mltp.s2e.params.LogMaxSlots(); logSlots < ct.LogSlots() {
		ct.LogDimensions.Cols = logSlots
	}

	// Sets LT(-sum(M_i) + x) * diffscale in the RNS domain
	// Positional -> RNS -> NTT
	ringQ.SetCoefficientsBigint(mask, ciphertextOut.Value[0])
	rlwe.NTTSparseAndMontgomery(ringQ, ct.MetaData, ciphertextOut.Value[0])

	// LT(-sum(M_i) + x) * diffscale + [-a*s + LT(M_i) * diffscale + e] = [-a*s + LT(x) * diffscale + e]
	ringQ.Add(ciphertextOut.Value[0], share.ShareToEncShare.Value, ciphertextOut.Value[0])

	// Copies the result on the out ciphertext
	if err = mltp.s2e.GetEncryption(multiparty.KeySwitchShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut); err != nil {
		return
	}

	*ciphertextOut.MetaData = *ct.MetaData

	if transform != nil {
		ciphertextOut.IsBatched = transform.Encode
	}

	ciphertextOut.Scale = mltp.s2e.params.DefaultScale()

	return
}

func (mltp MaskedLinearTransformationProtocol) applyTransformAndScale(transform *MaskedLinearTransformationFunc, metadata rlwe.MetaData, mask []*big.Int) (err error) {

	slots := metadata.Slots()

	if transform != nil {

		bigComplex := make([]*bignum.Complex, slots)

		for i := range bigComplex {
			bigComplex[i] = bignum.NewComplex()
			bigComplex[i][0].SetPrec(mltp.prec)
			bigComplex[i][1].SetPrec(mltp.prec)
		}

		// Extracts sparse coefficients
		for i := 0; i < slots; i++ {
			bigComplex[i][0].SetInt(mask[i])
		}

		switch mltp.e2s.params.RingType() {
		case ring.Standard:
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				bigComplex[i][1].SetInt(mask[j])
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
			if err := mltp.encoder.FFT(bigComplex, metadata.LogSlots()); err != nil {
				return err
			}
		}

		// Applies the linear transform
		transform.Func(bigComplex)

		// Recodes if asked to
		if transform.Encode {
			if err := mltp.encoder.IFFT(bigComplex, metadata.LogSlots()); err != nil {
				return err
			}
		}

		// Puts the coefficient back
		for i := 0; i < slots; i++ {
			bigComplex[i].Real().Int(mask[i])
		}

		if mltp.e2s.params.RingType() == ring.Standard {
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				bigComplex[i].Imag().Int(mask[j])
			}
		}
	}

	// Applies LT(M_i) * diffscale
	inputScaleInt, d := new(big.Float).SetPrec(256).Set(&metadata.Scale.Value).Int(nil)

	// .Int truncates (i.e. does not round to the nearest integer)
	// Thus we check if we are below, and if yes add 1, which acts as rounding to the nearest integer
	if d == big.Below {
		inputScaleInt.Add(inputScaleInt, new(big.Int).SetInt64(1))
	}

	// Scales the mask by the ratio between the two scales
	for i := range mask {
		mask[i].Mul(mask[i], mltp.defaultScale)
		mask[i].Quo(mask[i], inputScaleInt)
	}

	return
}
