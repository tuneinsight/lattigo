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

// MaskedLinearTransformationProtocol is a struct storing the parameters for the MaskedLinearTransformationProtocol protocol.
type MaskedLinearTransformationProtocol struct {
	e2s EncToShareProtocol
	s2e ShareToEncProtocol

	noise ring.DistributionParameters

	defaultScale *big.Int
	prec         uint

	tmpMaskIn  []*big.Int
	tmpMaskOut []*big.Int
	encoder    *ckks.Encoder
}

// ShallowCopy creates a shallow copy of MaskedLinearTransformationProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// MaskedLinearTransformationProtocol can be used concurrently.
func (mltp MaskedLinearTransformationProtocol) ShallowCopy() MaskedLinearTransformationProtocol {

	tmpMaskIn := make([]*big.Int, mltp.e2s.params.N())
	for i := range tmpMaskIn {
		tmpMaskIn[i] = new(big.Int)
	}

	tmpMaskOut := make([]*big.Int, mltp.s2e.params.N())
	for i := range tmpMaskOut {
		tmpMaskOut[i] = new(big.Int)
	}

	return MaskedLinearTransformationProtocol{
		e2s:          mltp.e2s.ShallowCopy(),
		s2e:          mltp.s2e.ShallowCopy(),
		prec:         mltp.prec,
		defaultScale: mltp.defaultScale,
		tmpMaskIn:    tmpMaskIn,
		tmpMaskOut:   tmpMaskOut,
		encoder:      mltp.encoder.ShallowCopy(),
	}
}

// WithParams creates a shallow copy of the target MaskedLinearTransformationProtocol but with new output parameters.
// The expected input parameters remain unchanged.
func (mltp MaskedLinearTransformationProtocol) WithParams(paramsOut ckks.Parameters) MaskedLinearTransformationProtocol {

	s2e, err := NewShareToEncProtocol(paramsOut, mltp.noise)

	if err != nil {
		panic(err)
	}

	tmpMaskIn := make([]*big.Int, mltp.e2s.params.N())
	for i := range tmpMaskIn {
		tmpMaskIn[i] = new(big.Int)
	}

	tmpMaskOut := make([]*big.Int, mltp.s2e.params.N())
	for i := range tmpMaskOut {
		tmpMaskOut[i] = new(big.Int)
	}

	scale := paramsOut.DefaultScale().Value

	defaultScale, _ := new(big.Float).SetPrec(mltp.prec).Set(&scale).Int(nil)

	return MaskedLinearTransformationProtocol{
		e2s:          mltp.e2s.ShallowCopy(),
		s2e:          s2e,
		prec:         mltp.prec,
		defaultScale: defaultScale,
		tmpMaskIn:    tmpMaskIn,
		tmpMaskOut:   tmpMaskOut,
		encoder:      ckks.NewEncoder(paramsOut, mltp.prec),
	}
}

// MaskedLinearTransformationFunc represents a user-defined in-place function that can be evaluated on masked CKKS plaintexts, as a part of the
// Masked Transform Protocol.
// The function is called with a vector of *Complex modulo ckks.Parameters.Slots() as input, and must write
// its output on the same buffer.
// Transform can be the identity.
// Decode: if true, then the masked CKKS plaintext will be decoded before applying Transform.
// Recode: if true, then the masked CKKS plaintext will be recoded after applying Transform.
// i.e. : Decode (true/false) -> Transform -> Recode (true/false).
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

	mltp.tmpMaskIn = make([]*big.Int, paramsIn.N())
	for i := range mltp.tmpMaskIn {
		mltp.tmpMaskIn[i] = new(big.Int)
	}

	mltp.tmpMaskOut = make([]*big.Int, paramsOut.N())
	for i := range mltp.tmpMaskOut {
		mltp.tmpMaskOut[i] = new(big.Int)
	}

	mltp.encoder = ckks.NewEncoder(paramsOut, prec)

	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (mltp MaskedLinearTransformationProtocol) AllocateShare(levelDecrypt, levelRecrypt int) drlwe.RefreshShare {
	return drlwe.RefreshShare{EncToShareShare: mltp.e2s.AllocateShare(levelDecrypt), ShareToEncShare: mltp.s2e.AllocateShare(levelRecrypt)}
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string. The CRP is considered to be in the NTT domain.
func (mltp MaskedLinearTransformationProtocol) SampleCRP(level int, crs sampling.PRNG) drlwe.KeySwitchCRP {
	return mltp.s2e.SampleCRP(level, crs)
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
func (mltp MaskedLinearTransformationProtocol) GenShare(skIn, skOut *rlwe.SecretKey, logBound uint, ct *rlwe.Ciphertext, crs drlwe.KeySwitchCRP, transform *MaskedLinearTransformationFunc, shareOut *drlwe.RefreshShare) (err error) {

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

	// Generates the decryption share
	// Returns [M_i] on mltp.tmpMask and [a*s_i -M_i + e] on EncToShareShare
	if err = mltp.e2s.GenShare(skIn, logBound, ct, &drlwe.AdditiveShareBigint{Value: mltp.tmpMaskIn}, &shareOut.EncToShareShare); err != nil {
		return
	}

	// Changes ring if necessary:
	// X -> X or Y -> X or X -> Y for Y = X^(2^s)
	maskOut := mltp.changeRing(mltp.tmpMaskIn)

	// Applies LT(M_i)
	if err = mltp.applyTransformAndScale(transform, ct.Scale, maskOut); err != nil {
		return
	}

	// Returns [-a*s_i + LT(M_i) * diffscale + e] on ShareToEncShare
	return mltp.s2e.GenShare(skOut, crs, ct.MetaData, drlwe.AdditiveShareBigint{Value: maskOut}, &shareOut.ShareToEncShare)
}

// AggregateShares sums share1 and share2 on shareOut.
func (mltp MaskedLinearTransformationProtocol) AggregateShares(share1, share2, shareOut *drlwe.RefreshShare) (err error) {

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

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
// The ciphertext scale is reset to the default scale.
func (mltp MaskedLinearTransformationProtocol) Transform(ct *rlwe.Ciphertext, transform *MaskedLinearTransformationFunc, crs drlwe.KeySwitchCRP, share drlwe.RefreshShare, ciphertextOut *rlwe.Ciphertext) (err error) {

	if ct.Level() < share.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot Transform: input ciphertext level must be at least equal to e2s level")
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

	// Returns -sum(M_i) + x (outside of the NTT domain)

	mltp.e2s.GetShare(nil, share.EncToShareShare, ct, &drlwe.AdditiveShareBigint{Value: mltp.tmpMaskIn})

	// Changes ring if necessary:
	// X -> X or Y -> X or X -> Y for Y = X^(2^s)
	maskOut := mltp.changeRing(mltp.tmpMaskIn)

	// Returns LT(-sum(M_i) + x)
	if err = mltp.applyTransformAndScale(transform, ct.Scale, maskOut); err != nil {
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
	ringQ.SetCoefficientsBigint(maskOut, ciphertextOut.Value[0])
	ringQ.NTT(ciphertextOut.Value[0], ciphertextOut.Value[0])

	// LT(-sum(M_i) + x) * diffscale + [-a*s + LT(M_i) * diffscale + e] = [-a*s + LT(x) * diffscale + e]
	ringQ.Add(ciphertextOut.Value[0], share.ShareToEncShare.Value, ciphertextOut.Value[0])

	// Copies the result on the out ciphertext
	if err = mltp.s2e.GetEncryption(drlwe.KeySwitchShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut); err != nil {
		return
	}

	*ciphertextOut.MetaData = *ct.MetaData

	if transform != nil {
		ciphertextOut.IsBatched = transform.Encode
	}

	ciphertextOut.Scale = mltp.s2e.params.DefaultScale()

	return
}

func (mltp MaskedLinearTransformationProtocol) changeRing(maskIn []*big.Int) (maskOut []*big.Int) {

	NIn := mltp.e2s.params.N()
	NOut := mltp.s2e.params.N()

	if NIn == NOut {
		maskOut = maskIn
	} else if NIn < NOut {

		maskOut = mltp.tmpMaskOut

		gap := NOut / NIn

		for i := 0; i < NIn; i++ {
			maskOut[i*gap].Set(maskIn[i])
		}

	} else {

		maskOut = mltp.tmpMaskOut

		gap := NIn / NOut

		for i := 0; i < NOut; i++ {
			maskOut[i].Set(maskIn[i*gap])
		}
	}

	return
}

func (mltp MaskedLinearTransformationProtocol) applyTransformAndScale(transform *MaskedLinearTransformationFunc, inputScale rlwe.Scale, mask []*big.Int) (err error) {

	slots := mltp.s2e.params.MaxSlots()

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
			if err := mltp.encoder.FFT(bigComplex, mltp.s2e.params.LogMaxSlots()); err != nil {
				return err
			}
		}

		// Applies the linear transform
		transform.Func(bigComplex)

		// Recodes if asked to
		if transform.Encode {
			if err := mltp.encoder.IFFT(bigComplex, mltp.s2e.params.LogMaxSlots()); err != nil {
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
	inputScaleInt, _ := new(big.Float).SetPrec(256).Set(&inputScale.Value).Int(nil)

	// Scales the mask by the ratio between the two scales
	for i := range mask {
		mask[i].Mul(mask[i], mltp.defaultScale)
		mask[i].Quo(mask[i], inputScaleInt)
	}

	return
}
