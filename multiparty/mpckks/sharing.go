package mpckks

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// EncToShareProtocol is the structure storing the parameters and temporary buffers
// required by the encryption-to-shares protocol.
type EncToShareProtocol struct {
	multiparty.KeySwitchProtocol

	params     ckks.Parameters
	zero       *rlwe.SecretKey
	maskBigint []*big.Int
	buff       ring.Poly
}

func NewAdditiveShare(params ckks.Parameters, logSlots int) multiparty.AdditiveShareBigint {

	nValues := 1 << logSlots
	if params.RingType() == ring.Standard {
		nValues <<= 1
	}

	return multiparty.NewAdditiveShareBigint(nValues)
}

// ShallowCopy creates a shallow copy of [EncToShareProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [EncToShareProtocol] can be used concurrently.
func (e2s EncToShareProtocol) ShallowCopy() EncToShareProtocol {

	maskBigint := make([]*big.Int, len(e2s.maskBigint))
	for i := range maskBigint {
		maskBigint[i] = new(big.Int)
	}

	return EncToShareProtocol{
		KeySwitchProtocol: e2s.KeySwitchProtocol.ShallowCopy(),
		params:            e2s.params,
		zero:              e2s.zero,
		maskBigint:        maskBigint,
		buff:              e2s.params.RingQ().NewPoly(),
	}
}

// NewEncToShareProtocol creates a new EncToShareProtocol struct from the passed parameters.
func NewEncToShareProtocol(params ckks.Parameters, noise ring.DistributionParameters) (EncToShareProtocol, error) {
	e2s := EncToShareProtocol{}

	var err error
	if e2s.KeySwitchProtocol, err = multiparty.NewKeySwitchProtocol(params.Parameters, noise); err != nil {
		return EncToShareProtocol{}, err
	}

	e2s.params = params
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.maskBigint = make([]*big.Int, params.N())
	for i := range e2s.maskBigint {
		e2s.maskBigint[i] = new(big.Int)
	}
	e2s.buff = e2s.params.RingQ().NewPoly()
	return e2s, nil
}

// AllocateShare allocates a share of the EncToShare protocol
func (e2s EncToShareProtocol) AllocateShare(level int) (share multiparty.KeySwitchShare) {
	return e2s.KeySwitchProtocol.AllocateShare(level)
}

// GenShare generates a party's share in the encryption-to-shares protocol. This share consist in the additive secret-share of the party
// which is written in secretShareOut and in the public masked-decryption share written in publicShareOut.
// This protocol requires additional inputs which are:
//
//   - logBound : the bit length of the masks
//   - ct: the ciphertext to share
//
// publicShareOut is always returned in the NTT domain.
// The method [GetMinimumLevelForRefresh] should be used to get the minimum level at which EncToShare can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (e2s EncToShareProtocol) GenShare(sk *rlwe.SecretKey, logBound uint, ct *rlwe.Ciphertext, secretShareOut *multiparty.AdditiveShareBigint, publicShareOut *multiparty.KeySwitchShare) (err error) {

	levelQ := utils.Min(ct.Value[1].Level(), publicShareOut.Value.Level())

	ringQ := e2s.params.RingQ().AtLevel(levelQ)

	// Get the upperbound on the norm
	// Ensures that bound >= 2^{128+logbound}
	bound := bignum.NewInt(1)
	bound.Lsh(bound, uint(logBound))

	boundMax := new(big.Int).Set(ringQ.ModulusAtLevel[levelQ])

	var sign int

	sign = bound.Cmp(boundMax)

	if sign == 1 || bound.Cmp(boundMax) == 1 {
		return fmt.Errorf("cannot GenShare: ciphertext level is not large enough for refresh correctness")
	}

	boundHalf := new(big.Int).Rsh(bound, 1)

	prng, _ := sampling.NewPRNG()

	dslots := ct.Slots()
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	// Generate the mask in Z[Y] for Y = X^{N/(2*slots)}
	for i := 0; i < dslots; i++ {
		e2s.maskBigint[i] = bignum.RandInt(prng, bound)
		sign = e2s.maskBigint[i].Cmp(boundHalf)
		if sign == 1 || sign == 0 {
			e2s.maskBigint[i].Sub(e2s.maskBigint[i], bound)
		}

		secretShareOut.Value[i].Set(e2s.maskBigint[i])
	}

	// Encrypt the mask
	// Generates an encryption of zero and subtracts the mask
	e2s.KeySwitchProtocol.GenShare(sk, e2s.zero, ct, publicShareOut)

	// Positional -> RNS -> NTT
	ringQ.SetCoefficientsBigint(secretShareOut.Value[:dslots], e2s.buff)
	rlwe.NTTSparseAndMontgomery(ringQ, ct.MetaData, e2s.buff)

	// Subtracts the mask to the encryption of zero
	ringQ.Sub(publicShareOut.Value, e2s.buff, publicShareOut.Value)

	return
}

// GetShare is the final step of the encryption-to-share protocol. It performs the masked decryption of the target ciphertext followed by a
// the removal of the caller's secretShare as generated in the GenShare method.
// If the caller is not secret-key-share holder (i.e., didn't generate a decryption share), secretShare can be set to nil.
// Therefore, in order to obtain an additive sharing of the message, only one party should call this method, and the other parties should use
// the secretShareOut output of the GenShare method.
func (e2s EncToShareProtocol) GetShare(secretShare *multiparty.AdditiveShareBigint, aggregatePublicShare multiparty.KeySwitchShare, ct *rlwe.Ciphertext, secretShareOut *multiparty.AdditiveShareBigint) {

	levelQ := utils.Min(ct.Level(), aggregatePublicShare.Value.Level())

	ringQ := e2s.params.RingQ().AtLevel(levelQ)

	// Adds the decryption share on the ciphertext and stores the result in a buff
	ringQ.Add(aggregatePublicShare.Value, ct.Value[0], e2s.buff)

	// INTT -> RNS -> Positional
	ringQ.INTT(e2s.buff, e2s.buff)

	dslots := ct.Slots()
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	gap := ringQ.N() / dslots

	ringQ.PolyToBigintCentered(e2s.buff, gap, e2s.maskBigint)

	// Subtracts the last mask
	if secretShare != nil {
		a := secretShareOut.Value
		b := e2s.maskBigint
		c := secretShare.Value
		for i := range secretShareOut.Value[:dslots] {
			a[i].Add(c[i], b[i])
		}
	} else {
		a := secretShareOut.Value
		b := e2s.maskBigint
		for i := range secretShareOut.Value[:dslots] {
			a[i].Set(b[i])
		}
	}
}

// ShareToEncProtocol is the structure storing the parameters and temporary buffers
// required by the shares-to-encryption protocol.
type ShareToEncProtocol struct {
	multiparty.KeySwitchProtocol
	params   ckks.Parameters
	tmp      ring.Poly
	ssBigint []*big.Int
	zero     *rlwe.SecretKey
}

// ShallowCopy creates a shallow copy of [ShareToEncProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [ShareToEncProtocol] can be used concurrently.
func (s2e ShareToEncProtocol) ShallowCopy() ShareToEncProtocol {
	return ShareToEncProtocol{
		KeySwitchProtocol: s2e.KeySwitchProtocol.ShallowCopy(),
		params:            s2e.params,
		tmp:               s2e.params.RingQ().NewPoly(),
		ssBigint:          make([]*big.Int, s2e.params.N()),
		zero:              s2e.zero,
	}
}

// NewShareToEncProtocol creates a new ShareToEncProtocol struct from the passed parameters.
func NewShareToEncProtocol(params ckks.Parameters, noise ring.DistributionParameters) (ShareToEncProtocol, error) {
	s2e := ShareToEncProtocol{}

	var err error
	if s2e.KeySwitchProtocol, err = multiparty.NewKeySwitchProtocol(params.Parameters, noise); err != nil {
		return ShareToEncProtocol{}, err
	}

	s2e.params = params
	s2e.tmp = s2e.params.RingQ().NewPoly()
	s2e.ssBigint = make([]*big.Int, s2e.params.N())
	s2e.zero = rlwe.NewSecretKey(params.Parameters)
	return s2e, nil
}

// AllocateShare allocates a share of the ShareToEnc protocol
func (s2e ShareToEncProtocol) AllocateShare(level int) (share multiparty.KeySwitchShare) {
	return s2e.KeySwitchProtocol.AllocateShare(level)
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `crs` and the party's secret share of the message.
func (s2e ShareToEncProtocol) GenShare(sk *rlwe.SecretKey, crs multiparty.KeySwitchCRP, metadata *rlwe.MetaData, secretShare multiparty.AdditiveShareBigint, c0ShareOut *multiparty.KeySwitchShare) (err error) {

	if crs.Value.Level() != c0ShareOut.Value.Level() {
		return fmt.Errorf("cannot GenShare: crs and c0ShareOut level must be equal")
	}

	ringQ := s2e.params.RingQ().AtLevel(crs.Value.Level())

	// Generates an encryption share
	ct := &rlwe.Ciphertext{}
	ct.Value = []ring.Poly{{}, crs.Value}
	ct.MetaData = &rlwe.MetaData{}
	ct.IsNTT = true
	s2e.KeySwitchProtocol.GenShare(s2e.zero, sk, ct, c0ShareOut)

	dslots := metadata.Slots()
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	// Positional -> RNS -> NTT
	ringQ.SetCoefficientsBigint(secretShare.Value[:dslots], s2e.tmp)

	rlwe.NTTSparseAndMontgomery(ringQ, metadata, s2e.tmp)

	ringQ.Add(c0ShareOut.Value, s2e.tmp, c0ShareOut.Value)

	return
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// share in the protocol and with the common, CRS-sampled polynomial `crs`.
func (s2e ShareToEncProtocol) GetEncryption(c0Agg multiparty.KeySwitchShare, crs multiparty.KeySwitchCRP, opOut *rlwe.Ciphertext) (err error) {
	if opOut.Degree() != 1 {
		return fmt.Errorf("cannot GetEncryption: opOut must have degree 1")
	}

	if c0Agg.Value.Level() != crs.Value.Level() {
		return fmt.Errorf("cannot GetEncryption: c0Agg level must be equal to crs level")
	}

	if opOut.Level() != crs.Value.Level() {
		return fmt.Errorf("cannot GetEncryption: opOut level must be equal to crs level")
	}

	opOut.Value[0].Copy(c0Agg.Value)
	opOut.Value[1].Copy(crs.Value)

	return
}
