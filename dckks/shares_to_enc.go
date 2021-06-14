package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math/big"
)

// S2EProtocol is the structure storing the parameters and temporary buffers
// required by the shares-to-encryption protocol.
type S2EProtocol struct {
	CKSProtocol
	ringQ    *ring.Ring
	tmp      *ring.Poly
	ssBigint []*big.Int
	zero     *rlwe.SecretKey
}

// NewS2EProtocol creates a new S2EProtocol struct from the passed CKKS parameters.
func NewS2EProtocol(params ckks.Parameters, sigmaSmudging float64) *S2EProtocol {
	s2e := new(S2EProtocol)
	s2e.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	s2e.ringQ = params.RingQ()
	s2e.tmp = s2e.ringQ.NewPoly()
	s2e.ssBigint = make([]*big.Int, s2e.ringQ.N)
	s2e.zero = rlwe.NewSecretKey(params.Parameters)
	return s2e
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `c1` and the party's secret share of the message.
func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, c1 *ring.Poly, secretShare rlwe.AdditiveShare, c0ShareOut *drlwe.CKSShare) {
	s2e.CKSProtocol.GenShare(s2e.zero, sk, &rlwe.Element{Value: []*ring.Poly{c1, c1}, IsNTT: true}, c0ShareOut)

	if c1.Level() != c0ShareOut.Value.Level() {
		panic("c1 and c0ShareOut level must be equal")
	}

	if secretShare.Value.Level() > c1.Level() {
		panic("secretShare level cannot be greater than c1 or c0ShareOut level")
	}

	ringQ := s2e.ringQ

	minLevel := utils.MinInt(c1.Level(), secretShare.Value.Level())
	maxLevel := utils.MaxInt(c1.Level(), secretShare.Value.Level())

	centerAndExtendBasisLargeNorm(minLevel, maxLevel, ringQ, &secretShare.Value, s2e.ssBigint, s2e.tmp)

	ringQ.AddLvl(maxLevel, c0ShareOut.Value, s2e.tmp, c0ShareOut.Value)
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// share in the protocol and with the common, CRS-sampled polynomial `c1`.
func (s2e *S2EProtocol) GetEncryption(c0Agg *drlwe.CKSShare, c1 *ring.Poly, ctOut ckks.Ciphertext) {

	if ctOut.Degree() != 1 {
		panic("ctOut must have degree 1.")
	}

	if c0Agg.Value.Level() != c1.Level() {
		panic("c0Agg level must be equal to c1 level")
	}

	if ctOut.Level() != c1.Level() {
		panic("ctOut level must be equal to c1 level")
	}

	ctOut.Value[0].Copy(c0Agg.Value)
	ctOut.Value[1].Copy(c1)
}
