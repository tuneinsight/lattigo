package dckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// E2SProtocol is the structure storing the parameters and temporary buffers
// required by the encryption-to-shares protocol.
type E2SProtocol struct {
	CKSProtocol

	ringQ      *ring.Ring
	zero       *rlwe.SecretKey
	maskBigint []*big.Int
	pool       *ring.Poly
}

// NewE2SProtocol creates a new E2SProtocol struct from the passed CKKS parameters.
func NewE2SProtocol(params ckks.Parameters, sigmaSmudging float64) *E2SProtocol {
	e2s := new(E2SProtocol)
	e2s.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	e2s.ringQ = params.RingQ()
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.maskBigint = make([]*big.Int, params.Slots())
	e2s.pool = e2s.ringQ.NewPoly()
	return e2s
}

// GenShare generates a party's share in the encryption-to-shares protocol. This share consist in the additive secret-share of the party
// which is written in secretShareOut and in the public masked-decryption share written in publicShareOut.
func (e2s *E2SProtocol) GenShare(sk *rlwe.SecretKey, nbParties int, ct *ckks.Ciphertext, secretShareOut rlwe.AdditiveShare, publicShareOut *drlwe.CKSShare) {

	// Generates the share
	e2s.CKSProtocol.GenShare(sk, e2s.zero, ct, publicShareOut)

	// Generates the mask
	ringQ := e2s.ringQ

	// Get the upperbound on the norm
	bound := ring.NewUint(ringQ.Modulus[0])
	for i := 1; i < ct.Level()+1; i++ {
		bound.Mul(bound, ring.NewUint(ringQ.Modulus[i]))
	}

	// Divide the upperbound by 2*#Parties
	bound.Quo(bound, ring.NewUint(uint64(2*nbParties)))

	boundHalf := new(big.Int).Rsh(bound, 1)

	var sign int
	for i := range e2s.maskBigint {
		e2s.maskBigint[i] = ring.RandInt(bound)
		sign = e2s.maskBigint[i].Cmp(boundHalf)
		if sign == 1 || sign == 0 {
			e2s.maskBigint[i].Sub(e2s.maskBigint[i], bound)
		}
	}

	// h0 = mask (at level min)
	ringQ.SetCoefficientsBigintLvl(ct.Level(), e2s.maskBigint, e2s.pool)

	e2s.ringQ.CopyLvl(ct.Level(), e2s.pool, &secretShareOut.Value)
	e2s.ringQ.Sub(publicShareOut.Value, e2s.pool, publicShareOut.Value)
}

// GetShare is the final step of the encryption-to-share protocol. It performs the masked decryption of the target ciphertext followed by a
// the removal of the caller's secretShare as generated in the GenShare method.
// If the caller is not secret-key-share holder (i.e., didn't generate a decryption share), `secretShare` can be set to nil.
// Therefore, in order to obtain an additive sharing of the message, only one party should call this method, and the other parties should use
// the secretShareOut output of the GenShare method.
func (e2s *E2SProtocol) GetShare(secretShare *rlwe.AdditiveShare, aggregatePublicShare *drlwe.CKSShare, ct *ckks.Ciphertext, secretShareOut *rlwe.AdditiveShare) {
	e2s.ringQ.AddLvl(ct.Level(), aggregatePublicShare.Value, ct.Value[0], e2s.pool)
	if secretShare != nil {
		e2s.ringQ.AddLvl(ct.Level(), &secretShare.Value, e2s.pool, &secretShareOut.Value)
	} else {
		secretShareOut.Value.Copy(e2s.pool)
	}
}
