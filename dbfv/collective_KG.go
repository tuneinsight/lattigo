//Package dbfv implements a distributed (or threshold) version of the BFV scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// CKG is the structure storing the parameters for the collective key generation protocol.
type CKG struct {
	context         *ring.Context
	gaussianSampler *ring.KYSampler

	share *ring.Poly
	cpk   [2]*ring.Poly
}

// NewCKG creates a new CKG instance that will be used to generate a
// collective public key among j parties, using a common reference poylnomial.
func NewCKG(context *ring.Context, crs *ring.Poly) *CKG {
	ckg := new(CKG)
	ckg.context = context
	ckg.gaussianSampler = context.NewKYSampler(3.19, 19)

	ckg.cpk[0] = context.NewPoly()
	ckg.cpk[1] = crs.CopyNew()

	return ckg
}

// GenShare is the first and unique round of the CKG protocol. Each party generates a secret share
// and computes from it a public-share of the form :
//
// [a*s_i + e_i, a]
//
// and broadcasts it to all other j-1 parties.
func (ckg *CKG) GenShare(sk *ring.Poly) error {

	// -(sk * crs) + e
	ckg.share = ckg.gaussianSampler.SampleNTTNew()
	ckg.context.MulCoeffsMontgomeryAndSub(sk, ckg.cpk[1], ckg.share)

	return nil
}

// GetShare returns the public-share stored in the CKG structure.
func (ckg *CKG) GetShare() *ring.Poly {
	return ckg.share
}

// AggregateShares is the second part of the first and unique round of the CKS protocol. Uppon receiving the j-1 shares,
// each party computes :
//
// [sum(a* s_i + e_i), sum(a)] = [b * s + e, b]
func (ckg *CKG) AggregateShares(shares []*ring.Poly) error {

	for i := 0; i < len(shares); i++ {
		ckg.context.Add(ckg.cpk[0], shares[i], ckg.cpk[0])
	}

	return nil
}

// Finalize transforms the aggregated shares into a new public-key element and returns it.
func (ckg *CKG) Finalize() (*bfv.PublicKey, error) {
	collectivePk := new(bfv.PublicKey)
	collectivePk.Set(ckg.cpk)
	return collectivePk, nil
}
