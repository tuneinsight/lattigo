package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKSecretKey is a type for ckks secret keys in a multi key context.
type MKSecretKey struct {
	key    *ckks.SecretKey
	peerID uint64
}

// MKPublicKey is a type for ckks public keys and ID in a multi key context. key[1] = a and key[0] = -s * a + e mod q
type MKPublicKey struct {
	key    [2]*MKDecomposedPoly
	peerID uint64
}

// MKDecomposedPoly is a type for vectors decomposed in a basis w (belong to Rq^d)(gadget decomposition)
type MKDecomposedPoly struct {
	poly []*ring.Poly
}

// MKEvaluationKey is a type for CKKS evaluation keys in a multi key context.
type MKEvaluationKey struct {
	key    []*MKDecomposedPoly
	peerID uint64
}

// MKSwitchingKey is a type for CKKS switching keys in a multi key context.
type MKSwitchingKey struct {
	key []*MKDecomposedPoly
}

// MKEvalGalKey is a type for CKKS rotation keys in a multi key context.
type MKEvalGalKey struct {
	key    []*MKDecomposedPoly
	peerID uint64
}

// MKRelinearizationKey is a type for CKKS relinearization keys in a multi key context.
type MKRelinearizationKey struct {
	key [][]*MKSwitchingKey
}

// MKKeys is a type that contains all keys necessary for the multi key protocol.
type MKKeys struct {
	secretKey *MKSecretKey
	publicKey *MKPublicKey
	evalKey   *MKEvaluationKey
	relinKey  *MKRelinearizationKey
}

// NewMKSwitchingKey allocate a MKSwitchingKey with zero polynomials in the ring r
func NewMKSwitchingKey(r *ring.Ring, params *ckks.Parameters) *MKSwitchingKey {

	key := new(MKSwitchingKey)
	key.key = make([]*MKDecomposedPoly, 3)

	key.key[0] = NewDecomposedPoly(r, params.Beta())
	key.key[1] = NewDecomposedPoly(r, params.Beta())
	key.key[2] = NewDecomposedPoly(r, params.Beta())

	return key
}

// NewMKEvaluationKey allocate a MKSwitchingKey with zero polynomials in the ring r adn with id = peerID
func NewMKEvaluationKey(r *ring.Ring, id uint64, params *ckks.Parameters) *MKEvaluationKey {

	key := new(MKEvaluationKey)
	key.key = make([]*MKDecomposedPoly, 3)

	key.key[0] = NewDecomposedPoly(r, params.Beta())
	key.key[1] = NewDecomposedPoly(r, params.Beta())
	key.key[2] = NewDecomposedPoly(r, params.Beta())

	key.peerID = id
	return key
}

// NewDecomposedPoly allocate a MKDecomposedPoly with zero polynomials in the ring r
func NewDecomposedPoly(r *ring.Ring, size uint64) *MKDecomposedPoly {

	res := new(MKDecomposedPoly)
	res.poly = make([]*ring.Poly, size)

	for i := uint64(0); i < size; i++ {
		res.poly[i] = r.NewPoly()
	}

	return res
}
