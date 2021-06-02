package mkrlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MKSecretKey is a type for rlwe secret keys in a multi key context.
type MKSecretKey struct {
	Key    *rlwe.SecretKey
	PeerID uint64
}

// MKPublicKey is a type for rlwe public keys and ID in a multi key context. key[1] = a and key[0] = -s * a + e mod q
type MKPublicKey struct {
	Key    [2]*MKDecomposedPoly
	PeerID uint64
}

// MKDecomposedPoly is a type for vectors decomposed in a basis w (belong to Rq^d)(gadget decomposition)
type MKDecomposedPoly struct {
	Poly []*ring.Poly
}

// MKEvaluationKey is a type for evaluation keys in a multi key context.
type MKEvaluationKey struct {
	Key01  *rlwe.SwitchingKey
	Key2   *MKDecomposedPoly
	PeerID uint64
}

// MKEvalGalKey is a type for rotation keys in a multi key context.
type MKEvalGalKey struct {
	Key    *rlwe.SwitchingKey
	PeerID uint64
}

// MKKeys is a type that contains all keys necessary for the multi key protocol.
type MKKeys struct {
	SecretKey *MKSecretKey
	PublicKey *MKPublicKey
	EvalKey   *MKEvaluationKey
}

// NewMKEvaluationKey allocate a MKSwitchingKey with zero polynomials in the ring r adn with id = peerID
func NewMKEvaluationKey(r *ring.Ring, id uint64, params *rlwe.Parameters) *MKEvaluationKey {

	key := new(MKEvaluationKey)
	key.Key01 = rlwe.NewSwitchingKey(*params)      //D0 and D1 in Chen et Al's notation
	key.Key2 = NewDecomposedPoly(r, params.Beta()) // D2 in Chen et Al's notation
	key.PeerID = id

	return key
}

// NewDecomposedPoly allocate a MKDecomposedPoly with zero polynomials in the ring r
func NewDecomposedPoly(r *ring.Ring, size uint64) *MKDecomposedPoly {

	res := new(MKDecomposedPoly)
	res.Poly = make([]*ring.Poly, size)

	for i := uint64(0); i < size; i++ {
		res.Poly[i] = r.NewPoly()
	}

	return res
}

// CopyNewDecomposed copy a decomposedd polynomial and return the copy
func CopyNewDecomposed(p *MKDecomposedPoly) *MKDecomposedPoly {

	res := new(MKDecomposedPoly)
	res.Poly = make([]*ring.Poly, len(p.Poly))

	for i := range p.Poly {
		res.Poly[i] = p.Poly[i].CopyNew()
	}

	return res
}

// FromDecomposedToSwitchingKey converts two decomposed Poly into an rlwe switching key
func FromDecomposedToSwitchingKey(p1 *MKDecomposedPoly, p2 *MKDecomposedPoly, params *rlwe.Parameters) *rlwe.SwitchingKey {

	decompSize := params.Beta()

	if len(p1.Poly) != int(decompSize) {
		panic("Decomposed Poly should have same size as decomposition basis size")
	}
	swk := new(rlwe.SwitchingKey)
	swk.Value = make([][2]*ring.Poly, int(decompSize))
	for i := uint64(0); i < decompSize; i++ {
		swk.Value[i][0] = p1.Poly[i]
		swk.Value[i][1] = p2.Poly[i]
	}
	return swk
}
