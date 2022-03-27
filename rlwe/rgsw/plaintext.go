package rgsw

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"math/big"
)

// Plaintext stores an RGSW plaintext value.
type Plaintext struct {
	Value []*ring.Poly
}

// NewPlaintext creates a new RGSW plaintext fron value, which can be either uint64, int64 or *ring.Poly.
// Plaintext is returned in the NTT and Mongtomery domain.
func NewPlaintext(value interface{}, levelQ, levelP, logBase2, decompBIT int, ringQP ringqp.Ring) (pt *Plaintext) {

	ringQ := ringQP.RingQ

	pt = new(Plaintext)
	pt.Value = make([]*ring.Poly, decompBIT)

	switch el := value.(type) {
	case uint64:
		pt.Value[0] = ringQ.NewPolyLvl(levelQ)
		for i := range ringQ.Modulus[:levelQ+1] {
			pt.Value[0].Coeffs[i][0] = el
		}
	case int64:
		pt.Value[0] = ringQ.NewPolyLvl(levelQ)
		if el < 0 {
			for i, qi := range ringQ.Modulus[:levelQ+1] {
				pt.Value[0].Coeffs[i][0] = qi - uint64(-el)
			}
		} else {
			for i := range ringQ.Modulus[:levelQ+1] {
				pt.Value[0].Coeffs[i][0] = uint64(el)
			}
		}
	case *ring.Poly:
		pt.Value[0] = el.CopyNew()
	default:
		panic("unsupported type, must be wither uint64 or *ring.Poly")
	}

	var pBigInt *big.Int
	if levelP > -1 {
		ringP := ringQP.RingP
		if levelP == len(ringP.Modulus)-1 {
			pBigInt = ringP.ModulusBigint
		} else {
			P := ringP.Modulus
			pBigInt = new(big.Int).SetUint64(P[0])
			for i := 1; i < levelP+1; i++ {
				pBigInt.Mul(pBigInt, ring.NewUint(P[i]))
			}
		}
		ringQ.MulScalarBigintLvl(levelQ, pt.Value[0], pBigInt, pt.Value[0])
	}

	ringQ.NTTLvl(levelQ, pt.Value[0], pt.Value[0])
	ringQ.MFormLvl(levelQ, pt.Value[0], pt.Value[0])

	for i := 1; i < len(pt.Value); i++ {
		pt.Value[i] = pt.Value[0].CopyNew()
		ringQ.MulByPow2Lvl(levelQ, pt.Value[i], i*logBase2, pt.Value[i])
	}

	return
}
