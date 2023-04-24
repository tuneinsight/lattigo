package bgv

import (
	"io"
	"math"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

// PowerBasis is a struct storing powers of a ciphertext.
type PowerBasis struct {
	*rlwe.PowerBasis
}

// NewPowerBasis creates a new PowerBasis.
func NewPowerBasis(ct *rlwe.Ciphertext) (p *PowerBasis) {
	return &PowerBasis{rlwe.NewPowerBasis(ct, polynomial.Monomial)}
}

func (p *PowerBasis) UnmarshalBinary(data []byte) (err error) {
	p.PowerBasis = &rlwe.PowerBasis{}
	return p.PowerBasis.UnmarshalBinary(data)
}

func (p *PowerBasis) ReadFrom(r io.Reader) (n int64, err error) {
	p.PowerBasis = &rlwe.PowerBasis{}
	return p.PowerBasis.ReadFrom(r)
}

func (p *PowerBasis) Decode(data []byte) (n int, err error) {
	p.PowerBasis = &rlwe.PowerBasis{}
	return p.PowerBasis.Decode(data)
}

// GenPower generates the n-th power of the power basis,
// as well as all the necessary intermediate powers if
// they are not yet present.
func (p *PowerBasis) GenPower(n int, lazy bool, eval Evaluator) (err error) {

	var rescale bool
	if rescale, err = p.genPower(n, n, lazy, true, eval); err != nil {
		return
	}

	if rescale {
		if err = eval.Rescale(p.Value[n], p.Value[n]); err != nil {
			return
		}
	}

	return nil
}

func (p *PowerBasis) genPower(target, n int, lazy, rescale bool, eval Evaluator) (rescaleN bool, err error) {

	if p.Value[n] == nil {

		isPow2 := n&(n-1) == 0

		// Computes the index required to compute the required ring evaluation
		var a, b int
		if isPow2 {
			a, b = n/2, n/2 // Necessary for optimal depth
		} else {
			// Maximize the number of odd terms
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)
		}

		var rescaleA, rescaleB bool

		// Recurses on the given indexes
		if rescaleA, err = p.genPower(target, a, lazy, rescale, eval); err != nil {
			return false, err
		}

		if rescaleB, err = p.genPower(target, b, lazy, rescale, eval); err != nil {
			return false, err
		}

		if p.Value[a].Degree() == 2 {
			eval.Relinearize(p.Value[a], p.Value[a])
		}

		if p.Value[b].Degree() == 2 {
			eval.Relinearize(p.Value[b], p.Value[b])
		}

		if rescaleA {
			if err = eval.Rescale(p.Value[a], p.Value[a]); err != nil {
				return false, err
			}
		}

		if rescaleB {
			if err = eval.Rescale(p.Value[b], p.Value[b]); err != nil {
				return false, err
			}
		}

		// Computes C[n] = C[a]*C[b]
		if lazy && !isPow2 {
			p.Value[n] = eval.MulNew(p.Value[a], p.Value[b])
			return true, nil
		}

		p.Value[n] = eval.MulRelinNew(p.Value[a], p.Value[b])
		if err = eval.Rescale(p.Value[n], p.Value[n]); err != nil {
			return false, err
		}

	}

	return false, nil
}
