package ckks

import (
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

// PowerBasis is a struct storing powers of a ciphertext.
type PowerBasis struct {
	*rlwe.PowerBasis
}

// NewPowerBasis creates a new PowerBasis.
func NewPowerBasis(ct *rlwe.Ciphertext, basis polynomial.Basis) (p *PowerBasis) {
	return &PowerBasis{rlwe.NewPowerBasis(ct, basis)}
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

// GenPower recursively computes X^{n}.
// If lazy = true, the final X^{n} will not be relinearized.
// Previous non-relinearized X^{n} that are required to compute the target X^{n} are automatically relinearized.
// Scale sets the threshold for rescaling (ciphertext won't be rescaled if the rescaling operation would make the scale go under this threshold).
func (p *PowerBasis) GenPower(n int, lazy bool, scale rlwe.Scale, eval *Evaluator) (err error) {

	if p.Value[n] == nil {
		if err = p.genPower(n, lazy, scale, eval); err != nil {
			return fmt.Errorf("genpower: p.Value[%d]: %w", n, err)
		}

		if err = eval.Rescale(p.Value[n], scale, p.Value[n]); err != nil {
			return fmt.Errorf("genpower: p.Value[%d]: rescale: %w", n, err)
		}
	}

	return nil
}

func (p *PowerBasis) genPower(n int, lazy bool, scale rlwe.Scale, eval *Evaluator) (err error) {

	if p.Value[n] == nil {

		a, b := rlwe.SplitDegree(n)

		// Recurses on the given indexes
		isPow2 := n&(n-1) == 0
		if err = p.genPower(a, lazy && !isPow2, scale, eval); err != nil {
			return fmt.Errorf("genpower: p.Value[%d]: %w", a, err)
		}
		if err = p.genPower(b, lazy && !isPow2, scale, eval); err != nil {
			return fmt.Errorf("genpower: p.Value[%d]: %w", b, err)
		}

		// Computes C[n] = C[a]*C[b]
		if lazy {
			if p.Value[a].Degree() == 2 {
				eval.Relinearize(p.Value[a], p.Value[a])
			}

			if p.Value[b].Degree() == 2 {
				eval.Relinearize(p.Value[b], p.Value[b])
			}

			if err = eval.Rescale(p.Value[a], scale, p.Value[a]); err != nil {
				return fmt.Errorf("genpower: rescale: p.Value[%d]: %w", a, err)
			}

			if err = eval.Rescale(p.Value[b], scale, p.Value[b]); err != nil {
				return fmt.Errorf("genpower: rescale: p.Value[%d]: %w", b, err)
			}

			p.Value[n] = eval.MulNew(p.Value[a], p.Value[b])

		} else {

			if err = eval.Rescale(p.Value[a], scale, p.Value[a]); err != nil {
				return fmt.Errorf("genpower: rescale: p.Value[%d]: %w", a, err)
			}

			if err = eval.Rescale(p.Value[b], scale, p.Value[b]); err != nil {
				return fmt.Errorf("genpower: rescale: p.Value[%d]: %w", b, err)
			}

			p.Value[n] = eval.MulRelinNew(p.Value[a], p.Value[b])
		}

		if p.Basis == polynomial.Chebyshev {

			// Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
			c := a - b
			if c < 0 {
				c = -c
			}

			// Computes C[n] = 2*C[a]*C[b]
			eval.Add(p.Value[n], p.Value[n], p.Value[n])

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				eval.Add(p.Value[n], -1, p.Value[n])
			} else {
				// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
				if err = p.GenPower(c, lazy, scale, eval); err != nil {
					return fmt.Errorf("genpower: p.Value[%d]: %w", c, err)
				}

				eval.Sub(p.Value[n], p.Value[c], p.Value[n])
			}
		}
	}
	return
}
