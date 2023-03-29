package bfv

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

func (p *PowerBasis) Write(data []byte) (n int, err error) {
	p.PowerBasis = &rlwe.PowerBasis{}
	return p.PowerBasis.Write(data)
}

// GenPower generates the n-th power of the power basis,
// as well as all the necessary intermediate powers if
// they are not yet present.
func (p *PowerBasis) GenPower(n int, eval Evaluator) {

	if p.Value[n] == nil {

		// Computes the index required to compute the required ring evaluation
		var a, b int
		if n&(n-1) == 0 {
			a, b = n/2, n/2 // Necessary for optimal depth
		} else {
			// Maximize the number of odd terms
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)
		}

		// Recurses on the given indexes
		p.GenPower(a, eval)
		p.GenPower(b, eval)

		// Computes C[n] = C[a]*C[b]
		p.Value[n] = eval.MulNew(p.Value[a], p.Value[b])
		eval.Relinearize(p.Value[n], p.Value[n])
	}
}
