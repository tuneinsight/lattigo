package polynomial

import (
	"bufio"
	"fmt"
	"io"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/buffer"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// PowerBasis is a struct storing powers of a ciphertext.
type PowerBasis struct {
	bignum.Basis
	Value structs.Map[int, rlwe.Ciphertext]
}

// NewPowerBasis creates a new [PowerBasis]. It takes as input a ciphertext
// and a basis type. The struct treats the input ciphertext as a monomial X and
// can be used to generates power of this monomial X^{n} in the given [BasisType].
func NewPowerBasis(ct *rlwe.Ciphertext, basis bignum.Basis) (p PowerBasis) {
	return PowerBasis{
		Value: map[int]*rlwe.Ciphertext{1: ct.CopyNew()},
		Basis: basis,
	}
}

// SplitDegree returns a * b = n such that |a-b| is minimized
// with a and/or b odd if possible.
func SplitDegree(n int) (a, b int) {

	if n&(n-1) == 0 {
		a, b = n/2, n/2 //Necessary for optimal depth
	} else {
		// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
		// Maximize the number of odd terms of Chebyshev basis
		k := bits.Len64(uint64(n-1)) - 1
		a = (1 << k) - 1
		b = n + 1 - (1 << k)
	}

	return
}

// GenPower recursively computes X^{n}.
// If lazy = true, the final X^{n} will not be relinearized.
// Previous non-relinearized X^{n} that are required to compute the target X^{n} are automatically relinearized.
func (p *PowerBasis) GenPower(n int, lazy bool, eval schemes.Evaluator) (err error) {

	if eval == nil {
		return fmt.Errorf("cannot GenPower: EvaluatorInterface is nil")
	}

	if p.Value[n] == nil {

		var rescale bool
		if rescale, err = p.genPower(n, lazy, true, eval); err != nil {
			return fmt.Errorf("genpower: p.Value[%d]: %w", n, err)
		}

		if rescale {
			if err = eval.Rescale(p.Value[n], p.Value[n]); err != nil {
				return fmt.Errorf("genpower: p.Value[%d]: final rescale: %w", n, err)
			}
		}
	}

	return nil
}

func (p *PowerBasis) genPower(n int, lazy, rescale bool, eval schemes.Evaluator) (rescaltOut bool, err error) {

	if p.Value[n] == nil {

		a, b := SplitDegree(n)

		// Recurses on the given indexes
		isPow2 := n&(n-1) == 0

		var rescaleA, rescaleB bool // Avoids calling rescale on already generated powers

		if rescaleA, err = p.genPower(a, lazy && !isPow2, rescale, eval); err != nil {
			return false, fmt.Errorf("genpower: p.Value[%d]: %w", a, err)
		}
		if rescaleB, err = p.genPower(b, lazy && !isPow2, rescale, eval); err != nil {
			return false, fmt.Errorf("genpower: p.Value[%d]: %w", b, err)
		}

		// Computes C[n] = C[a]*C[b]
		if lazy {

			if p.Value[a].Degree() == 2 {
				if err = eval.Relinearize(p.Value[a], p.Value[a]); err != nil {
					return false, fmt.Errorf("genpower (lazy): eval.Relinearize(p.Value[%d], p.Value[%d]): %w", a, a, err)
				}
			}

			if p.Value[b].Degree() == 2 {
				if err = eval.Relinearize(p.Value[b], p.Value[b]); err != nil {
					return false, fmt.Errorf("genpower (lazy): eval.Relinearize(p.Value[%d], p.Value[%d]): %w", b, b, err)
				}
			}

			if rescaleA {
				if err = eval.Rescale(p.Value[a], p.Value[a]); err != nil {
					return false, fmt.Errorf("genpower (lazy): rescale[a]: p.Value[%d]: %w", a, err)
				}
			}

			if rescaleB {
				if err = eval.Rescale(p.Value[b], p.Value[b]); err != nil {
					return false, fmt.Errorf("genpower (lazy): rescale[b]: p.Value[%d]: %w", b, err)
				}
			}

			if p.Value[n], err = eval.MulNew(p.Value[a], p.Value[b]); err != nil {
				return false, fmt.Errorf("genpower (lazy): Mulnew(p.Value[%d], p.Value[%d]): %w", a, b, err)
			}

		} else {

			if rescaleA {
				if err = eval.Rescale(p.Value[a], p.Value[a]); err != nil {
					return false, fmt.Errorf("genpower: rescale[a]: p.Value[%d]: %w", a, err)
				}
			}

			if rescaleB {
				if err = eval.Rescale(p.Value[b], p.Value[b]); err != nil {
					return false, fmt.Errorf("genpower: rescale[b]: p.Value[%d]: %w", b, err)
				}
			}

			if p.Value[n], err = eval.MulRelinNew(p.Value[a], p.Value[b]); err != nil {
				return false, fmt.Errorf("genpower: MulRelinNew(p.Value[%d], p.Value[%d]): %w", a, b, err)
			}
		}

		if p.Basis == bignum.Chebyshev {

			// Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
			c := a - b
			if c < 0 {
				c = -c
			}

			// Computes C[n] = 2*C[a]*C[b]
			if err = eval.Add(p.Value[n], p.Value[n], p.Value[n]); err != nil {
				return false, fmt.Errorf("genpower: Add(p.Value[%d], p.Value[%d], p.Value[%d]): %w", n, n, n, err)
			}

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				if err = eval.Add(p.Value[n], -1, p.Value[n]); err != nil {
					return false, fmt.Errorf("genpower: Add(p.Value[%d], -1, p.Value[%d]): %w", n, n, err)
				}
			} else {
				// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
				if err = p.GenPower(c, lazy, eval); err != nil {
					return false, fmt.Errorf("genpower: p.Value[%d]: %w", c, err)
				}

				if err = eval.Sub(p.Value[n], p.Value[c], p.Value[n]); err != nil {
					return false, fmt.Errorf("genpower: Add(p.Value[%d], p.Value[%d], p.Value[%d]): %w", n, c, n, err)
				}
			}
		}

		return true, nil
	}

	return false, nil
}

// BinarySize returns the serialized size of the object in bytes.
func (p PowerBasis) BinarySize() (size int) {
	return 1 + p.Value.BinarySize()
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the [buffer.Writer] interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a [bufio.Writer]. Since this requires allocations, it
// is preferable to pass a [buffer.Writer] directly:
//
//   - When writing multiple times to a [io.Writer], it is preferable to first wrap the
//     [io.Writer] in a pre-allocated [bufio.Writer].
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (p PowerBasis) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		var inc int64

		if inc, err = buffer.WriteUint8(w, uint8(p.Basis)); err != nil {
			return n + inc, err
		}

		n += inc

		inc, err = p.Value.WriteTo(w)

		return n + inc, err

	default:
		return p.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a [io.Reader], it is preferable to first
//     first wrap [io.Reader] in a pre-allocated [bufio.Reader].
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (p *PowerBasis) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:
		var inc int64

		var Basis uint8

		if inc, err = buffer.ReadUint8(r, &Basis); err != nil {
			return n + inc, err
		}

		n += inc

		p.Basis = bignum.Basis(Basis)

		if p.Value == nil {
			p.Value = map[int]*rlwe.Ciphertext{}
		}

		inc, err = p.Value.ReadFrom(r)

		return n + inc, err

	default:
		return p.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (p PowerBasis) MarshalBinary() (data []byte, err error) {
	buf := buffer.NewBufferSize(p.BinarySize())
	_, err = p.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// [PowerBasis.MarshalBinary] or [PowerBasis.WriteTo] on the object.
func (p *PowerBasis) UnmarshalBinary(data []byte) (err error) {
	_, err = p.ReadFrom(buffer.NewBuffer(data))
	return
}
