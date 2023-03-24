package ckks

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"math"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
)

// PolynomialBasis is a struct storing powers of a ciphertext.
type PolynomialBasis struct {
	BasisType
	Value map[int]*rlwe.Ciphertext
}

// NewPolynomialBasis creates a new PolynomialBasis. It takes as input a ciphertext
// and a basistype. The struct treats the input ciphertext as a monomial X and
// can be used to generates power of this monomial X^{n} in the given BasisType.
func NewPolynomialBasis(ct *rlwe.Ciphertext, basistype BasisType) (p *PolynomialBasis) {
	p = new(PolynomialBasis)
	p.Value = make(map[int]*rlwe.Ciphertext)
	p.Value[1] = ct.CopyNew()
	p.BasisType = basistype
	return
}

// GenPower recursively computes X^{n}.
// If lazy = true, the final X^{n} will not be relinearized.
// Previous non-relinearized X^{n} that are required to compute the target X^{n} are automatically relinearized.
// Scale sets the threshold for rescaling (ciphertext won't be rescaled if the rescaling operation would make the scale go under this threshold).
func (p *PolynomialBasis) GenPower(n int, lazy bool, scale rlwe.Scale, eval Evaluator) (err error) {

	if p.Value[n] == nil {
		if err = p.genPower(n, lazy, scale, eval); err != nil {
			return
		}

		if err = eval.Rescale(p.Value[n], scale, p.Value[n]); err != nil {
			return
		}
	}

	return nil
}

func (p *PolynomialBasis) genPower(n int, lazy bool, scale rlwe.Scale, eval Evaluator) (err error) {

	if p.Value[n] == nil {

		isPow2 := n&(n-1) == 0

		// Computes the index required to compute the asked ring evaluation
		var a, b, c int
		if isPow2 {
			a, b = n/2, n/2 //Necessary for optimal depth
		} else {
			// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
			// Maximize the number of odd terms of Chebyshev basis
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)

			if p.BasisType == Chebyshev {
				c = int(math.Abs(float64(a) - float64(b))) // Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
			}
		}

		// Recurses on the given indexes
		if err = p.genPower(a, lazy && !isPow2, scale, eval); err != nil {
			return err
		}
		if err = p.genPower(b, lazy && !isPow2, scale, eval); err != nil {
			return err
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
				return err
			}

			if err = eval.Rescale(p.Value[b], scale, p.Value[b]); err != nil {
				return err
			}

			p.Value[n] = eval.MulNew(p.Value[a], p.Value[b])

		} else {

			if err = eval.Rescale(p.Value[a], scale, p.Value[a]); err != nil {
				return err
			}

			if err = eval.Rescale(p.Value[b], scale, p.Value[b]); err != nil {
				return err
			}

			p.Value[n] = eval.MulRelinNew(p.Value[a], p.Value[b])
		}

		if p.BasisType == Chebyshev {

			// Computes C[n] = 2*C[a]*C[b]
			eval.Add(p.Value[n], p.Value[n], p.Value[n])

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				eval.AddConst(p.Value[n], -1, p.Value[n])
			} else {
				// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
				if err = p.GenPower(c, lazy, scale, eval); err != nil {
					return err
				}
				eval.Sub(p.Value[n], p.Value[c], p.Value[n])
			}
		}
	}
	return
}

func (p *PolynomialBasis) BinarySize() (size int) {
	size = 5 // Type & #Ct
	for _, ct := range p.Value {
		size += 4 + ct.BinarySize()
	}

	return
}

// MarshalBinary encodes the target on a slice of bytes.
func (p *PolynomialBasis) MarshalBinary() (data []byte, err error) {
	data = make([]byte, p.BinarySize())
	_, err = p.Read(data)
	return
}

func (p *PolynomialBasis) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		var inc1 int

		if inc1, err = buffer.WriteUint8(w, uint8(p.BasisType)); err != nil {
			return n + int64(inc1), err
		}

		n += int64(inc1)

		if inc1, err = buffer.WriteUint32(w, uint32(len(p.Value))); err != nil {
			return n + int64(inc1), err
		}

		n += int64(inc1)

		for _, key := range utils.GetSortedKeys(p.Value) {

			ct := p.Value[key]

			if inc1, err = buffer.WriteUint32(w, uint32(key)); err != nil {
				return n + int64(inc1), err
			}

			n += int64(inc1)

			var inc2 int64
			if inc2, err = ct.WriteTo(w); err != nil {
				return n + inc2, err
			}

			n += inc2
		}

		return

	default:
		return p.WriteTo(bufio.NewWriter(w))
	}
}

func (p *PolynomialBasis) Read(data []byte) (n int, err error) {

	if len(data) < p.BinarySize() {
		return n, fmt.Errorf("cannot Read: len(data)=%d < %d", len(data), p.BinarySize())
	}

	data[n] = uint8(p.BasisType)
	n++

	binary.LittleEndian.PutUint32(data[n:], uint32(len(p.Value)))
	n += 4

	for _, key := range utils.GetSortedKeys(p.Value) {

		ct := p.Value[key]

		binary.LittleEndian.PutUint32(data[n:], uint32(key))
		n += 4

		var inc int
		if inc, err = ct.Read(data[n:]); err != nil {
			return n + inc, err
		}

		n += inc
	}

	return
}

// UnmarshalBinary decodes a slice of bytes on the target.
func (p *PolynomialBasis) UnmarshalBinary(data []byte) (err error) {
	_, err = p.Write(data)
	return
}

func (p *PolynomialBasis) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:
		var inc1 int

		var BType uint8

		if inc1, err = buffer.ReadUint8(r, &BType); err != nil {
			return n + int64(inc1), err
		}

		n += int64(inc1)

		p.BasisType = BasisType(BType)

		var nbCts uint32
		if inc1, err = buffer.ReadUint32(r, &nbCts); err != nil {
			return n + int64(inc1), err
		}

		n += int64(inc1)

		p.Value = make(map[int]*rlwe.Ciphertext)

		for i := 0; i < int(nbCts); i++ {

			var key uint32

			if inc1, err = buffer.ReadUint32(r, &key); err != nil {
				return n + int64(inc1), err
			}

			n += int64(inc1)

			if p.Value[int(key)] == nil {
				p.Value[int(key)] = new(rlwe.Ciphertext)
			}

			var inc2 int64
			if inc2, err = p.Value[int(key)].ReadFrom(r); err != nil {
				return n + inc2, err
			}

			n += inc2
		}

		return

	default:
		return p.ReadFrom(bufio.NewReader(r))
	}
}

func (p *PolynomialBasis) Write(data []byte) (n int, err error) {

	p.BasisType = BasisType(data[n])
	n++

	nbCts := int(binary.LittleEndian.Uint32(data[n:]))
	n += 4

	p.Value = make(map[int]*rlwe.Ciphertext)

	for i := 0; i < nbCts; i++ {

		idx := int(binary.LittleEndian.Uint32(data[n:]))
		n += 4

		if p.Value[idx] == nil {
			p.Value[idx] = new(rlwe.Ciphertext)
		}

		var inc int
		if inc, err = p.Value[idx].Write(data[n:]); err != nil {
			return n + inc, err
		}

		n += inc
	}

	return
}
