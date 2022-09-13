package rlwe

import (
	"encoding/binary"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v3/ring"
)

type EvaluatorInterface interface {
	Add(op0, op1, op2 *Ciphertext)
	Sub(op0, op1, op2 *Ciphertext)
	AddConst(op0 *Ciphertext, constant interface{}, op1 *Ciphertext)
	Rescale(op0, op1 *Ciphertext) (err error)
	Mul(op0, op1, op2 *Ciphertext)
	MulNew(op0, op1 *Ciphertext) (op2 *Ciphertext)
	MulAndAdd(op0, op1, op2 *Ciphertext)
	MulRelinNew(op0, op1 *Ciphertext) (op2 *Ciphertext)
	Relinearize(op0 *Ciphertext, op1 *Ciphertext)
	Parameters() Parameters
}

// PolynomialBasis is a struct storing powers of a ciphertext.
type PolynomialBasis struct {
	PolynomialBasisType
	Value map[int]*Ciphertext
}

// NewPolynomialBasis creates a new PolynomialBasis. It takes as input a ciphertext
// and a basistype. The struct treates the input ciphertext as a monomial X and
// can be used to generates power of this monomial X^{n} in the given BasisType.
func NewPolynomialBasis(ct *Ciphertext, basistype PolynomialBasisType) (p *PolynomialBasis) {
	p = new(PolynomialBasis)
	p.Value = make(map[int]*Ciphertext)
	p.Value[1] = ct.CopyNew()
	p.PolynomialBasisType = basistype
	return
}

// GenPower recursively computes X^{n}.
// If lazy = true, the final X^{n} will not be relinearized.
// Previous non-relinearized X^{n} that are required to compute the target X^{n} are automatically relinearized.
// Scale sets the threshold for rescaling (ciphertext won't be rescaled if the rescaling operation would make the scale go under this threshold).
func (p *PolynomialBasis) GenPower(n int, lazy bool, eval EvaluatorInterface) (err error) {

	if p.Value[n] == nil {
		if err = p.genPower(n, lazy, eval); err != nil {
			return
		}

		if err = eval.Rescale(p.Value[n], p.Value[n]); err != nil {
			return
		}
	}

	return nil
}

func (p *PolynomialBasis) GenPowersFromPolynomial(pol Polynomial, eval EvaluatorInterface) (err error) {

	if levels, depth := p.Value[1].Level(), pol.Depth(); levels < depth {
		return fmt.Errorf("%d levels < %d log(d) -> cannot evaluate", levels, depth)
	}

	giant, baby := pol.BSGSSplit()

	odd, even := pol.OddEven()

	isRingStandard := (eval.Parameters().RingType() == ring.Standard) && !pol.IsEncrypted()

	for i := (1 << baby) - 1; i > 1; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = p.GenPower(i, isRingStandard, eval); err != nil {
				return
			}
		}
	}

	for i := baby; i < giant; i++ {
		if err = p.GenPower(1<<i, false, eval); err != nil {
			return
		}
	}

	return
}

func (p *PolynomialBasis) genPower(n int, lazy bool, eval EvaluatorInterface) (err error) {
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

			if p.PolynomialBasisType == Chebyshev {
				c = int(math.Abs(float64(a) - float64(b))) // Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
			}
		}

		// Recurses on the given indexes
		if err = p.genPower(a, lazy && !isPow2, eval); err != nil {
			return err
		}
		if err = p.genPower(b, lazy && !isPow2, eval); err != nil {
			return err
		}

		if p.Value[a].Degree() == 2 {
			eval.Relinearize(p.Value[a], p.Value[a])
		}

		if p.Value[b].Degree() == 2 {
			eval.Relinearize(p.Value[b], p.Value[b])
		}

		if err = eval.Rescale(p.Value[a], p.Value[a]); err != nil {
			return err
		}

		if err = eval.Rescale(p.Value[b], p.Value[b]); err != nil {
			return err
		}

		// Computes C[n] = C[a]*C[b]
		if lazy {
			p.Value[n] = eval.MulNew(p.Value[a], p.Value[b])

		} else {
			p.Value[n] = eval.MulRelinNew(p.Value[a], p.Value[b])
			if err = eval.Rescale(p.Value[n], p.Value[n]); err != nil {
				return err
			}
		}

		if p.PolynomialBasisType == Chebyshev {

			// Computes C[n] = 2*C[a]*C[b]
			eval.Add(p.Value[n], p.Value[n], p.Value[n])

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				eval.AddConst(p.Value[n], -1, p.Value[n])
			} else {
				// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
				if err = p.GenPower(c, lazy, eval); err != nil {
					return err
				}
				eval.Sub(p.Value[n], p.Value[c], p.Value[n])
			}
		}
	}
	return
}

// MarshalBinary encodes the target on a slice of bytes.
func (p *PolynomialBasis) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 16)
	binary.LittleEndian.PutUint64(data[0:8], uint64(len(p.Value)))
	binary.LittleEndian.PutUint64(data[8:16], uint64(p.Value[1].GetDataLen(true)))
	for key, ct := range p.Value {
		keyBytes := make([]byte, 8)
		binary.LittleEndian.PutUint64(keyBytes, uint64(key))
		data = append(data, keyBytes...)
		ctBytes, err := ct.MarshalBinary()
		if err != nil {
			return []byte{}, err
		}
		data = append(data, ctBytes...)
	}
	return
}

// UnmarshalBinary decodes a slice of bytes on the target.
func (p *PolynomialBasis) UnmarshalBinary(data []byte) (err error) {
	p.Value = make(map[int]*Ciphertext)
	nbct := int(binary.LittleEndian.Uint64(data[0:8]))
	dtLen := int(binary.LittleEndian.Uint64(data[8:16]))
	ptr := 16
	for i := 0; i < nbct; i++ {
		idx := int(binary.LittleEndian.Uint64(data[ptr : ptr+8]))
		ptr += 8
		p.Value[idx] = new(Ciphertext)
		if err = p.Value[idx].UnmarshalBinary(data[ptr : ptr+dtLen]); err != nil {
			return
		}
		ptr += dtLen
	}
	return
}
