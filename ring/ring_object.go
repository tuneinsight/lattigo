package ring

import (
	"encoding/binary"
	"errors"
	"math/bits"
)

// Poly is the structure containing the coefficients of a polynomial.
type Poly struct {
	Coeffs [][]uint64 //Coefficients in CRT representation
}

// GetDegree returns the number of coefficients (degree) of the polynomial.
func (Pol *Poly) GetDegree() int {
	return len(Pol.Coeffs[0])
}

// GetLenModuli returns the number of modulies
func (Pol *Poly) GetLenModuli() int {
	return len(Pol.Coeffs)
}

// Zero sets all coefficient of the target polynomial to 0.
func (Pol *Poly) Zero() {
	for i := range Pol.Coeffs {
		for j := range Pol.Coeffs[0] {
			Pol.Coeffs[i][j] = 0
		}
	}
}

// CopyNew creates a new polynomial p1 which is a copy of the target polynomial.
func (Pol *Poly) CopyNew() (p1 *Poly) {
	p1 = new(Poly)
	p1.Coeffs = make([][]uint64, len(Pol.Coeffs))
	for i := range Pol.Coeffs {
		p1.Coeffs[i] = make([]uint64, len(Pol.Coeffs[i]))
		for j := range Pol.Coeffs[i] {
			p1.Coeffs[i][j] = Pol.Coeffs[i][j]
		}
	}

	return p1
}

// Copy copies the coefficients of p0 on p1 within the given context. Requiers p1 to be as big as the target context.
func (context *Context) Copy(p0, p1 *Poly) {

	if p0 != p1 {
		for i := range context.Modulus {
			for j := uint64(0); j < context.N; j++ {
				p1.Coeffs[i][j] = p0.Coeffs[i][j]
			}
		}
	}
}

// Copy copies the coefficients of Pol on p1.
func (Pol *Poly) Copy(p1 *Poly) {

	if Pol != p1 {
		for i := range p1.Coeffs {
			for j := range p1.Coeffs[i] {
				Pol.Coeffs[i][j] = p1.Coeffs[i][j]
			}
		}
	}
}

// SetCoefficients sets the coefficients of polynomial directly from a CRT format (double slice).
func (Pol *Poly) SetCoefficients(coeffs [][]uint64) error {

	if len(coeffs) > len(Pol.Coeffs) {
		return errors.New("error : len(coeffs) > len(Pol.Coeffs")
	}

	if len(coeffs[0]) > len(Pol.Coeffs[0]) {
		return errors.New("error : len(coeffs[0]) > len(Pol.Coeffs[0]")
	}

	for i := range coeffs {
		for j := range coeffs[0] {
			Pol.Coeffs[i][j] = coeffs[i][j]
		}
	}

	return nil
}

// GetCoefficients returns a double slice containing the coefficients of the polynomial.
func (Pol *Poly) GetCoefficients() [][]uint64 {
	coeffs := make([][]uint64, len(Pol.Coeffs))

	for i := range Pol.Coeffs {
		coeffs[i] = make([]uint64, len(Pol.Coeffs[i]))
		for j := range Pol.Coeffs[i] {
			coeffs[i][j] = Pol.Coeffs[i][j]
		}
	}

	return coeffs
}

// WriteCoeffsTo converts a matrix of coefficients to a byte array.
func WriteCoeffsTo(pointer, N, numberModuli uint64, coeffs [][]uint64, data []byte) (uint64, error) {
	tmp := N << 3
	for i := uint64(0); i < numberModuli; i++ {
		for j := uint64(0); j < N; j++ {
			binary.BigEndian.PutUint64(data[pointer+(j<<3):pointer+((j+1)<<3)], coeffs[i][j])
		}
		pointer += tmp
	}

	return pointer, nil
}

// DecodeCoeffs converts a byte array to a matrix of coefficients.
func DecodeCoeffs(pointer, N, numberModuli uint64, coeffs [][]uint64, data []byte) (uint64, error) {
	tmp := N << 3
	for i := uint64(0); i < numberModuli; i++ {
		for j := uint64(0); j < N; j++ {
			coeffs[i][j] = binary.BigEndian.Uint64(data[pointer+(j<<3) : pointer+((j+1)<<3)])
		}
		pointer += tmp
	}

	return pointer, nil
}

// DecodeCoeffs converts a byte array to a matrix of coefficients.
func DecodeCoeffsNew(pointer, N, numberModuli uint64, coeffs [][]uint64, data []byte) (uint64, error) {
	tmp := N << 3
	for i := uint64(0); i < numberModuli; i++ {
		coeffs[i] = make([]uint64, N)
		for j := uint64(0); j < N; j++ {
			coeffs[i][j] = binary.BigEndian.Uint64(data[pointer+(j<<3) : pointer+((j+1)<<3)])
		}
		pointer += tmp
	}

	return pointer, nil
}

func (Pol *Poly) MarshalBinary() ([]byte, error) {

	N := uint64(len(Pol.Coeffs[0]))
	numberModulies := uint64(len(Pol.Coeffs))

	data := make([]byte, 2+((N*numberModulies)<<3))

	if numberModulies > 0xFF {
		return nil, errors.New("error : poly max modulies uint16 overflow")
	}

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModulies)

	var pointer uint64

	pointer = 2

	if _, err := WriteCoeffsTo(pointer, N, numberModulies, Pol.Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

func (Pol *Poly) UnMarshalBinary(data []byte) (*Poly, error) {

	N := uint64(int(1 << data[0]))
	numberModulies := uint64(int(data[1]))

	var pointer uint64

	pointer = 2

	if ((uint64(len(data)) - pointer) >> 3) != N*numberModulies {
		return nil, errors.New("error : invalid polynomial encoding")
	}

	if _, err := DecodeCoeffs(pointer, N, numberModulies, Pol.Coeffs, data); err != nil {
		return nil, err
	}

	return Pol, nil
}
