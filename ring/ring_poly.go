package ring

import (
	"encoding/binary"
	"errors"
	"math/bits"
)

// Poly is the structure that contains the coefficients of a polynomial.
type Poly struct {
	Coeffs [][]uint64 // Coefficients in CRT representation
}

// NewPoly creates a new polynomial with N coefficients set to zero and nbModuli moduli.
func NewPoly(N, nbModuli int) (pol *Poly) {
	pol = new(Poly)
	pol.Coeffs = make([][]uint64, nbModuli)
	for i := 0; i < nbModuli; i++ {
		pol.Coeffs[i] = make([]uint64, N)
	}
	return
}

// Degree returns the number of coefficients of the polynomial, which equals the degree of the Ring cyclotomic polynomial.
func (pol *Poly) Degree() int {
	return len(pol.Coeffs[0])
}

// LenModuli returns the current number of moduli.
func (pol *Poly) LenModuli() int {
	return len(pol.Coeffs)
}

// Level returns the current number of moduli minus 1.
func (pol *Poly) Level() int {
	return len(pol.Coeffs) - 1
}

// Zero sets all coefficients of the target polynomial to 0.
func (pol *Poly) Zero() {
	for i := range pol.Coeffs {
		p0tmp := pol.Coeffs[i]
		for j := range p0tmp {
			p0tmp[j] = 0
		}
	}
}

// CopyNew creates a new polynomial p1 which is a copy of the target polynomial.
func (pol *Poly) CopyNew() (p1 *Poly) {
	p1 = new(Poly)
	p1.Coeffs = make([][]uint64, len(pol.Coeffs))
	for i := range pol.Coeffs {
		p1.Coeffs[i] = make([]uint64, len(pol.Coeffs[i]))
		p0tmp, p1tmp := pol.Coeffs[i], p1.Coeffs[i]
		for j := range pol.Coeffs[i] {
			p1tmp[j] = p0tmp[j]
		}
	}

	return p1
}

// Copy copies the coefficients of p0 on p1 within the given Ring. It requires p1 to be at least as big p0.
func (r *Ring) Copy(p0, p1 *Poly) {

	if p0 != p1 {
		for i := range r.Modulus {
			p0tmp, p1tmp := p0.Coeffs[i], p1.Coeffs[i]
			for j := 0; j < r.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		}
	}
}

// CopyLvl copies the coefficients of p0 on p1 within the given Ring for the moduli from 0 to level.
// Requires p1 to be as big as the target Ring.
func (r *Ring) CopyLvl(level int, p0, p1 *Poly) {

	if p0 != p1 {
		for i := 0; i < level+1; i++ {
			p0tmp, p1tmp := p0.Coeffs[i], p1.Coeffs[i]
			for j := 0; j < r.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		}
	}
}

// Copy copies the coefficients of p1 on the target polynomial.
func (pol *Poly) Copy(p1 *Poly) {

	if pol != p1 {
		for i := range p1.Coeffs {
			p0tmp, p1tmp := pol.Coeffs[i], p1.Coeffs[i]
			for j := range p1.Coeffs[i] {
				p0tmp[j] = p1tmp[j]
			}
		}
	}
}

// SetCoefficients sets the coefficients of the polynomial directly from a CRT format (double slice).
func (pol *Poly) SetCoefficients(coeffs [][]uint64) {
	for i := range coeffs {
		for j := range coeffs[i] {
			pol.Coeffs[i][j] = coeffs[i][j]
		}
	}
}

// GetCoefficients returns a new double slice that contains the coefficients of the polynomial.
func (pol *Poly) GetCoefficients() (coeffs [][]uint64) {
	coeffs = make([][]uint64, len(pol.Coeffs))

	for i := range pol.Coeffs {
		coeffs[i] = make([]uint64, len(pol.Coeffs[i]))
		for j := range pol.Coeffs[i] {
			coeffs[i][j] = pol.Coeffs[i][j]
		}
	}

	return
}

// WriteCoeffsTo converts a matrix of coefficients to a byte array.
func WriteCoeffsTo(pointer, N, numberModuli int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 3
	for i := 0; i < numberModuli; i++ {
		for j := 0; j < N; j++ {
			binary.BigEndian.PutUint64(data[pointer+(j<<3):pointer+((j+1)<<3)], coeffs[i][j])
		}
		pointer += tmp
	}

	return pointer, nil
}

// WriteTo writes the given poly to the data array.
// It returns the number of written bytes, and the corresponding error, if it occurred.
func (pol *Poly) WriteTo(data []byte) (int, error) {

	N := pol.Degree()
	numberModuli := pol.LenModuli()

	if len(data) < pol.GetDataLen(true) {
		// The data is not big enough to write all the information
		return 0, errors.New("data array is too small to write ring.Poly")
	}
	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModuli)

	cnt, err := WriteCoeffsTo(2, N, numberModuli, pol.Coeffs, data)

	return cnt, err
}

// WriteTo32 writes the given poly to the data array.
// It returns the number of written bytes, and the corresponding error, if it occurred.
func (pol *Poly) WriteTo32(data []byte) (int, error) {

	N := pol.Degree()
	numberModuli := pol.LenModuli()

	if len(data) < pol.GetDataLen32(true) {
		//The data is not big enough to write all the information
		return 0, errors.New("data array is too small to write ring.Poly")
	}
	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModuli)

	cnt, err := WriteCoeffsTo32(2, N, numberModuli, pol.Coeffs, data)

	return cnt, err
}

// WriteCoeffsTo32 converts a matrix of coefficients to a byte array.
func WriteCoeffsTo32(pointer, N, numberModuli int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 2
	for i := 0; i < numberModuli; i++ {
		for j := 0; j < N; j++ {
			binary.BigEndian.PutUint32(data[pointer+(j<<2):pointer+((j+1)<<2)], uint32(coeffs[i][j]))
		}
		pointer += tmp
	}

	return pointer, nil
}

// GetDataLen32 returns the number of bytes the polynomial will take when written to data.
// It can take into account meta data if necessary.
func (pol *Poly) GetDataLen32(WithMetadata bool) (cnt int) {
	cnt = (pol.LenModuli() * pol.Degree()) << 2

	if WithMetadata {
		cnt += 2
	}
	return
}

// WriteCoeffs writes the coefficients to the given data array.
// It fails if the data array is not big enough to contain the ring.Poly
func (pol *Poly) WriteCoeffs(data []byte) (int, error) {
	return WriteCoeffsTo(0, pol.Degree(), pol.LenModuli(), pol.Coeffs, data)

}

// GetDataLen returns the number of bytes the polynomial will take when written to data.
// It can take into account meta data if necessary.
func (pol *Poly) GetDataLen(WithMetadata bool) (cnt int) {
	cnt = (pol.LenModuli() * pol.Degree()) << 3

	if WithMetadata {
		cnt += 2
	}
	return
}

// DecodeCoeffs converts a byte array to a matrix of coefficients.
func DecodeCoeffs(pointer, N, numberModuli int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 3
	for i := 0; i < numberModuli; i++ {
		for j := 0; j < N; j++ {
			coeffs[i][j] = binary.BigEndian.Uint64(data[pointer+(j<<3) : pointer+((j+1)<<3)])
		}
		pointer += tmp
	}

	return pointer, nil
}

// DecodeCoeffsNew converts a byte array to a matrix of coefficients.
func DecodeCoeffsNew(pointer, N, numberModuli int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 3
	for i := 0; i < numberModuli; i++ {
		coeffs[i] = make([]uint64, N)
		for j := 0; j < N; j++ {
			coeffs[i][j] = binary.BigEndian.Uint64(data[pointer+(j<<3) : pointer+((j+1)<<3)])
		}
		pointer += tmp
	}

	return pointer, nil
}

// MarshalBinary encodes the target polynomial on a slice of bytes.
func (pol *Poly) MarshalBinary() (data []byte, err error) {
	data = make([]byte, pol.GetDataLen(true))
	_, err = pol.WriteTo(data)
	return
}

// UnmarshalBinary decodes a slice of byte on the target polynomial.
func (pol *Poly) UnmarshalBinary(data []byte) (err error) {

	N := uint64(1 << data[0])
	numberModulies := uint64(data[1])
	pointer := uint64(2)

	if ((uint64(len(data)) - pointer) >> 3) != N*numberModulies {
		return errors.New("invalid polynomial encoding")
	}

	if _, err = pol.DecodePolyNew(data); err != nil {
		return err
	}

	return nil
}

// DecodePolyNew decodes a slice of bytes in the target polynomial returns the number of bytes
// decoded.
func (pol *Poly) DecodePolyNew(data []byte) (pointer int, err error) {

	N := int(1 << data[0])
	numberModulies := int(data[1])
	pointer = 2

	if pol.Coeffs == nil {
		pol.Coeffs = make([][]uint64, numberModulies)
	}

	if pointer, err = DecodeCoeffsNew(pointer, N, numberModulies, pol.Coeffs, data); err != nil {
		return pointer, err
	}

	return pointer, nil
}

// DecodePolyNew32 decodes a slice of bytes in the target polynomial returns the number of bytes
// decoded.
func (pol *Poly) DecodePolyNew32(data []byte) (pointer int, err error) {

	N := int(1 << data[0])
	numberModulies := int(data[1])
	pointer = 2

	if pol.Coeffs == nil {
		pol.Coeffs = make([][]uint64, numberModulies)
	}

	if pointer, err = DecodeCoeffsNew32(pointer, N, numberModulies, pol.Coeffs, data); err != nil {
		return pointer, err
	}

	return pointer, nil
}

// DecodeCoeffsNew32 converts a byte array to a matrix of coefficients.
func DecodeCoeffsNew32(pointer, N, numberModuli int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 2
	for i := 0; i < numberModuli; i++ {
		coeffs[i] = make([]uint64, N)
		for j := 0; j < N; j++ {
			coeffs[i][j] = uint64(binary.BigEndian.Uint32(data[pointer+(j<<2) : pointer+((j+1)<<2)]))
		}
		pointer += tmp
	}

	return pointer, nil
}
