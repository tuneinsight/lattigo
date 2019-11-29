package ring

import (
	"crypto/rand"
	"encoding/binary"
	"errors"
	"math/bits"
)

// Poly is the structure containing the coefficients of a polynomial.
type Poly struct {
	Coeffs [][]uint64 //Coefficients in CRT representation
}

func NewPoly(N, nbModuli uint64) (pol *Poly) {
	pol = new(Poly)
	pol.Coeffs = make([][]uint64, nbModuli)
	for i := uint64(0); i < nbModuli; i++ {
		pol.Coeffs[i] = make([]uint64, N)
	}
	return
}

func NewPolyUniform(N, nbModuli uint64) (pol *Poly) {
	pol = new(Poly)

	randomBytes := make([]byte, N<<3)

	pol.Coeffs = make([][]uint64, nbModuli)
	for i := uint64(0); i < nbModuli; i++ {
		pol.Coeffs[i] = make([]uint64, N)

		tmp := pol.Coeffs[i]

		if _, err := rand.Read(randomBytes); err != nil {
			panic("crypto rand error")
		}

		for j := uint64(0); j < N; j++ {
			tmp[j] = binary.BigEndian.Uint64(randomBytes[j<<3 : (j+1)<<3])
		}
	}

	return
}

// GetDegree returns the number of coefficients (degree) of the polynomial.
func (pol *Poly) GetDegree() int {
	return len(pol.Coeffs[0])
}

// GetLenModuli returns the number of modulie.
func (pol *Poly) GetLenModuli() int {
	return len(pol.Coeffs)
}

// Zero sets all coefficients of the target polynomial to 0.
func (pol *Poly) Zero() {
	for i := range pol.Coeffs {
		p0tmp := pol.Coeffs[i]
		for j := range pol.Coeffs[0] {
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

// Copy copies the coefficients of p0 on p1 within the given context. Requires p1 to be as big as the target context.
func (context *Context) Copy(p0, p1 *Poly) {

	if p0 != p1 {
		for i := range context.Modulus {
			p0tmp, p1tmp := p0.Coeffs[i], p1.Coeffs[i]
			for j := uint64(0); j < context.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		}
	}
}

// CopyLvl copies the coefficients of p0 on p1 within the given context. Requiers p1 to be as big as the target context.
func (context *Context) CopyLvl(level uint64, p0, p1 *Poly) {

	if p0 != p1 {
		for i := uint64(0); i < level+1; i++ {
			p0tmp, p1tmp := p0.Coeffs[i], p1.Coeffs[i]
			for j := uint64(0); j < context.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		}
	}
}

// Copy copies the receiver's coefficients from p1.
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

// SetCoefficients sets the coefficients of polynomial directly from a CRT format (double slice).
func (pol *Poly) SetCoefficients(coeffs [][]uint64) {
	for i := range coeffs {
		for j := range coeffs[0] {
			pol.Coeffs[i][j] = coeffs[i][j]
		}
	}
}

// GetCoefficients returns a double slice containing the coefficients of the polynomial.
func (pol *Poly) GetCoefficients() [][]uint64 {
	coeffs := make([][]uint64, len(pol.Coeffs))

	for i := range pol.Coeffs {
		coeffs[i] = make([]uint64, len(pol.Coeffs[i]))
		for j := range pol.Coeffs[i] {
			coeffs[i][j] = pol.Coeffs[i][j]
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

// WriteTo writes the given poly to the data array
// returns the number of bytes written and error if it occured.
func (pol *Poly) WriteTo(data []byte) (uint64, error) {

	N := uint64(pol.GetDegree())
	numberModulies := uint64(pol.GetLenModuli())

	if uint64(len(data)) < pol.GetDataLen(true) {
		//the data is not big enough to write all the information
		return 0, errors.New("Data array is too small to write ring.Poly")
	}
	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModulies)

	cnt, err := WriteCoeffsTo(2, N, numberModulies, pol.Coeffs, data)

	return cnt, err
}

// WriteCoeffs write the coefficient to the given data array.
// Fails if the data array is not big enough to contain the ring.Poly
func (pol *Poly) WriteCoeffs(data []byte) (uint64, error) {

	cnt, err := WriteCoeffsTo(0, uint64(pol.GetDegree()), uint64(pol.GetLenModuli()), pol.Coeffs, data)
	return cnt, err

}

// GetDataLen returns the lenght the poly will take when written to data.
// Can take into account meta data if necessary.
func (pol *Poly) GetDataLen(WithMetadata bool) uint64 {
	cnt := 0
	if WithMetadata {
		cnt = 2
	}
	return uint64(cnt + (pol.GetLenModuli()*pol.GetDegree())<<3)
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

// DecodeCoeffsNew converts a byte array to a matrix of coefficients.
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

// MarshalBinary encodes the target polynomial on a slice of bytes.
func (pol *Poly) MarshalBinary() ([]byte, error) {

	//N := uint64(len(pol.Coeffs[0]))
	//numberModulies := uint64(len(pol.Coeffs))
	data := make([]byte, pol.GetDataLen(true))

	_, err := pol.WriteTo(data)
	return data, err
	//if numberModulies > 0xFF {
	//	return nil, errors.New("error : poly max modulies uint16 overflow")
	//}
	//
	//data[0] = uint8(bits.Len64(uint64(N)) - 1)
	//data[1] = uint8(numberModulies)
	//
	//var pointer uint64
	//
	//pointer = 2
	////test to check if our method works.
	//
	//if _, err := WriteCoeffsTo(pointer, N, numberModulies, pol.Coeffs, data); err != nil {
	//	return nil, err
	//}

	//return data, nil
}

// UnmarshalBinary decodes a slice of byte on the target polynomial.
func (pol *Poly) UnmarshalBinary(data []byte) (err error) {

	N := uint64(1 << data[0])
	numberModulies := uint64(data[1])

	var pointer uint64

	pointer = 2

	if ((uint64(len(data)) - pointer) >> 3) != N*numberModulies {
		return errors.New("error : invalid polynomial encoding")
	}
	if _, err = pol.DecodePolyNew(data); err != nil {
		return err
	}

	return nil
}

// DecodePolyNew decodes a slice of bytes in the target polynomial returns the number of bytes
// decoded.
func (pol *Poly) DecodePolyNew(data []byte) (pointer uint64, err error) {

	N := uint64(1 << data[0])
	numberModulies := uint64(data[1])

	pointer = 2

	if pol.Coeffs == nil {
		pol.Coeffs = make([][]uint64, numberModulies)
	}
	if pointer, err = DecodeCoeffsNew(pointer, N, numberModulies, pol.Coeffs, data); err != nil {
		return pointer, err
	}

	return pointer, nil
}
