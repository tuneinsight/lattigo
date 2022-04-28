package ring

import (
	"encoding/binary"
	"errors"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Poly is the structure that contains the coefficients of a polynomial.
type Poly struct {
	Coeffs  [][]uint64 // Coefficients in CRT representation
	IsNTT   bool
	IsMForm bool
}

// NewPoly creates a new polynomial with N coefficients set to zero and Level+1 moduli.
func NewPoly(N, Level int) (pol *Poly) {
	pol = new(Poly)
	pol.Coeffs = make([][]uint64, Level+1)
	for i := 0; i < Level+1; i++ {
		pol.Coeffs[i] = make([]uint64, N)
	}
	return
}

// N returns the number of coefficients of the polynomial, which equals the degree of the Ring cyclotomic polynomial.
func (pol *Poly) N() int {
	return len(pol.Coeffs[0])
}

// Level returns the current number of moduli minus 1.
func (pol *Poly) Level() int {
	return len(pol.Coeffs) - 1
}

// Zero sets all coefficients of the target polynomial to 0.
func (pol *Poly) Zero() {
	for _, coeffs := range pol.Coeffs {
		ZeroVec(coeffs)
	}
}

// CopyNew creates an exact copy of the target polynomial.
func (pol *Poly) CopyNew() (p1 *Poly) {
	p1 = new(Poly)
	p1.Coeffs = make([][]uint64, pol.Level()+1)
	for i := range pol.Coeffs {
		p1.Coeffs[i] = make([]uint64, len(pol.Coeffs[i]))
		copy(p1.Coeffs[i], pol.Coeffs[i])
	}
	p1.IsNTT = pol.IsNTT
	p1.IsMForm = pol.IsMForm

	return
}

// CopyValues copies the coefficients of p0 on p1 within the given Ring. It requires p1 to be at least as big p0.
// Expects the degree of both polynomials to be identical.
// Does not transfer the IsNTT and IsMForm flags.
func CopyValues(p0, p1 *Poly) {
	CopyValuesLvl(utils.MinInt(p0.Level(), p1.Level()), p0, p1)

}

// Copy copies the coefficients of p0 on p1 within the given Ring. It requires p1 to be at least as big p0.
// Expects the degree of both polynomials to be identical.
// Transfers the IsNTT and IsMForm flags.
func Copy(p0, p1 *Poly) {
	CopyValuesLvl(utils.MinInt(p0.Level(), p1.Level()), p0, p1)
	p1.IsNTT = p0.IsNTT
	p1.IsMForm = p0.IsMForm
}

// CopyValuesLvl copies the coefficients of p0 on p1 within the given Ring for the moduli from 0 to level.
// Expects the degree of both polynomials to be identical.
// Does not transfer the IsNTT and IsMForm flags.
func CopyValuesLvl(level int, p0, p1 *Poly) {
	if p0 != p1 {
		for i := 0; i < level+1; i++ {
			copy(p1.Coeffs[i], p0.Coeffs[i])
		}
	}
}

// CopyLvl copies the coefficients of p0 on p1 within the given Ring for the moduli from 0 to level.
// Expects the degree of both polynomials to be identical.
// Transfers the IsNTT and IsMForm flags.
func CopyLvl(level int, p0, p1 *Poly) {
	CopyValuesLvl(level, p0, p1)
	p1.IsNTT = p0.IsNTT
	p1.IsMForm = p0.IsMForm
}

// CopyValues copies the coefficients of p1 on the target polynomial.
// Onyl copies minLevel(pol, p1) levels.
// Expects the degree of both polynomials to be identical.
// Does not transfer the IsNTT and IsMForm flags.
func (pol *Poly) CopyValues(p1 *Poly) {
	if pol != p1 {
		minLevel := utils.MinInt(pol.Level(), p1.Level())
		for i := range p1.Coeffs[:minLevel+1] {
			copy(pol.Coeffs[i], p1.Coeffs[i])
		}
	}
}

// Copy copies the coefficients of p1 on the target polynomial.
// Onyl copies minLevel(pol, p1) levels.
// Transfers the IsNTT and IsMForm flags.
func (pol *Poly) Copy(p1 *Poly) {
	pol.CopyValues(p1)
	pol.IsNTT = p1.IsNTT
	pol.IsMForm = p1.IsMForm
}

// Equals returns true if the receiver Poly is equal to the provided other Poly.
// This function checks for strict equality between the polynomial coefficients
// (i.e., it does not consider congruence as equality within the ring like
// `Ring.Equals` does).
// Will not check if IsNTT and IsMForm flags are equal
func (pol *Poly) Equals(other *Poly) bool {
	if pol == other {
		return true
	}
	if pol != nil && other != nil && len(pol.Coeffs) == len(other.Coeffs) {
		for i := range pol.Coeffs {
			if len(other.Coeffs[i]) != len(pol.Coeffs[i]) {
				return false
			}
			for j := range pol.Coeffs[i] {
				if other.Coeffs[i][j] != pol.Coeffs[i][j] {
					return false
				}
			}
		}
		return true
	}
	return false
}

// SetCoefficients sets the coefficients of the polynomial directly from a CRT format (double slice).
func (pol *Poly) SetCoefficients(coeffs [][]uint64) {
	for i := range coeffs {
		copy(pol.Coeffs[i], coeffs[i])
	}
}

// GetCoefficients returns a new double slice that contains the coefficients of the polynomial.
func (pol *Poly) GetCoefficients() (coeffs [][]uint64) {
	coeffs = make([][]uint64, len(pol.Coeffs))
	for i := range pol.Coeffs {
		coeffs[i] = make([]uint64, len(pol.Coeffs[i]))
		copy(coeffs[i], pol.Coeffs[i])
	}

	return
}

// GetDataLen64 returns the number of bytes a polynomial of N coefficients
// with Level+1 moduli will take when converted to a slice of bytes.
// Assumes that each coefficient will be encoded on 8 bytes.
// It can take into account meta data if necessary.
func GetDataLen64(N, Level int, WithMetadata bool) (cnt int) {
	cnt += N * (Level + 1) << 3

	if WithMetadata {
		cnt += 7
	}

	return
}

// GetDataLen64 returns the number of bytes the polynomial will take when written to data.
// Assumes that each coefficient takes 8 bytes.
// It can take into account meta data if necessary.
func (pol *Poly) GetDataLen64(WithMetadata bool) (cnt int) {
	return GetDataLen64(pol.N(), pol.Level(), WithMetadata)
}

// MarshalBinary encodes the target polynomial on a slice of bytes.
// Encodes each coefficient on 8 bytes.
func (pol *Poly) MarshalBinary() (data []byte, err error) {
	data = make([]byte, pol.GetDataLen64(true))
	_, err = pol.WriteTo64(data)
	return
}

// UnmarshalBinary decodes a slice of byte on the target polynomial.
// Assumes each coefficient is encoded on 8 bytes.
func (pol *Poly) UnmarshalBinary(data []byte) (err error) {

	N := int(binary.BigEndian.Uint32(data))
	Level := int(data[4])

	if data[5] == 1 {
		pol.IsNTT = true
	}

	if data[6] == 1 {
		pol.IsMForm = true
	}

	pointer := 6

	if ((len(data) - pointer) >> 3) != N*(Level+1) {
		return errors.New("invalid polynomial encoding")
	}

	if _, err = pol.DecodePoly64New(data); err != nil {
		return err
	}

	return nil
}

// WriteTo64 writes the given poly to the data array, using 8 bytes per coefficient.
// It returns the number of written bytes, and the corresponding error, if it occurred.
func (pol *Poly) WriteTo64(data []byte) (int, error) {

	N := pol.N()
	Level := pol.Level()

	if len(data) < pol.GetDataLen64(true) {
		// The data is not big enough to write all the information
		return 0, errors.New("data array is too small to write ring.Poly")
	}

	binary.BigEndian.PutUint32(data, uint32(N))

	data[4] = uint8(Level)

	if pol.IsNTT {
		data[5] = 1
	}

	if pol.IsMForm {
		data[6] = 1
	}

	return WriteCoeffsTo64(7, N, Level, pol.Coeffs, data)
}

// WriteCoeffsTo64 converts a matrix of coefficients to a byte array, using 8 bytes per coefficient.
func WriteCoeffsTo64(pointer, N, Level int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 3
	for i := 0; i < Level+1; i++ {
		for j, k := 0, pointer; j < N; j, k = j+1, k+8 {
			binary.BigEndian.PutUint64(data[k:], coeffs[i][j])
		}
		pointer += tmp
	}

	return pointer, nil
}

// DecodePoly64 decodes a slice of bytes in the target polynomial returns the number of bytes decoded.
// Assumes that each coefficient is encoded on 8 bytes.
func (pol *Poly) DecodePoly64(data []byte) (pointer int, err error) {

	N := int(binary.BigEndian.Uint32(data))
	Level := int(data[4])

	if data[5] == 1 {
		pol.IsNTT = true
	}

	if data[6] == 1 {
		pol.IsMForm = true
	}

	pointer = 7

	if pointer, err = DecodeCoeffs64(pointer, N, Level, pol.Coeffs, data); err != nil {
		return pointer, err
	}

	return pointer, nil
}

// DecodePoly64New decodes a slice of bytes in the target polynomial returns the number of bytes decoded.
// Assumes that each coefficient is encoded on 8 bytes.
func (pol *Poly) DecodePoly64New(data []byte) (pointer int, err error) {

	N := int(binary.BigEndian.Uint32(data))
	Level := int(data[4])

	if data[5] == 1 {
		pol.IsNTT = true
	}

	if data[6] == 1 {
		pol.IsMForm = true
	}

	pointer = 7

	if pol.Coeffs == nil {
		pol.Coeffs = make([][]uint64, Level+1)
	}

	if pointer, err = DecodeCoeffs64New(pointer, N, Level, pol.Coeffs, data); err != nil {
		return pointer, err
	}

	return pointer, nil
}

// DecodeCoeffs64 converts a byte array to a matrix of coefficients.
// Assumes that each coefficient is encoded on 8 bytes.
func DecodeCoeffs64(pointer, N, Level int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 3
	for i := 0; i < Level+1; i++ {
		for j, k := 0, pointer; j < N; j, k = j+1, k+8 {
			coeffs[i][j] = binary.BigEndian.Uint64(data[k:])
		}
		pointer += tmp
	}

	return pointer, nil
}

// DecodeCoeffs64New converts a byte array to a matrix of coefficients.
// Assumes taht each coefficient is encoded on 8 bytes.
func DecodeCoeffs64New(pointer, N, Level int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 3
	for i := 0; i < Level+1; i++ {
		coeffs[i] = make([]uint64, N)
		for j, k := 0, pointer; j < N; j, k = j+1, k+8 {
			coeffs[i][j] = binary.BigEndian.Uint64(data[k:])
		}
		pointer += tmp
	}

	return pointer, nil
}

// WriteTo32 writes the given poly to the data array.
// Encodes each coefficient on 4 bytes.
// It returns the number of written bytes, and the corresponding error, if it occurred.
func (pol *Poly) WriteTo32(data []byte) (int, error) {

	N := pol.N()
	Level := pol.Level()

	if len(data) < pol.GetDataLen32(true) {
		//The data is not big enough to write all the information
		return 0, errors.New("data array is too small to write ring.Poly")
	}
	binary.BigEndian.PutUint32(data, uint32(N))
	data[4] = uint8(Level)
	if pol.IsNTT {
		data[5] = 1
	}

	if pol.IsMForm {
		data[6] = 1
	}

	cnt, err := WriteCoeffsTo32(7, N, Level, pol.Coeffs, data)

	return cnt, err
}

// WriteCoeffsTo32 converts a matrix of coefficients to a byte array.
// Encodes each coefficient on 4 bytes.
func WriteCoeffsTo32(pointer, N, Level int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 2
	for i := 0; i < Level+1; i++ {
		for j := 0; j < N; j++ {
			binary.BigEndian.PutUint32(data[pointer+(j<<2):pointer+((j+1)<<2)], uint32(coeffs[i][j]))
		}
		pointer += tmp
	}

	return pointer, nil
}

// GetPolyDataLen32 returns the number of bytes a polynomial of N coefficients
// with Level+1 moduli will take when converted to a slice of bytes.
// Assumes that each coefficient will be encoded on 4 bytes.
// It can take into account meta data if necessary.
func GetPolyDataLen32(N, Level int, WithMetadata bool) (cnt int) {
	cnt += N * (Level + 1) << 4

	if WithMetadata {
		cnt += 7
	}

	return
}

// GetDataLen32 returns the number of bytes the polynomial will take when written to data.
// Assumes that each coefficient is encoded on 4 bytes.
// It can take into account meta data if necessary.
func (pol *Poly) GetDataLen32(WithMetadata bool) (cnt int) {
	return GetPolyDataLen32(pol.N(), pol.Level(), WithMetadata)
}

// DecodePoly32 decodes a slice of bytes in the target polynomial returns the number of bytes decoded.
// Assumes that each coefficient is encoded on 4 bytes.
// Writes on pre-allocated coefficients.
func (pol *Poly) DecodePoly32(data []byte) (pointer int, err error) {

	N := int(binary.BigEndian.Uint32(data))
	Level := int(data[4])

	if data[5] == 1 {
		pol.IsNTT = true
	}

	if data[6] == 1 {
		pol.IsMForm = true
	}

	pointer = 7

	if pointer, err = DecodeCoeffs32(pointer, N, Level, pol.Coeffs, data); err != nil {
		return pointer, err
	}

	return pointer, nil
}

// DecodePoly32New decodes a slice of bytes in the target polynomial returns the number of bytes decoded.
// Assumes that each coefficient is encoded on 4 bytes.
// Allocates the coefficients.
func (pol *Poly) DecodePoly32New(data []byte) (pointer int, err error) {

	N := int(binary.BigEndian.Uint32(data))
	Level := int(data[4])

	if data[5] == 1 {
		pol.IsNTT = true
	}

	if data[6] == 1 {
		pol.IsMForm = true
	}

	pointer = 7

	if pol.Coeffs == nil {
		pol.Coeffs = make([][]uint64, Level+1)
	}

	if pointer, err = DecodeCoeffs32New(pointer, N, Level, pol.Coeffs, data); err != nil {
		return pointer, err
	}

	return pointer, nil
}

// DecodeCoeffs32 converts a byte array to a matrix of coefficients.
// Assumes that each coefficient is encoded on 4 bytes.
func DecodeCoeffs32(pointer, N, Level int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 2
	for i := 0; i < Level+1; i++ {
		for j, k := 0, pointer; j < N; j, k = j+1, k+4 {
			coeffs[i][j] = uint64(binary.BigEndian.Uint32(data[k:]))
		}
		pointer += tmp
	}

	return pointer, nil
}

// DecodeCoeffs32New converts a byte array to a matrix of coefficients.
// Assumes that each coefficient is encoded on 4 bytes.
func DecodeCoeffs32New(pointer, N, Level int, coeffs [][]uint64, data []byte) (int, error) {
	tmp := N << 2
	for i := 0; i < Level+1; i++ {
		coeffs[i] = make([]uint64, N)
		for j, k := 0, pointer; j < N; j, k = j+1, k+4 {
			coeffs[i][j] = uint64(binary.BigEndian.Uint32(data[k:]))
		}
		pointer += tmp
	}

	return pointer, nil
}
