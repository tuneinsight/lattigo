package ring

import (
	"encoding/binary"
	"errors"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// Poly is the structure that contains the coefficients of a polynomial.
type Poly struct {
	Coeffs [][]uint64 // Dimension-2 slice of coefficients (re-slice of Buff)
	Buff   []uint64   // Dimension-1 slice of coefficient
}

// NewPoly creates a new polynomial with N coefficients set to zero and Level+1 moduli.
func NewPoly(N, Level int) (pol *Poly) {
	pol = new(Poly)

	pol.Buff = make([]uint64, N*(Level+1))
	pol.Coeffs = make([][]uint64, Level+1)
	for i := 0; i < Level+1; i++ {
		pol.Coeffs[i] = pol.Buff[i*N : (i+1)*N]
	}

	return
}

// Resize resizes the level of the target polynomial to the provided level.
// If the provided level is larger than the current level, then allocates zero
// coefficients, otherwise dereferences the coefficients above the provided level.
func (pol *Poly) Resize(level int) {
	N := pol.N()
	if pol.Level() > level {
		pol.Buff = pol.Buff[:N*(level+1)]
		pol.Coeffs = pol.Coeffs[:level+1]
	} else if level > pol.Level() {
		pol.Buff = append(pol.Buff, make([]uint64, N*(level-pol.Level()))...)
		pol.Coeffs = append(pol.Coeffs, make([][]uint64, level-pol.Level())...)
		for i := 0; i < level+1; i++ {
			pol.Coeffs[i] = pol.Buff[i*N : (i+1)*N]
		}
	}
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
	ZeroVec(pol.Buff)
}

// CopyNew creates an exact copy of the target polynomial.
func (pol *Poly) CopyNew() (p1 *Poly) {
	p1 = NewPoly(pol.N(), pol.Level())
	copy(p1.Buff, pol.Buff)
	return
}

// Copy copies the coefficients of p0 on p1 within the given Ring. It requires p1 to be at least as big p0.
// Expects the degree of both polynomials to be identical.
func Copy(p0, p1 *Poly) {
	copy(p1.Buff, p0.Buff)
}

// CopyLvl copies the coefficients of p0 on p1 within the given Ring.
// Copies for up to level+1 moduli.
// Expects the degree of both polynomials to be identical.
func CopyLvl(level int, p0, p1 *Poly) {
	copy(p1.Buff[:p1.N()*(level+1)], p0.Buff)
}

// CopyValues copies the coefficients of p1 on the target polynomial.
// Onyl copies minLevel(pol, p1) levels.
// Expects the degree of both polynomials to be identical.
func (pol *Poly) CopyValues(p1 *Poly) {
	if pol != p1 {
		copy(pol.Buff, p1.Buff)
	}
}

// Copy copies the coefficients of p1 on the target polynomial.
// Onyl copies minLevel(pol, p1) levels.
func (pol *Poly) Copy(p1 *Poly) {
	pol.CopyValues(p1)
}

// Equals returns true if the receiver Poly is equal to the provided other Poly.
// This function checks for strict equality between the polynomial coefficients
// (i.e., it does not consider congruence as equality within the ring like
// `Ring.Equals` does).
func (pol *Poly) Equals(other *Poly) bool {
	if pol == other {
		return true
	}

	if pol != nil && other != nil && len(pol.Buff) == len(other.Buff) {
		for i := range pol.Buff {
			if other.Buff[i] != pol.Buff[i] {
				return false
			}
		}
		return true
	}
	return false
}

// MarshalBinarySize64 returns the number of bytes a polynomial of N coefficients
// with Level+1 moduli will take when converted to a slice of bytes.
// Assumes that each coefficient will be encoded on 8 bytes.
func MarshalBinarySize64(N, Level int) (cnt int) {
	return 5 + N*(Level+1)<<3
}

// MarshalBinarySize64 returns the number of bytes the polynomial will take when written to data.
// Assumes that each coefficient takes 8 bytes.
func (pol *Poly) MarshalBinarySize64() (cnt int) {
	return MarshalBinarySize64(pol.N(), pol.Level())
}

// MarshalBinary encodes the target polynomial on a slice of bytes.
// Encodes each coefficient on 8 bytes.
func (pol *Poly) MarshalBinary() (data []byte, err error) {
	data = make([]byte, pol.MarshalBinarySize64())
	_, err = pol.Encode64(data)
	return
}

// UnmarshalBinary decodes a slice of byte on the target polynomial.
// Assumes each coefficient is encoded on 8 bytes.
func (pol *Poly) UnmarshalBinary(data []byte) (err error) {

	N := int(binary.BigEndian.Uint32(data))
	Level := int(data[4])

	ptr := 5

	if ((len(data) - ptr) >> 3) != N*(Level+1) {
		return errors.New("invalid polynomial encoding")
	}

	if _, err = pol.Decode64(data); err != nil {
		return err
	}

	return nil
}

func (pol *Poly) Write(w *utils.Writer) (n int, err error) {

	var inc int
	if inc, err = w.WriteUint32(uint32(pol.N())); err != nil {
		return n + inc, err
	}

	n += inc

	if inc, err = w.WriteUint8(uint8(pol.Level())); err != nil {
		return n + inc, err
	}

	n += inc

	if inc, err = w.WriteUint64Slice(pol.Buff); err != nil {
		return n + inc, err
	}

	n += inc

	return n, w.Flush()
}

// Encode64 writes the given poly to the data array, using 8 bytes per coefficient.
// It returns the number of written bytes, and the corresponding error, if it occurred.
func (pol *Poly) Encode64(data []byte) (int, error) {

	N := pol.N()
	Level := pol.Level()

	if len(data) < pol.MarshalBinarySize64() {
		// The data is not big enough to write all the information
		return 0, errors.New("data array is too small to write ring.Poly")
	}

	binary.BigEndian.PutUint32(data, uint32(N))

	data[4] = uint8(Level)

	return Encode64(5, pol.Buff, data)
}

// Encode64 converts a matrix of coefficients to a byte array, using 8 bytes per coefficient.
func Encode64(ptr int, coeffs []uint64, data []byte) (int, error) {
	for i, j := 0, ptr; i < len(coeffs); i, j = i+1, j+8 {
		binary.BigEndian.PutUint64(data[j:], coeffs[i])
	}

	return ptr + len(coeffs)*8, nil
}

// Decode64 decodes a slice of bytes in the target polynomial and returns the number of bytes decoded.
// The method will first try to write on the buffer. If this step fails, either because the buffer isn't
// allocated or because it is of the wrong size, the method will allocate the correct buffer.
// Assumes that each coefficient is encoded on 8 bytes.
func (pol *Poly) Decode64(data []byte) (ptr int, err error) {

	N := int(binary.BigEndian.Uint32(data))
	Level := int(data[4])

	ptr = 5

	if pol.Buff == nil || len(pol.Buff) != N*(Level+1) {
		pol.Buff = make([]uint64, N*(Level+1))
	}

	if ptr, err = Decode64(ptr, pol.Buff, data); err != nil {
		return ptr, err
	}

	// Reslice
	pol.Coeffs = make([][]uint64, Level+1)
	for i := 0; i < Level+1; i++ {
		pol.Coeffs[i] = pol.Buff[i*N : (i+1)*N]
	}

	return ptr, nil
}

// Decode64 converts a byte array to a matrix of coefficients.
// Assumes that each coefficient is encoded on 8 bytes.
func Decode64(ptr int, coeffs []uint64, data []byte) (int, error) {
	for i, j := 0, ptr; i < len(coeffs); i, j = i+1, j+8 {
		coeffs[i] = binary.BigEndian.Uint64(data[j:])
	}

	return ptr + len(coeffs)*8, nil
}

// Encode32 writes the given poly to the data array.
// Encodes each coefficient on 4 bytes.
// It returns the number of written bytes, and the corresponding error, if it occurred.
func (pol *Poly) Encode32(data []byte) (int, error) {

	N := pol.N()
	Level := pol.Level()

	if len(data) < pol.MarshalBinarySize32() {
		//The data is not big enough to write all the information
		return 0, errors.New("data array is too small to write ring.Poly")
	}

	binary.BigEndian.PutUint32(data, uint32(N))
	data[4] = uint8(Level)

	return Encode32(5, pol.Buff, data)
}

// Encode32 converts a matrix of coefficients to a byte array, using 4 bytes per coefficient.
func Encode32(ptr int, coeffs []uint64, data []byte) (int, error) {
	for i, j := 0, ptr; i < len(coeffs); i, j = i+1, j+4 {
		binary.BigEndian.PutUint32(data[j:], uint32(coeffs[i]))
	}

	return ptr + len(coeffs)*4, nil
}

// MarshalBinarySize32 returns the number of bytes a polynomial of N coefficients
// with Level+1 moduli will take when converted to a slice of bytes.
// Assumes that each coefficient will be encoded on 4 bytes.
func MarshalBinarySize32(N, Level int) (cnt int) {
	return 5 + N*(Level+1)<<2
}

// MarshalBinarySize32 returns the number of bytes the polynomial will take when written to data.
// Assumes that each coefficient is encoded on 4 bytes.
func (pol *Poly) MarshalBinarySize32() (cnt int) {
	return MarshalBinarySize32(pol.N(), pol.Level())
}

// Decode32 decodes a slice of bytes in the target polynomial returns the number of bytes decoded.
// The method will first try to write on the buffer. If this step fails, either because the buffer isn't
// allocated or because it is of the wrong size, the method will allocate the correct buffer.
// Assumes that each coefficient is encoded on 8 bytes.
func (pol *Poly) Decode32(data []byte) (ptr int, err error) {

	N := int(binary.BigEndian.Uint32(data))
	Level := int(data[4])

	ptr = 5

	if pol.Buff == nil || len(pol.Buff) != N*(Level+1) {
		pol.Buff = make([]uint64, N*(Level+1))
	}

	if ptr, err = Decode32(ptr, pol.Buff, data); err != nil {
		return ptr, err
	}

	pol.Coeffs = make([][]uint64, Level+1)

	for i := 0; i < Level+1; i++ {
		pol.Coeffs[i] = pol.Buff[i*N : (i+1)*N]
	}

	return ptr, nil
}

// Decode32 converts a byte array to a matrix of coefficients.
// Assumes that each coefficient is encoded on 4 bytes.
func Decode32(ptr int, coeffs []uint64, data []byte) (int, error) {
	for i, j := 0, ptr; i < len(coeffs); i, j = i+1, j+4 {
		coeffs[i] = uint64(binary.BigEndian.Uint32(data[j:]))
	}

	return ptr + len(coeffs)*4, nil
}
