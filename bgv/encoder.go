package bgv

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = ring.GaloisGen

// Encoder is an interface for plaintext encoding and decoding operations.
// It provides methods to embed []uint64 and []int64 types into plaintext
// polynomials and the inverse operations.
type Encoder interface {
	Encode(values interface{}, pt *rlwe.Plaintext)
	EncodeNew(values interface{}, level int, scale rlwe.Scale) (pt *rlwe.Plaintext)
	EncodeCoeffs(values []uint64, pt *rlwe.Plaintext)
	EncodeCoeffsNew(values []uint64, level int, scale rlwe.Scale) (pt *rlwe.Plaintext)

	RingT2Q(level int, pT, pQ *ring.Poly)
	RingQ2T(level int, pQ, pT *ring.Poly)

	ScaleUp(level int, pIn, pOut *ring.Poly)
	ScaleDown(level int, pIn, pOut *ring.Poly)

	EncodeRingT(values interface{}, scale rlwe.Scale, pT *ring.Poly)
	DecodeRingT(pT *ring.Poly, scale rlwe.Scale, values interface{})

	DecodeUint(pt *rlwe.Plaintext, values []uint64)
	DecodeInt(pt *rlwe.Plaintext, values []int64)
	DecodeUintNew(pt *rlwe.Plaintext) (values []uint64)
	DecodeIntNew(pt *rlwe.Plaintext) (values []int64)
	DecodeCoeffs(pt *rlwe.Plaintext, values []uint64)
	DecodeCoeffsNew(pt *rlwe.Plaintext) (values []uint64)

	ShallowCopy() Encoder
}

// encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type encoder struct {
	params Parameters

	indexMatrix []uint64

	buffQ *ring.Poly
	buffT *ring.Poly

	paramsQP []ring.ModupParams
	qHalf    []*big.Int

	tInvModQ []*big.Int
}

// NewEncoder creates a new encoder from the provided parameters.
func NewEncoder(params Parameters) Encoder {

	var N, logN, pow, pos uint64 = uint64(params.N()), uint64(params.LogN()), 1, 0

	mask := 2*N - 1

	indexMatrix := make([]uint64, N)

	for i, j := 0, int(N>>1); i < int(N>>1); i, j = i+1, j+1 {

		pos = utils.BitReverse64(pow>>1, logN)

		indexMatrix[i] = pos
		indexMatrix[j] = N - pos - 1

		pow *= GaloisGen
		pow &= mask
	}

	ringQ := params.RingQ()
	ringT := params.RingT()

	paramsQP := make([]ring.ModupParams, ringQ.NbModuli())

	qHalf := make([]*big.Int, ringQ.NbModuli())

	moduli := ringQ.Moduli()
	T := ringT.Tables[0].Modulus

	for i := 1; i < ringQ.NbModuli(); i++ {
		paramsQP[i] = ring.GenModUpParams(moduli[:i+1], []uint64{T})
		qHalf[i] = new(big.Int).Set(ringQ.ModulusAtLevel[i])
		qHalf[i].Rsh(qHalf[i], 1)
	}

	tInvModQ := make([]*big.Int, ringQ.NbModuli())
	for i := range moduli {
		tInvModQ[i] = ring.NewUint(T)
		tInvModQ[i].ModInverse(tInvModQ[i], ringQ.ModulusAtLevel[i])
	}

	return &encoder{
		params:      params,
		indexMatrix: indexMatrix,
		buffQ:       ringQ.NewPoly(),
		buffT:       ringT.NewPoly(),
		paramsQP:    paramsQP,
		qHalf:       qHalf,
		tInvModQ:    tInvModQ,
	}
}

// EncodeNew encodes a slice of integers of type []uint64 or []int64 of size at most N on a newly allocated plaintext.
func (ecd *encoder) EncodeNew(values interface{}, level int, scale rlwe.Scale) (pt *rlwe.Plaintext) {
	pt = NewPlaintext(ecd.params, level)
	pt.Scale = scale
	ecd.Encode(values, pt)
	return
}

// Encode encodes a slice of integers of type []uint64 or []int64 of size at most N into a pre-allocated plaintext.
func (ecd *encoder) Encode(values interface{}, pt *rlwe.Plaintext) {
	ecd.EncodeRingT(values, pt.Scale, ecd.buffT)
	ecd.RingT2Q(pt.Level(), ecd.buffT, pt.Value)

	if pt.IsNTT {
		ecd.params.RingQ().NTTLvl(pt.Level(), pt.Value, pt.Value)
	}

	ecd.ScaleUp(pt.Level(), pt.Value, pt.Value)
}

// EncodeCoeffs encodes a slice of []uint64 of size at most N on a pre-allocated plaintext.
// The encoding is done coefficient wise, i.e. [1, 2, 3, 4] -> 1 + 2X + 3X^2 + 4X^3.
func (ecd *encoder) EncodeCoeffs(values []uint64, pt *rlwe.Plaintext) {
	copy(ecd.buffT.Coeffs[0], values)

	for i := len(values); i < len(ecd.buffT.Coeffs[0]); i++ {
		ecd.buffT.Coeffs[0][i] = 0
	}

	ringT := ecd.params.RingT()

	ringT.MulScalar(ecd.buffT, pt.Scale.Uint64(), ecd.buffT)
	ecd.RingT2Q(pt.Level(), ecd.buffT, pt.Value)

	if pt.IsNTT {
		ecd.params.RingQ().NTTLvl(pt.Level(), pt.Value, pt.Value)
	}

	ecd.ScaleUp(pt.Level(), pt.Value, pt.Value)
}

// EncodeCoeffsNew encodes a slice of []uint64 of size at most N on a newly allocated plaintext.
// The encoding is done coefficient wise, i.e. [1, 2, 3, 4] -> 1 + 2X + 3X^2 + 4X^3.}
func (ecd *encoder) EncodeCoeffsNew(values []uint64, level int, scale rlwe.Scale) (pt *rlwe.Plaintext) {
	pt = NewPlaintext(ecd.params, level)
	pt.Scale = scale
	ecd.EncodeCoeffs(values, pt)
	return
}

// EncodeRingT encodes a slice of []uint64 or []int64 on a polynomial in basis T.
func (ecd *encoder) EncodeRingT(values interface{}, scale rlwe.Scale, pT *ring.Poly) {

	if len(pT.Coeffs[0]) != len(ecd.indexMatrix) {
		panic("cannot EncodeRingT: invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	pt := pT.Coeffs[0]

	ringT := ecd.params.RingT()

	var valLen int
	switch values := values.(type) {
	case []uint64:
		for i, c := range values {
			pt[ecd.indexMatrix[i]] = c
		}
		ringT.Reduce(pT, pT)
		valLen = len(values)
	case []int64:

		T := ringT.Tables[0].Modulus
		bredparamsT := ringT.Tables[0].BRedParams

		var sign, abs uint64
		for i, c := range values {
			sign = uint64(c) >> 63
			abs = ring.BRedAdd(uint64(c*((int64(sign)^1)-int64(sign))), T, bredparamsT)
			pt[ecd.indexMatrix[i]] = sign*(T-abs) | (sign^1)*abs
		}
		valLen = len(values)
	default:
		panic("cannot EncodeRingT: values must be either []uint64 or []int64")
	}

	for i := valLen; i < len(ecd.indexMatrix); i++ {
		pt[ecd.indexMatrix[i]] = 0
	}

	ringT.InvNTT(pT, pT)
	ringT.MulScalar(pT, scale.Uint64(), pT)
}

// EncodeRingT decodes a pT in basis T on a slice of []uint64 or []int64.
func (ecd *encoder) DecodeRingT(pT *ring.Poly, scale rlwe.Scale, values interface{}) {
	ringT := ecd.params.RingT()
	ringT.MulScalar(pT, ring.ModExp(scale.Uint64(), ringT.Tables[0].Modulus-2, ringT.Tables[0].Modulus), ecd.buffT)
	ringT.NTT(ecd.buffT, ecd.buffT)

	tmp := ecd.buffT.Coeffs[0]

	switch values := values.(type) {
	case []uint64:
		for i := 0; i < ecd.params.N(); i++ {
			values[i] = tmp[ecd.indexMatrix[i]]
		}
	case []int64:
		modulus := int64(ecd.params.T())
		modulusHalf := modulus >> 1
		var value int64
		for i := 0; i < ecd.params.N(); i++ {
			if value = int64(tmp[ecd.indexMatrix[i]]); value >= modulusHalf {
				values[i] = value - modulus
			} else {
				values[i] = value
			}
		}
	default:
		panic("cannot DecodeRingT: values must be either []uint64 or []int64")
	}
}

// RingT2Q takes pT in base T and returns it in base Q on pQ.
func (ecd *encoder) RingT2Q(level int, pT, pQ *ring.Poly) {
	for i := 0; i < level+1; i++ {
		copy(pQ.Coeffs[i], pT.Coeffs[0])
	}
}

// ScaleUp scales pIn up T^1 mod Q and returns the result in pOut.
func (ecd *encoder) ScaleUp(level int, pIn, pOut *ring.Poly) {
	ecd.params.RingQ().MulScalarBigintLvl(level, pIn, ecd.tInvModQ[level], pOut)
}

// RingQ2T takes pQ in base Q and returns it in base T on pT.
func (ecd *encoder) RingQ2T(level int, pQ, pT *ring.Poly) {

	ringQ := ecd.params.RingQ()
	ringT := ecd.params.RingT()

	if level > 0 {
		ringQ.AddScalarBigintLvl(level, pQ, ecd.qHalf[level], ecd.buffQ)
		ring.ModUpExact(ecd.buffQ.Coeffs[:level+1], pT.Coeffs, ringQ, ringT, ecd.paramsQP[level])
		ringT.SubScalarBigint(pT, ecd.qHalf[level], pT)
	} else {
		ringQ.AddScalarLvl(level, pQ, ringQ.Tables[0].Modulus>>1, ecd.buffQ)
		ringT.Reduce(ecd.buffQ, pT)
		ringT.SubScalar(pT, ring.BRedAdd(ringQ.Tables[0].Modulus>>1, ringT.Tables[0].Modulus, ringT.Tables[0].BRedParams), pT)
	}
}

// ScaleDown scales pIn down by T and returns the result in pOut.
func (ecd *encoder) ScaleDown(level int, pIn, pOut *ring.Poly) {
	ecd.params.RingQ().MulScalarLvl(level, pIn, ecd.params.T(), pOut)
}

// DecodeUint decodes a any plaintext type and write the coefficients on an pre-allocated uint64 slice.
func (ecd *encoder) DecodeUint(pt *rlwe.Plaintext, values []uint64) {

	if pt.IsNTT {
		ecd.params.RingQ().InvNTTLvl(pt.Level(), pt.Value, ecd.buffQ)
		ecd.ScaleDown(pt.Level(), ecd.buffQ, ecd.buffQ)
	} else {
		ecd.ScaleDown(pt.Level(), pt.Value, ecd.buffQ)
	}

	ecd.RingQ2T(pt.Level(), ecd.buffQ, ecd.buffT)
	ecd.DecodeRingT(ecd.buffT, pt.Scale, values)
}

// DecodeUintNew decodes any plaintext type and returns the coefficients on a new []uint64 slice.
func (ecd *encoder) DecodeUintNew(pt *rlwe.Plaintext) (values []uint64) {
	values = make([]uint64, ecd.params.N())
	ecd.DecodeUint(pt, values)
	return
}

// DecodeInt decodes a any plaintext type and write the coefficients on an pre-allocated int64 slice.
// Values are centered between [t/2, t/2).
func (ecd *encoder) DecodeInt(pt *rlwe.Plaintext, values []int64) {

	if pt.IsNTT {
		ecd.params.RingQ().InvNTTLvl(pt.Level(), pt.Value, ecd.buffQ)
		ecd.ScaleDown(pt.Level(), ecd.buffQ, ecd.buffQ)
	} else {
		ecd.ScaleDown(pt.Level(), pt.Value, ecd.buffQ)
	}

	ecd.RingQ2T(pt.Level(), ecd.buffQ, ecd.buffT)
	ecd.DecodeRingT(ecd.buffT, pt.Scale, values)
}

// DecodeInt decodes a any plaintext type and write the coefficients on an new int64 slice.
// Values are centered between [t/2, t/2).
func (ecd *encoder) DecodeIntNew(pt *rlwe.Plaintext) (values []int64) {
	values = make([]int64, ecd.params.N())
	ecd.DecodeInt(pt, values)
	return
}

func (ecd *encoder) DecodeCoeffs(pt *rlwe.Plaintext, values []uint64) {

	if pt.IsNTT {
		ecd.params.RingQ().InvNTTLvl(pt.Level(), pt.Value, ecd.buffQ)
		ecd.ScaleDown(pt.Level(), ecd.buffQ, ecd.buffQ)
	} else {
		ecd.ScaleDown(pt.Level(), pt.Value, ecd.buffQ)
	}

	ecd.RingQ2T(pt.Level(), ecd.buffQ, ecd.buffT)
	ringT := ecd.params.RingT()
	ringT.MulScalar(ecd.buffT, ring.ModExp(pt.Scale.Uint64(), ringT.Tables[0].Modulus-2, ringT.Tables[0].Modulus), ecd.buffT)
	copy(values, ecd.buffT.Coeffs[0])
}

func (ecd *encoder) DecodeCoeffsNew(pt *rlwe.Plaintext) (values []uint64) {
	values = make([]uint64, ecd.params.N())
	ecd.DecodeCoeffs(pt, values)
	return
}

// ShallowCopy creates a shallow copy of Encoder in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encoder can be used concurrently.
func (ecd *encoder) ShallowCopy() Encoder {
	return &encoder{
		params:      ecd.params,
		indexMatrix: ecd.indexMatrix,
		buffQ:       ecd.params.RingQ().NewPoly(),
		buffT:       ecd.params.RingT().NewPoly(),
		paramsQP:    ecd.paramsQP,
		qHalf:       ecd.qHalf,
		tInvModQ:    ecd.tInvModQ,
	}
}
