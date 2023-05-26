package bgv

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = ring.GaloisGen

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type Encoder struct {
	params Parameters

	indexMatrix []uint64

	buffQ *ring.Poly
	buffT *ring.Poly

	paramsQP []ring.ModUpConstants
	qHalf    []*big.Int

	tInvModQ []*big.Int
}

// NewEncoder creates a new Encoder from the provided parameters.
func NewEncoder(params Parameters) *Encoder {

	var N, pow, pos uint64 = uint64(params.N()), 1, 0

	logN := params.LogN()

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

	paramsQP := make([]ring.ModUpConstants, ringQ.ModuliChainLength())

	qHalf := make([]*big.Int, ringQ.ModuliChainLength())

	moduli := ringQ.ModuliChain()
	T := ringT.SubRings[0].Modulus

	for i := 1; i < ringQ.ModuliChainLength(); i++ {
		paramsQP[i] = ring.GenModUpConstants(moduli[:i+1], []uint64{T})
		qHalf[i] = new(big.Int).Set(ringQ.ModulusAtLevel[i])
		qHalf[i].Rsh(qHalf[i], 1)
	}

	tInvModQ := make([]*big.Int, ringQ.ModuliChainLength())
	for i := range moduli {
		tInvModQ[i] = bignum.NewInt(T)
		tInvModQ[i].ModInverse(tInvModQ[i], ringQ.ModulusAtLevel[i])
	}

	return &Encoder{
		params:      params,
		indexMatrix: indexMatrix,
		buffQ:       ringQ.NewPoly(),
		buffT:       ringT.NewPoly(),
		paramsQP:    paramsQP,
		qHalf:       qHalf,
		tInvModQ:    tInvModQ,
	}
}

// Parameters returns the underlying parameters of the Encoder.
func (ecd *Encoder) Parameters() Parameters {
	return ecd.params
}

// EncodeNew encodes a slice of integers of type []uint64 or []int64 of size at most N on a newly allocated plaintext.
func (ecd *Encoder) EncodeNew(values interface{}, level int, scale rlwe.Scale) (pt *rlwe.Plaintext) {
	pt = NewPlaintext(ecd.params, level)
	pt.Scale = scale
	ecd.Encode(values, pt)
	return
}

// Encode encodes a slice of integers of type []uint64 or []int64 of size at most N into a pre-allocated plaintext.
func (ecd *Encoder) Encode(values interface{}, pt *rlwe.Plaintext) {
	ecd.EncodeRingT(values, pt.Scale, ecd.buffT)
	ecd.RingT2Q(pt.Level(), true, ecd.buffT, pt.Value)

	if pt.IsNTT {
		ecd.params.RingQ().AtLevel(pt.Level()).NTT(pt.Value, pt.Value)
	}
}

// EncodeCoeffs encodes a slice of []uint64 of size at most N on a pre-allocated plaintext.
// The encoding is done coefficient wise, i.e. [1, 2, 3, 4] -> 1 + 2X + 3X^2 + 4X^3.
func (ecd *Encoder) EncodeCoeffs(values []uint64, pt *rlwe.Plaintext) {
	copy(ecd.buffT.Coeffs[0], values)

	N := len(ecd.buffT.Coeffs[0])

	for i := len(values); i < N; i++ {
		ecd.buffT.Coeffs[0][i] = 0
	}

	ringT := ecd.params.RingT()

	ringT.MulScalar(ecd.buffT, pt.Scale.Uint64(), ecd.buffT)
	ecd.RingT2Q(pt.Level(), true, ecd.buffT, pt.Value)

	if pt.IsNTT {
		ecd.params.RingQ().AtLevel(pt.Level()).NTT(pt.Value, pt.Value)
	}
}

// EncodeCoeffsNew encodes a slice of []uint64 of size at most N on a newly allocated plaintext.
// The encoding is done coefficient wise, i.e. [1, 2, 3, 4] -> 1 + 2X + 3X^2 + 4X^3.}
func (ecd *Encoder) EncodeCoeffsNew(values []uint64, level int, scale rlwe.Scale) (pt *rlwe.Plaintext) {
	pt = NewPlaintext(ecd.params, level)
	pt.Scale = scale
	ecd.EncodeCoeffs(values, pt)
	return
}

// EncodeRingT encodes a slice of []uint64 or []int64 on a polynomial in basis T.
func (ecd *Encoder) EncodeRingT(values interface{}, scale rlwe.Scale, pT *ring.Poly) {

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

		T := ringT.SubRings[0].Modulus
		BRC := ringT.SubRings[0].BRedConstant

		var sign, abs uint64
		for i, c := range values {
			sign = uint64(c) >> 63
			abs = ring.BRedAdd(uint64(c*((int64(sign)^1)-int64(sign))), T, BRC)
			pt[ecd.indexMatrix[i]] = sign*(T-abs) | (sign^1)*abs
		}
		valLen = len(values)
	default:
		panic("cannot EncodeRingT: values must be either []uint64 or []int64")
	}

	N := len(ecd.indexMatrix)
	for i := valLen; i < N; i++ {
		pt[ecd.indexMatrix[i]] = 0
	}

	ringT.INTT(pT, pT)
	ringT.MulScalar(pT, scale.Uint64(), pT)
}

// DecodeRingT decodes a pT in basis T on a slice of []uint64 or []int64.
func (ecd *Encoder) DecodeRingT(pT *ring.Poly, scale rlwe.Scale, values interface{}) {
	ringT := ecd.params.RingT()
	ringT.MulScalar(pT, ring.ModExp(scale.Uint64(), ringT.SubRings[0].Modulus-2, ringT.SubRings[0].Modulus), ecd.buffT)
	ringT.NTT(ecd.buffT, ecd.buffT)

	tmp := ecd.buffT.Coeffs[0]

	N := ringT.N()

	switch values := values.(type) {
	case []uint64:
		for i := 0; i < N; i++ {
			values[i] = tmp[ecd.indexMatrix[i]]
		}
	case []int64:
		modulus := int64(ecd.params.T())
		modulusHalf := modulus >> 1
		var value int64
		for i := 0; i < N; i++ {
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
// If scaleUp, then scales pQ by T^-1 mod Q (or Q/T if T|Q).
func (ecd *Encoder) RingT2Q(level int, scaleUp bool, pT, pQ *ring.Poly) {

	for i := 0; i < level+1; i++ {
		copy(pQ.Coeffs[i], pT.Coeffs[0])
	}

	if scaleUp {
		ecd.params.RingQ().AtLevel(level).MulScalarBigint(pQ, ecd.tInvModQ[level], pQ)
	}
}

// RingQ2T takes pQ in base Q and returns it in base T on pT.
// If scaleDown, scales first pQ by T.
func (ecd *Encoder) RingQ2T(level int, scaleDown bool, pQ, pT *ring.Poly) {

	ringQ := ecd.params.RingQ().AtLevel(level)
	ringT := ecd.params.RingT()

	var poly *ring.Poly
	if scaleDown {
		ringQ.MulScalar(pQ, ecd.params.T(), ecd.buffQ)
		poly = ecd.buffQ
	} else {
		poly = pQ
	}

	if level > 0 {
		ringQ.AddScalarBigint(poly, ecd.qHalf[level], ecd.buffQ)
		ring.ModUpExact(ecd.buffQ.Coeffs[:level+1], pT.Coeffs, ringQ, ringT, ecd.paramsQP[level])
		ringT.SubScalarBigint(pT, ecd.qHalf[level], pT)
	} else {
		ringQ.AddScalar(poly, ringQ.SubRings[0].Modulus>>1, ecd.buffQ)
		ringT.Reduce(ecd.buffQ, pT)
		ringT.SubScalar(pT, ring.BRedAdd(ringQ.SubRings[0].Modulus>>1, ringT.SubRings[0].Modulus, ringT.SubRings[0].BRedConstant), pT)
	}
}

// DecodeUint decodes a any plaintext type and write the coefficients on an pre-allocated uint64 slice.
func (ecd *Encoder) DecodeUint(pt *rlwe.Plaintext, values []uint64) {

	if pt.IsNTT {
		ecd.params.RingQ().AtLevel(pt.Level()).INTT(pt.Value, ecd.buffQ)
	}

	ecd.RingQ2T(pt.Level(), true, ecd.buffQ, ecd.buffT)
	ecd.DecodeRingT(ecd.buffT, pt.Scale, values)
}

// DecodeUintNew decodes any plaintext type and returns the coefficients on a new []uint64 slice.
func (ecd *Encoder) DecodeUintNew(pt *rlwe.Plaintext) (values []uint64) {
	values = make([]uint64, ecd.params.N())
	ecd.DecodeUint(pt, values)
	return
}

// DecodeInt decodes a any plaintext type and write the coefficients on an pre-allocated int64 slice.
// Values are centered between [t/2, t/2).
func (ecd *Encoder) DecodeInt(pt *rlwe.Plaintext, values []int64) {

	if pt.IsNTT {
		ecd.params.RingQ().AtLevel(pt.Level()).INTT(pt.Value, ecd.buffQ)
	}

	ecd.RingQ2T(pt.Level(), true, ecd.buffQ, ecd.buffT)
	ecd.DecodeRingT(ecd.buffT, pt.Scale, values)
}

// DecodeIntNew decodes a any plaintext type and write the coefficients on an new int64 slice.
// Values are centered between [t/2, t/2).
func (ecd *Encoder) DecodeIntNew(pt *rlwe.Plaintext) (values []int64) {
	values = make([]int64, ecd.params.N())
	ecd.DecodeInt(pt, values)
	return
}

func (ecd *Encoder) DecodeCoeffs(pt *rlwe.Plaintext, values []uint64) {

	if pt.IsNTT {
		ecd.params.RingQ().AtLevel(pt.Level()).INTT(pt.Value, ecd.buffQ)
	}

	ecd.RingQ2T(pt.Level(), true, ecd.buffQ, ecd.buffT)
	ringT := ecd.params.RingT()
	ringT.MulScalar(ecd.buffT, ring.ModExp(pt.Scale.Uint64(), ringT.SubRings[0].Modulus-2, ringT.SubRings[0].Modulus), ecd.buffT)
	copy(values, ecd.buffT.Coeffs[0])
}

func (ecd *Encoder) DecodeCoeffsNew(pt *rlwe.Plaintext) (values []uint64) {
	values = make([]uint64, ecd.params.N())
	ecd.DecodeCoeffs(pt, values)
	return
}

// ShallowCopy creates a shallow copy of Encoder in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encoder can be used concurrently.
func (ecd *Encoder) ShallowCopy() *Encoder {
	return &Encoder{
		params:      ecd.params,
		indexMatrix: ecd.indexMatrix,
		buffQ:       ecd.params.RingQ().NewPoly(),
		buffT:       ecd.params.RingT().NewPoly(),
		paramsQP:    ecd.paramsQP,
		qHalf:       ecd.qHalf,
		tInvModQ:    ecd.tInvModQ,
	}
}
