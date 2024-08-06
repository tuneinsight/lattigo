package bgv

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
)

type Integer interface {
	int64 | uint64
}

// IntegerSlice is an empty interface whose goal is to
// indicate that the expected input should be []Integer.
// See Integer for information on the type constraint.
type IntegerSlice interface {
}

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = ring.GaloisGen

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type Encoder struct {
	parameters Parameters

	indexMatrix []uint64

	bufQ ring.Poly
	bufT ring.Poly

	// bufB is allocated in the case when the degree of RingT is smaller
	// than the degree of RingQ (gap > 1), hence a more involved conversion
	// between the two structures is necessary. The size of bufB is then
	// MaxSlots() elements.
	bufB []*big.Int

	paramsQP []ring.ModUpConstants
	qHalf    []*big.Int

	tInvModQ []*big.Int
}

// NewEncoder creates a new [Encoder] from the provided parameters.
func NewEncoder(parameters Parameters) *Encoder {

	ringQ := parameters.RingQ()
	ringT := parameters.RingT()

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
	TBig := new(big.Int).SetUint64(T)
	for i := range moduli {
		tInvModQ[i] = new(big.Int).ModInverse(TBig, ringQ.ModulusAtLevel[i])
	}

	var bufB []*big.Int

	if parameters.LogMaxDimensions().Cols < parameters.LogN()-1 {

		slots := parameters.MaxSlots()

		bufB = make([]*big.Int, slots)

		for i := 0; i < slots; i++ {
			bufB[i] = new(big.Int)
		}
	}

	return &Encoder{
		parameters:  parameters,
		indexMatrix: permuteMatrix(parameters.LogMaxSlots()),
		bufQ:        ringQ.NewPoly(),
		bufT:        ringT.NewPoly(),
		bufB:        bufB,
		paramsQP:    paramsQP,
		qHalf:       qHalf,
		tInvModQ:    tInvModQ,
	}
}

func permuteMatrix(logN int) (perm []uint64) {

	var N, pow, pos uint64 = uint64(1 << logN), 1, 0

	mask := 2*N - 1

	perm = make([]uint64, N)

	halfN := int(N >> 1)

	for i, j := 0, halfN; i < halfN; i, j = i+1, j+1 {

		pos = utils.BitReverse64(pow>>1, logN) // = (pow-1)/2

		perm[i] = pos
		perm[j] = N - pos - 1

		pow *= GaloisGen
		pow &= mask
	}

	return perm
}

// GetRLWEParameters returns the underlying [rlwe.Parameters] of the target object.
func (ecd Encoder) GetRLWEParameters() *rlwe.Parameters {
	return &ecd.parameters.Parameters
}

// Encode encodes an [IntegerSlice] of size at most N, where N is the smallest value satisfying PlaintextModulus = 1 mod 2N,
// on a pre-allocated plaintext.
func (ecd Encoder) Encode(values interface{}, pt *rlwe.Plaintext) (err error) {

	if pt.IsBatched {
		return ecd.EmbedScale(values, true, pt.MetaData, pt.Value)
	} else {

		ringT := ecd.parameters.RingT()
		N := ringT.N()
		T := ringT.SubRings[0].Modulus
		BRC := ringT.SubRings[0].BRedConstant

		ptT := ecd.bufT.Coeffs[0]

		var valLen int
		switch values := values.(type) {
		case []uint64:

			if len(values) > N {
				return fmt.Errorf("cannot Encode (TimeDomain): len(values)=%d > N=%d", len(values), N)
			}

			copy(ptT, values)
			valLen = len(values)
		case []int64:

			if len(values) > N {
				return fmt.Errorf("cannot Encode (TimeDomain: len(values)=%d > N=%d", len(values), N)
			}

			var sign, abs uint64
			for i, c := range values {
				sign = uint64(c) >> 63
				abs = ring.BRedAdd(uint64(c*((int64(sign)^1)-int64(sign))), T, BRC)
				ptT[i] = sign*(T-abs) | (sign^1)*abs
			}

			valLen = len(values)
		}

		for i := valLen; i < N; i++ {
			ptT[i] = 0
		}

		ringT.MulScalar(ecd.bufT, pt.Scale.Uint64(), ecd.bufT)
		ecd.RingT2Q(pt.Level(), true, ecd.bufT, pt.Value)

		if pt.IsNTT {
			ecd.parameters.RingQ().AtLevel(pt.Level()).NTT(pt.Value, pt.Value)
		}

		return
	}
}

// EncodeRingT encodes an [IntegerSlice] at the given scale on a polynomial pT with coefficients modulo the plaintext modulus PlaintextModulus.
func (ecd Encoder) EncodeRingT(values IntegerSlice, scale rlwe.Scale, pT ring.Poly) (err error) {
	perm := ecd.indexMatrix

	pt := pT.Coeffs[0]

	ringT := ecd.parameters.RingT()

	slots := pT.N()

	var valLen int
	switch values := values.(type) {
	case []uint64:

		if len(values) > slots {
			return fmt.Errorf("cannot EncodeRingT (FrequencyDomain): len(values)=%d > slots=%d", len(values), slots)
		}

		for i, c := range values {
			pt[perm[i]] = c
		}

		ringT.Reduce(pT, pT)

		valLen = len(values)

	case []int64:

		if len(values) > slots {
			return fmt.Errorf("cannot EncodeRingT (FrequencyDomain): len(values)=%d > slots=%d", len(values), slots)
		}

		T := ringT.SubRings[0].Modulus
		BRC := ringT.SubRings[0].BRedConstant

		var sign, abs uint64
		for i, c := range values {
			sign = uint64(c) >> 63
			abs = ring.BRedAdd(uint64(c*((int64(sign)^1)-int64(sign))), T, BRC)
			pt[perm[i]] = sign*(T-abs) | (sign^1)*abs
		}

		valLen = len(values)
	default:
		return fmt.Errorf("cannot EncodeRingT: values.(type) must be either []uint64 or []int64 but is %T", values)
	}

	// Zeroes the non-mapped coefficients
	N := len(ecd.indexMatrix)
	for i := valLen; i < N; i++ {
		pt[perm[i]] = 0
	}

	// INTT on the Y = X^{N/n}
	ringT.INTT(pT, pT)
	ringT.MulScalar(pT, scale.Uint64(), pT)

	return nil
}

// EmbedScale is a generic method to encode an IntegerSlice on [ringqp.Poly] or *[ring.Poly].
// If scaleUp is true, then the values will to be multiplied by PlaintextModulus^{-1} mod Q after being encoded on the polynomial.
// Encoding is done according to the metadata.
// Accepted polyOut.(type) are a ringqp.Poly and *ring.Poly
func (ecd Encoder) EmbedScale(values IntegerSlice, scaleUp bool, metadata *rlwe.MetaData, polyOut interface{}) (err error) {

	pT := ecd.bufT

	if err = ecd.EncodeRingT(values, metadata.Scale, pT); err != nil {
		return
	}

	// Maps Y = X^{N/n} -> X and quantizes.
	switch p := polyOut.(type) {
	case ringqp.Poly:

		levelQ := p.Q.Level()

		ecd.RingT2Q(levelQ, scaleUp, pT, p.Q)

		ringQ := ecd.parameters.RingQ().AtLevel(levelQ)

		if metadata.IsNTT {
			ringQ.NTT(p.Q, p.Q)
		}

		if metadata.IsMontgomery {
			ringQ.MForm(p.Q, p.Q)
		}

		if p.P.Level() > -1 {

			levelP := p.P.Level()

			ecd.RingT2Q(levelP, scaleUp, pT, p.P)

			ringP := ecd.parameters.RingP().AtLevel(levelP)

			if metadata.IsNTT {
				ringP.NTT(p.P, p.P)
			}

			if metadata.IsMontgomery {
				ringP.MForm(p.P, p.P)
			}
		}

	case ring.Poly:

		level := p.Level()

		ecd.RingT2Q(level, scaleUp, pT, p)

		ringQ := ecd.parameters.RingQ().AtLevel(level)

		if metadata.IsNTT {
			ringQ.NTT(p, p)
		}

		if metadata.IsMontgomery {
			ringQ.MForm(p, p)
		}

	default:
		return fmt.Errorf("cannot embed: invalid polyOut.(Type) must be ringqp.Poly or *ring.Poly")
	}

	return
}

func (ecd Encoder) Embed(values interface{}, metadata *rlwe.MetaData, polyOut interface{}) (err error) {
	return ecd.EmbedScale(values, false, metadata, polyOut)
}

// DecodeRingT decodes a polynomial pT with coefficients modulo the plaintext modulu PlaintextModulus on an InterSlice at the given scale.
func (ecd Encoder) DecodeRingT(pT ring.Poly, scale rlwe.Scale, values IntegerSlice) (err error) {
	ringT := ecd.parameters.RingT()
	ringT.MulScalar(pT, ring.ModExp(scale.Uint64(), ringT.SubRings[0].Modulus-2, ringT.SubRings[0].Modulus), ecd.bufT)
	ringT.NTT(ecd.bufT, ecd.bufT)

	tmp := ecd.bufT.Coeffs[0]

	switch values := values.(type) {
	case []uint64:
		for i := range values {
			values[i] = tmp[ecd.indexMatrix[i]]
		}
	case []int64:
		modulus := int64(ecd.parameters.PlaintextModulus())
		modulusHalf := modulus >> 1
		var value int64
		for i := range values {
			if value = int64(tmp[ecd.indexMatrix[i]]); value >= modulusHalf {
				values[i] = value - modulus
			} else {
				values[i] = value
			}
		}
	default:
		return fmt.Errorf("cannot DecodeRingT: values must be either []uint64 or []int64 but is %T", values)
	}

	return
}

// RingT2Q takes pT in base PlaintextModulus and writes it in base Q[level] on pQ.
// If scaleUp is true, multiplies the values of pQ by PlaintextModulus^{-1} mod Q[level].
func (ecd Encoder) RingT2Q(level int, scaleUp bool, pT, pQ ring.Poly) {

	N := pQ.N()
	n := pT.N()

	gap := N / n

	for i := 0; i < level+1; i++ {

		coeffs := pQ.Coeffs[i]

		copy(coeffs, pT.Coeffs[0])

		if gap > 1 {

			for j := n; j < N; j++ {
				coeffs[j] = 0
			}

			for j := n - 1; j > 0; j-- {
				coeffs[j*gap] = coeffs[j]
				coeffs[j] = 0
			}
		}
	}

	if scaleUp {
		ecd.parameters.RingQ().AtLevel(level).MulScalarBigint(pQ, ecd.tInvModQ[level], pQ)
	}
}

// RingQ2T takes pQ in base Q[level] and writes it in base PlaintextModulus on pT.
// If scaleUp is true, the values of pQ are multiplied by PlaintextModulus mod Q[level]
// before being converted into the base PlaintextModulus.
func (ecd Encoder) RingQ2T(level int, scaleDown bool, pQ, pT ring.Poly) {

	ringQ := ecd.parameters.RingQ().AtLevel(level)
	ringT := ecd.parameters.RingT()

	var poly ring.Poly
	if scaleDown {
		ringQ.MulScalar(pQ, ecd.parameters.PlaintextModulus(), ecd.bufQ)
		poly = ecd.bufQ
	} else {
		poly = pQ
	}

	gap := pQ.N() / pT.N()

	if level > 0 {

		if gap == 1 {
			ringQ.AddScalarBigint(poly, ecd.qHalf[level], ecd.bufQ)
			ring.ModUpExact(ecd.bufQ.Coeffs[:level+1], pT.Coeffs, ringQ, ringT, ecd.paramsQP[level])
			ringT.SubScalarBigint(pT, ecd.qHalf[level], pT)
		} else {
			ringQ.PolyToBigintCentered(poly, gap, ecd.bufB)
			ringT.SetCoefficientsBigint(ecd.bufB, pT)
		}

	} else {

		if gap == 1 {
			ringQ.AddScalar(poly, ringQ.SubRings[0].Modulus>>1, ecd.bufQ)
			ringT.Reduce(ecd.bufQ, pT)
		} else {

			n := pT.N()

			pQCoeffs := poly.Coeffs[0]
			bufQCoeffs := ecd.bufQ.Coeffs[0]

			for i := 0; i < n; i++ {
				bufQCoeffs[i] = pQCoeffs[i*gap]
			}

			ringQ.SubRings[0].AddScalar(bufQCoeffs[:n], ringQ.SubRings[0].Modulus>>1, bufQCoeffs[:n])
			ringT.SubRings[0].Reduce(bufQCoeffs[:n], pT.Coeffs[0])
		}

		ringT.SubScalar(pT, ring.BRedAdd(ringQ.SubRings[0].Modulus>>1, ringT.SubRings[0].Modulus, ringT.SubRings[0].BRedConstant), pT)
	}
}

// Decode decodes a plaintext on an IntegerSlice mod PlaintextModulus of size at most N, where N is the smallest value satisfying PlaintextModulus = 1 mod 2N.
func (ecd Encoder) Decode(pt *rlwe.Plaintext, values interface{}) (err error) {

	bufT := ecd.bufT

	if pt.IsNTT {
		ecd.parameters.RingQ().AtLevel(pt.Level()).INTT(pt.Value, ecd.bufQ)
		ecd.RingQ2T(pt.Level(), true, ecd.bufQ, bufT)
	} else {
		ecd.RingQ2T(pt.Level(), true, pt.Value, bufT)
	}

	if pt.IsBatched {
		return ecd.DecodeRingT(ecd.bufT, pt.Scale, values)
	} else {
		ringT := ecd.parameters.RingT()
		ringT.MulScalar(bufT, ring.ModExp(pt.Scale.Uint64(), ringT.SubRings[0].Modulus-2, ringT.SubRings[0].Modulus), bufT)

		switch values := values.(type) {
		case []uint64:
			copy(values, ecd.bufT.Coeffs[0])
		case []int64:

			ptT := bufT.Coeffs[0]

			N := ecd.parameters.RingT().N()
			modulus := int64(ecd.parameters.PlaintextModulus())
			modulusHalf := modulus >> 1

			var value int64
			for i := 0; i < N; i++ {
				if value = int64(ptT[i]); value >= modulusHalf {
					values[i] = value - modulus
				} else {
					values[i] = value
				}
			}

		default:
			return fmt.Errorf("cannot Decode: values must be either []uint64 or []int64 but is %T", values)
		}

		return
	}
}

// ShallowCopy returns a lightweight copy of the target object
// that can be used concurrently with the original object.
func (ecd Encoder) ShallowCopy() (e *Encoder) {
	e = &Encoder{
		parameters:  ecd.parameters,
		indexMatrix: ecd.indexMatrix,
		bufQ:        ecd.parameters.RingQ().NewPoly(),
		bufT:        ecd.parameters.RingT().NewPoly(),
		paramsQP:    ecd.paramsQP,
		qHalf:       ecd.qHalf,
		tInvModQ:    ecd.tInvModQ,
	}
	for i := 0; ecd.parameters.LogMaxDimensions().Cols < ecd.parameters.LogN()-1 && i < ecd.parameters.MaxSlots(); i++ {
		e.bufB = append(e.bufB, new(big.Int))
	}

	return
}
