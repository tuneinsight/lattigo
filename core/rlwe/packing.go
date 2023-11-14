package rlwe

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
)

// Trace maps X -> sum((-1)^i * X^{i*n+1}) for n <= i < N
// Monomial X^k vanishes if k is not divisible by (N/n), otherwise it is multiplied by (N/n).
// Ciphertext is pre-multiplied by (N/n)^-1 to remove the (N/n) factor.
// Examples of full Trace for [0 + 1X + 2X^2 + 3X^3 + 4X^4 + 5X^5 + 6X^6 + 7X^7]
//
// 1.
//
//	  [1 + 2X + 3X^2 + 4X^3 + 5X^4 + 6X^5 + 7X^6 + 8X^7]
//	+ [1 - 6X - 3X^2 + 8X^3 + 5X^4 + 2X^5 - 7X^6 - 4X^7]  {X-> X^(i * 5^1)}
//	= [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//
// 2.
//
//	  [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//	+ [2 + 4X + 0X^2 -12X^3 +10X^4 - 8X^5 + 0X^6 - 4X^7]  {X-> X^(i * 5^2)}
//	= [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//
// 3.
//
//	  [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//	+ [4 + 0X + 0X^2 - 0X^3 -20X^4 + 0X^5 + 0X^6 - 0X^7]  {X-> X^(i * -1)}
//	= [8 + 0X + 0X^2 - 0X^3 + 0X^4 + 0X^5 + 0X^6 - 0X^7]
//
// The method will return an error if the input and output ciphertexts degree is not one.
func (eval Evaluator) Trace(ctIn *Ciphertext, logN int, opOut *Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("ctIn.Degree() != 1 or opOut.Degree() != 1")
	}

	params := eval.GetRLWEParameters()

	level := utils.Min(ctIn.Level(), opOut.Level())

	opOut.Resize(opOut.Degree(), level)

	*opOut.MetaData = *ctIn.MetaData

	gap := 1 << (params.LogN() - logN - 1)

	if logN == 0 {
		gap <<= 1
	}

	if gap > 1 {

		ringQ := params.RingQ().AtLevel(level)

		if ringQ.Type() == ring.ConjugateInvariant {
			gap >>= 1 // We skip the last step that applies phi(5^{-1})
		}

		NInv := new(big.Int).SetUint64(uint64(gap))
		NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

		// pre-multiplication by (N/n)^-1
		ringQ.MulScalarBigint(ctIn.Value[0], NInv, opOut.Value[0])
		ringQ.MulScalarBigint(ctIn.Value[1], NInv, opOut.Value[1])

		if !ctIn.IsNTT {
			ringQ.NTT(opOut.Value[0], opOut.Value[0])
			ringQ.NTT(opOut.Value[1], opOut.Value[1])
			opOut.IsNTT = true
		}

		buff, err := NewCiphertextAtLevelFromPoly(level, []ring.Poly{eval.BuffQP[3].Q, eval.BuffQP[4].Q})

		// Sanity check, this error should not happen unless the
		// evaluator's buffer thave been improperly tempered with.
		if err != nil {
			panic(err)
		}

		buff.IsNTT = true

		for i := logN; i < params.LogN()-1; i++ {

			if err = eval.Automorphism(opOut, params.GaloisElement(1<<i), buff); err != nil {
				return err
			}

			ringQ.Add(opOut.Value[0], buff.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], buff.Value[1], opOut.Value[1])
		}

		if logN == 0 && ringQ.Type() == ring.Standard {

			if err = eval.Automorphism(opOut, ringQ.NthRoot()-1, buff); err != nil {
				return err
			}

			ringQ.Add(opOut.Value[0], buff.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], buff.Value[1], opOut.Value[1])
		}

		if !ctIn.IsNTT {
			ringQ.INTT(opOut.Value[0], opOut.Value[0])
			ringQ.INTT(opOut.Value[1], opOut.Value[1])
			opOut.IsNTT = false
		}

	} else {
		if ctIn != opOut {
			opOut.Copy(ctIn)
		}
	}

	return
}

// GaloisElementsForTrace returns the list of Galois elements requored for the for the `Trace` operation.
// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 2^{LogN} <= i < N.
func GaloisElementsForTrace(params ParameterProvider, logN int) (galEls []uint64) {

	p := params.GetRLWEParameters()

	galEls = []uint64{}
	for i, j := logN, 0; i < p.LogN()-1; i, j = i+1, j+1 {
		galEls = append(galEls, p.GaloisElement(1<<i))
	}

	if logN == 0 {
		switch p.RingType() {
		case ring.Standard:
			galEls = append(galEls, p.GaloisElementOrderTwoOrthogonalSubgroup())
		case ring.ConjugateInvariant:
			panic("cannot GaloisElementsForTrace: Galois element GaloisGen^-1 is undefined in ConjugateInvariant Ring")
		default:
			panic("cannot GaloisElementsForTrace: invalid ring type")
		}
	}

	return
}

// Expand expands a RLWE Ciphertext encrypting sum ai * X^i to 2^logN ciphertexts,
// each encrypting ai * X^0 for 0 <= i < 2^LogN. That is, it extracts the first 2^logN
// coefficients, whose degree is a multiple of 2^logGap, of ctIn and returns an RLWE
// Ciphertext for each coefficient extracted.
//
// The method will return an error if:
//   - The input ciphertext degree is not one
//   - The ring type is not ring.Standard
func (eval Evaluator) Expand(ctIn *Ciphertext, logN, logGap int) (opOut []*Ciphertext, err error) {

	if ctIn.Degree() != 1 {
		return nil, fmt.Errorf("cannot Expand: ctIn.Degree() != 1")
	}

	params := eval.GetRLWEParameters()

	if params.RingType() != ring.Standard {
		return nil, fmt.Errorf("cannot Expand: method is only supported for ring.Type = ring.Standard (X^{-2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	level := ctIn.Level()

	ringQ := params.RingQ().AtLevel(level)

	// Compute X^{-2^{i}} from 1 to LogN
	xPow2 := GenXPow2(ringQ, logN, true)

	opOut = make([]*Ciphertext, 1<<(logN-logGap))
	opOut[0] = ctIn.CopyNew()
	opOut[0].LogDimensions = ring.Dimensions{Rows: 0, Cols: 0}

	if ct := opOut[0]; !ctIn.IsNTT {
		ringQ.NTT(ct.Value[0], ct.Value[0])
		ringQ.NTT(ct.Value[1], ct.Value[1])
		ct.IsNTT = true
	}

	// Multiplies by 2^{-logN} mod Q
	NInv := new(big.Int).SetUint64(1 << logN)
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	ringQ.MulScalarBigint(opOut[0].Value[0], NInv, opOut[0].Value[0])
	ringQ.MulScalarBigint(opOut[0].Value[1], NInv, opOut[0].Value[1])

	gap := 1 << logGap

	tmp, err := NewCiphertextAtLevelFromPoly(level, []ring.Poly{eval.BuffCt.Value[0], eval.BuffCt.Value[1]})

	// Sanity check, this error should not happen unless the
	// evaluator's buffer thave been improperly tempered with.
	if err != nil {
		panic(err)
	}

	tmp.MetaData = ctIn.MetaData

	for i := 0; i < logN; i++ {

		n := 1 << i

		galEl := uint64(ringQ.N()/n + 1)

		half := n / gap

		for j := 0; j < (n+gap-1)/gap; j++ {

			c0 := opOut[j]

			// X -> X^{N/n + 1}
			//[a, b, c, d] -> [a, -b, c, -d]
			if err = eval.Automorphism(c0, galEl, tmp); err != nil {
				return
			}

			if j+half > 0 {

				c1 := opOut[j].CopyNew()

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])

				// Zeroes even coeffs: [a, b, c, d] - [a, -b, c, -d] -> [0, 2b, 0, 2d]
				ringQ.Sub(c1.Value[0], tmp.Value[0], c1.Value[0])
				ringQ.Sub(c1.Value[1], tmp.Value[1], c1.Value[1])

				// c1 * X^{-2^{i}}: [0, 2b, 0, 2d] * X^{-n} -> [2b, 0, 2d, 0]
				ringQ.MulCoeffsMontgomery(c1.Value[0], xPow2[i], c1.Value[0])
				ringQ.MulCoeffsMontgomery(c1.Value[1], xPow2[i], c1.Value[1])

				opOut[j+half] = c1

			} else {

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])
			}
		}
	}

	for _, ct := range opOut {
		if ct != nil && !ctIn.IsNTT {
			ringQ.INTT(ct.Value[0], ct.Value[0])
			ringQ.INTT(ct.Value[1], ct.Value[1])
			ct.IsNTT = false
		}
	}
	return
}

// GaloisElementsForExpand returns the list of Galois elements required
// to perform the `Expand` operation with parameter `logN`.
func GaloisElementsForExpand(params ParameterProvider, logN int) (galEls []uint64) {
	galEls = make([]uint64, logN)

	NthRoot := params.GetRLWEParameters().RingQ().NthRoot()

	for i := 0; i < logN; i++ {
		galEls[i] = uint64(NthRoot/(2<<i) + 1)
	}

	return
}

// Pack packs a batch of RLWE ciphertexts, packing the batch of ciphertexts into a single ciphertext.
// The number of key-switching operations is inputLogGap - log2(gap) + len(cts), where log2(gap) is the
// minimum distance between two keys of the map cts[int]*Ciphertext.
//
// Input:
//
//		cts: a map of Ciphertext, where the index in the map is the future position of the first coefficient
//		      of the indexed ciphertext in the final ciphertext (see example). Ciphertexts can be in or out of the NTT domain.
//		logGap: all coefficients of the input ciphertexts that are not a multiple of X^{2^{logGap}} will be zeroed
//		        during the merging (see example). This is equivalent to skipping the first 2^{logGap} steps of the
//		        algorithm, i.e. having as input ciphertexts that are already partially packed.
//	 zeroGarbageSlots: if set to true, slots which are not multiples of X^{2^{logGap}} will be zeroed during the procedure.
//	                   this will greatly increase the noise and increase the number of key-switching operations to inputLogGap + len(cts).
//
// Output: a ciphertext packing all input ciphertexts
//
// Example: we want to pack 4 ciphertexts into one, and keep only coefficients which are a multiple of X^{4}.
//
//	To do so, we must set logGap = 2.
//	Here the `X` slots are treated as garbage slots that we want to discard during the procedure.
//
//	input: map[int]{
//	   0: [x00, X, X, X, x01, X, X, X],   with logGap = 2
//	   1: [x10, X, X, X, x11, X, X, X],
//	   2: [x20, X, X, X, x21, X, X, X],
//	   3: [x30, X, X, X, x31, X, X, X],
//		}
//
//	 Step 1:
//	         map[0]: 2^{-1} * (map[0] + X^2 * map[2] + phi_{5^2}(map[0] - X^2 * map[2]) = [x00, X, x20, X, x01, X, x21, X]
//	         map[1]: 2^{-1} * (map[1] + X^2 * map[3] + phi_{5^2}(map[1] - X^2 * map[3]) = [x10, X, x30, X, x11, X, x31, X]
//	 Step 2:
//	         map[0]: 2^{-1} * (map[0] + X^1 * map[1] + phi_{5^4}(map[0] - X^1 * map[1]) = [x00, x10, x20, x30, x01, x11, x21, x22]
func (eval Evaluator) Pack(cts map[int]*Ciphertext, inputLogGap int, zeroGarbageSlots bool) (ct *Ciphertext, err error) {

	params := eval.GetRLWEParameters()

	if params.RingType() != ring.Standard {
		return nil, fmt.Errorf("cannot Pack: procedure is only supported for ring.Type = ring.Standard (X^{2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	if len(cts) < 2 {
		return nil, fmt.Errorf("cannot Pack: #cts must be at least 2")
	}

	keys := utils.GetSortedKeys(cts)

	gap := keys[1] - keys[0]
	level := cts[keys[0]].Level()

	for i, key := range keys[1:] {
		level = utils.Min(level, cts[key].Level())

		if i < len(keys)-1 {
			gap = utils.Min(gap, keys[i+1]-keys[i])
		}
	}

	logN := params.LogN()
	ringQ := params.RingQ().AtLevel(level)

	logStart := logN - inputLogGap
	logEnd := logN

	if !zeroGarbageSlots {
		if gap > 0 {
			logEnd -= bits.Len64(uint64(gap - 1))
		}
	}

	if logStart >= logEnd {
		return nil, fmt.Errorf("cannot Pack: gaps between ciphertexts is smaller than inputLogGap > N")
	}

	xPow2 := GenXPow2(ringQ.AtLevel(level), params.LogN(), false) // log(N) polynomial to generate, quick

	NInv := new(big.Int).SetUint64(uint64(1 << (logEnd - logStart)))
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	for _, key := range keys {

		ct := cts[key]

		if ct.Degree() != 1 {
			return nil, fmt.Errorf("cannot Pack: cts[%d].Degree() != 1", key)
		}

		if !ct.IsNTT {
			ringQ.NTT(ct.Value[0], ct.Value[0])
			ringQ.NTT(ct.Value[1], ct.Value[1])
			ct.IsNTT = true
		}

		ringQ.MulScalarBigint(ct.Value[0], NInv, ct.Value[0])
		ringQ.MulScalarBigint(ct.Value[1], NInv, ct.Value[1])
	}

	tmpa := &Ciphertext{}
	tmpa.Value = []ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}
	tmpa.MetaData = &MetaData{}
	tmpa.MetaData.IsNTT = true

	for i := logStart; i < logEnd; i++ {

		t := 1 << (logN - 1 - i)

		for jx, jy := 0, t; jx < t; jx, jy = jx+1, jy+1 {

			a := cts[jx]
			b := cts[jy]

			if b != nil {

				//X^(N/2^L)
				ringQ.MulCoeffsMontgomery(b.Value[0], xPow2[len(xPow2)-i-1], b.Value[0])
				ringQ.MulCoeffsMontgomery(b.Value[1], xPow2[len(xPow2)-i-1], b.Value[1])

				if a != nil {

					// tmpa = phi(a - b * X^{N/2^{i}}, 2^{i-1})
					ringQ.Sub(a.Value[0], b.Value[0], tmpa.Value[0])
					ringQ.Sub(a.Value[1], b.Value[1], tmpa.Value[1])

					// a = a + b * X^{N/2^{i}}
					ringQ.Add(a.Value[0], b.Value[0], a.Value[0])
					ringQ.Add(a.Value[1], b.Value[1], a.Value[1])

				} else {
					// if ct[jx] == nil, then simply re-assigns
					cts[jx] = cts[jy]

					// Required for correctness, since each log step is expected
					// to double the values, which are pre-scaled by N^{-1} mod Q
					// Maybe this can be omitted by doing an individual pre-scaling.
					ringQ.Add(cts[jx].Value[0], cts[jx].Value[0], cts[jx].Value[0])
					ringQ.Add(cts[jx].Value[1], cts[jx].Value[1], cts[jx].Value[1])
				}
			}

			if a != nil {

				var galEl uint64

				if i == 0 {
					galEl = ringQ.NthRoot() - 1
				} else {
					galEl = params.GaloisElement(1 << (i - 1))
				}

				if b != nil {
					if err = eval.Automorphism(tmpa, galEl, tmpa); err != nil {
						return
					}
				} else {
					if err = eval.Automorphism(a, galEl, tmpa); err != nil {
						return
					}
				}

				// a + b * X^{N/2^{i}} + phi(a - b * X^{N/2^{i}}, 2^{i-1})
				ringQ.Add(a.Value[0], tmpa.Value[0], a.Value[0])
				ringQ.Add(a.Value[1], tmpa.Value[1], a.Value[1])
			}
		}
	}

	return cts[0], nil
}

// GaloisElementsForPack returns the list of Galois elements required to perform the `Pack` operation.
func GaloisElementsForPack(params ParameterProvider, logGap int) (galEls []uint64) {

	p := params.GetRLWEParameters()

	// Sanity check
	if logGap > p.LogN() || logGap < 0 {
		panic(fmt.Errorf("cannot GaloisElementsForPack: logGap > logN || logGap < 0"))
	}

	galEls = make([]uint64, 0, logGap)
	for i := 0; i < logGap; i++ {
		galEls = append(galEls, p.GaloisElement(1<<i))
	}

	switch p.RingType() {
	case ring.Standard:
		if logGap == p.LogN() {
			galEls = append(galEls, p.GaloisElementOrderTwoOrthogonalSubgroup())
		}
	default:
		panic("cannot GaloisElementsForPack: invalid ring type")
	}

	return
}

func GenXPow2(r *ring.Ring, logN int, div bool) (xPow []ring.Poly) {

	// Compute X^{-n} from 0 to LogN
	xPow = make([]ring.Poly, logN)

	moduli := r.ModuliChain()[:r.Level()+1]
	BRC := r.BRedConstants()

	var idx int
	for i := 0; i < logN; i++ {

		idx = 1 << i

		if div {
			idx = r.N() - idx
		}

		xPow[i] = r.NewPoly()

		if i == 0 {

			for j := range moduli {
				xPow[i].Coeffs[j][idx] = ring.MForm(1, moduli[j], BRC[j])
			}

			r.NTT(xPow[i], xPow[i])

		} else {
			r.MulCoeffsMontgomery(xPow[i-1], xPow[i-1], xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	if div {
		r.Neg(xPow[0], xPow[0])
	}

	return
}
