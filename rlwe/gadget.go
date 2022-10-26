package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// GadgetCiphertext is a struct for storing an encrypted
// plaintext times the gadget power matrix.
type GadgetCiphertext struct {
	Value [][]CiphertextQP
}

// NewGadgetCiphertext returns a new Ciphertext key with pre-allocated zero-value.
// Ciphertext is always in the NTT domain.
func NewGadgetCiphertext(levelQ, levelP, decompRNS, decompBIT int, ringQP ringqp.Ring) (ct *GadgetCiphertext) {
	ct = new(GadgetCiphertext)
	ct.Value = make([][]CiphertextQP, decompRNS)
	for i := 0; i < decompRNS; i++ {
		ct.Value[i] = make([]CiphertextQP, decompBIT)
		for j := 0; j < decompBIT; j++ {
			ct.Value[i][j].Value[0] = ringQP.NewPolyLvl(levelQ, levelP)
			ct.Value[i][j].Value[1] = ringQP.NewPolyLvl(levelQ, levelP)
			ct.Value[i][j].IsNTT = true
		}
	}

	return ct
}

// LevelQ returns the level of the modulus Q of the target Ciphertext.
func (ct *GadgetCiphertext) LevelQ() int {
	return ct.Value[0][0].Value[0].Q.Level()
}

// LevelP returns the level of the modulus P of the target Ciphertext.
func (ct *GadgetCiphertext) LevelP() int {
	if ct.Value[0][0].Value[0].P != nil {
		return ct.Value[0][0].Value[0].P.Level()
	}

	return -1
}

// Equals checks two Ciphertexts for equality.
func (ct *GadgetCiphertext) Equals(other *GadgetCiphertext) bool {
	if ct == other {
		return true
	}
	if (ct == nil) != (other == nil) {
		return false
	}
	if len(ct.Value) != len(other.Value) {
		return false
	}

	if len(ct.Value[0]) != len(other.Value[0]) {
		return false
	}

	for i := range ct.Value {
		for j, pol := range ct.Value[i] {
			if !pol.Value[0].Equals(other.Value[i][j].Value[0]) && !pol.Value[1].Equals(other.Value[i][j].Value[1]) {
				return false
			}
		}
	}
	return true
}

// CopyNew creates a deep copy of the receiver Ciphertext and returns it.
func (ct *GadgetCiphertext) CopyNew() (ctCopy *GadgetCiphertext) {
	if ct == nil || len(ct.Value) == 0 {
		return nil
	}
	v := make([][]CiphertextQP, len(ct.Value))
	for i := range ct.Value {
		v[i] = make([]CiphertextQP, len(ct.Value[0]))
		for j, el := range ct.Value[i] {
			v[i][j] = *el.CopyNew()
		}
	}
	return &GadgetCiphertext{Value: v}
}

// MarshalBinarySize returns the length in bytes of the target GadgetCiphertext.
func (ct *GadgetCiphertext) MarshalBinarySize() (dataLen int) {

	dataLen = 2

	for i := range ct.Value {
		for _, el := range ct.Value[i] {
			dataLen += el.Value[0].MarshalBinarySize64()
			dataLen += el.Value[1].MarshalBinarySize64()
		}
	}

	return
}

// MarshalBinary encodes the target Ciphertext on a slice of bytes.
func (ct *GadgetCiphertext) MarshalBinary() (data []byte, err error) {
	data = make([]byte, ct.MarshalBinarySize())
	if _, err = ct.Encode(0, data); err != nil {
		return
	}

	return
}

// UnmarshalBinary decodes a slice of bytes on the target Ciphertext.
func (ct *GadgetCiphertext) UnmarshalBinary(data []byte) (err error) {
	if _, err = ct.Decode(data); err != nil {
		return
	}

	return
}

// Encode encodes the target ciphertext on a pre-allocated slice of bytes.
func (ct *GadgetCiphertext) Encode(pointer int, data []byte) (int, error) {

	var err error
	var inc int

	data[pointer] = uint8(len(ct.Value))
	pointer++
	data[pointer] = uint8(len(ct.Value[0]))
	pointer++

	for i := range ct.Value {
		for _, el := range ct.Value[i] {

			if inc, err = el.Value[0].Encode64(data[pointer:]); err != nil {
				return pointer, err
			}
			pointer += inc

			if inc, err = el.Value[1].Encode64(data[pointer:]); err != nil {
				return pointer, err
			}
			pointer += inc
		}
	}

	return pointer, nil
}

// Decode decodes a slice of bytes on the target ciphertext.
func (ct *GadgetCiphertext) Decode(data []byte) (pointer int, err error) {

	decompRNS := int(data[0])
	decompBIT := int(data[1])

	pointer = 2

	ct.Value = make([][]CiphertextQP, decompRNS)

	var inc int

	for i := range ct.Value {

		ct.Value[i] = make([]CiphertextQP, decompBIT)

		for j := range ct.Value[i] {

			if inc, err = ct.Value[i][j].Value[0].Decode64(data[pointer:]); err != nil {
				return
			}
			pointer += inc

			if inc, err = ct.Value[i][j].Value[1].Decode64(data[pointer:]); err != nil {
				return
			}
			pointer += inc
		}
	}

	return
}

// AddPolyTimesGadgetVectorToGadgetCiphertext takes a plaintext polynomial and a list of Ciphertexts and adds the
// plaintext times the RNS and BIT decomposition to the i-th element of the i-th Ciphertexts. This method panics if
// len(cts) > 2.
func AddPolyTimesGadgetVectorToGadgetCiphertext(pt *ring.Poly, cts []GadgetCiphertext, ringQP ringqp.Ring, logbase2 int, buff *ring.Poly) {

	ringQ := ringQP.RingQ
	levelQ := cts[0].LevelQ()
	levelP := cts[0].LevelP()

	if len(cts) > 2 {
		panic("cannot AddPolyTimesGadgetVectorToGadgetCiphertext: len(cts) should be <= 2")
	}

	if levelP != -1 {
		ringQ.MulScalarBigintLvl(levelQ, pt, ringQP.RingP.ModulusAtLevel[levelP], buff) // P * pt
	} else {
		levelP = 0
		if pt != buff {
			ring.CopyLvl(levelQ, pt, buff) // 1 * pt
		}
	}

	RNSDecomp := len(cts[0].Value)
	BITDecomp := len(cts[0].Value[0])

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e + (m * P * w^2j) * (q_star * q_tild) mod QP
			//
			// q_prod = prod(q[i*#Pi+j])
			// q_star = Q/qprod
			// q_tild = q_star^-1 mod q_prod
			//
			// Therefore : (pt * P * w^2j) * (q_star * q_tild) = pt*P*w^2j mod q[i*#Pi+j], else 0
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handle cases where #pj does not divide #qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.Modulus[index]
				p0tmp := buff.Coeffs[index]

				for u, ct := range cts {
					p1tmp := ct.Value[i][j].Value[u].Q.Coeffs[index]
					for w := 0; w < ringQ.N; w++ {
						p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
					}
				}

			}
		}

		// w^2j
		ringQ.MulScalar(buff, 1<<logbase2, buff)
	}
}

// AddPolyToGadgetMatrix takes a plaintext polynomial and a list of ringqp.Poly and adds the
// plaintext times the RNS and BIT decomposition to the list of ringqp.Poly.
func AddPolyToGadgetMatrix(pt *ring.Poly, gm [][]ringqp.Poly, ringQP ringqp.Ring, logbase2 int, buff *ring.Poly) {

	ringQ := ringQP.RingQ

	levelQ := gm[0][0].LevelQ()
	levelP := gm[0][0].LevelP()

	if levelP != -1 {
		ringQ.MulScalarBigintLvl(levelQ, pt, ringQP.RingP.ModulusAtLevel[levelP], buff) // P * pt
	} else {
		levelP = 0
		if pt != buff {
			ring.CopyLvl(levelQ, pt, buff) // 1 * pt
		}
	}

	RNSDecomp := len(gm)
	BITDecomp := len(gm[0])

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e + (m * P * w^2j) * (q_star * q_tild) mod QP
			//
			// q_prod = prod(q[i*#Pi+j])
			// q_star = Q/qprod
			// q_tild = q_star^-1 mod q_prod
			//
			// Therefore : (pt * P * w^2j) * (q_star * q_tild) = pt*P*w^2j mod q[i*#Pi+j], else 0
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handle cases where #pj does not divide #qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.Modulus[index]
				p0tmp := buff.Coeffs[index]

				p1tmp := gm[i][j].Q.Coeffs[index]
				for w := 0; w < ringQ.N; w++ {
					p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
				}
			}
		}

		// w^2j
		ringQ.MulScalar(buff, 1<<logbase2, buff)
	}
}

// GadgetPlaintext stores a RGSW plaintext value.
type GadgetPlaintext struct {
	Value []*ring.Poly
}

// NewGadgetPlaintext creates a new gadget plaintext from value, which can be either uint64, int64 or *ring.Poly.
// Plaintext is returned in the NTT and Mongtomery domain.
func NewGadgetPlaintext(value interface{}, levelQ, levelP, logBase2, decompBIT int, ringQP ringqp.Ring) (pt *GadgetPlaintext) {

	ringQ := ringQP.RingQ

	pt = new(GadgetPlaintext)
	pt.Value = make([]*ring.Poly, decompBIT)

	switch el := value.(type) {
	case uint64:
		pt.Value[0] = ringQ.NewPolyLvl(levelQ)
		for i := range ringQ.Modulus[:levelQ+1] {
			pt.Value[0].Coeffs[i][0] = el
		}
	case int64:
		pt.Value[0] = ringQ.NewPolyLvl(levelQ)
		if el < 0 {
			for i, qi := range ringQ.Modulus[:levelQ+1] {
				pt.Value[0].Coeffs[i][0] = qi - uint64(-el)
			}
		} else {
			for i := range ringQ.Modulus[:levelQ+1] {
				pt.Value[0].Coeffs[i][0] = uint64(el)
			}
		}
	case *ring.Poly:
		pt.Value[0] = el.CopyNew()
	default:
		panic("cannot NewGadgetPlaintext: unsupported type, must be wither uint64 or *ring.Poly")
	}

	if levelP > -1 {
		ringQ.MulScalarBigintLvl(levelQ, pt.Value[0], ringQP.RingP.ModulusAtLevel[levelP], pt.Value[0])
	}

	ringQ.NTTLvl(levelQ, pt.Value[0], pt.Value[0])
	ringQ.MFormLvl(levelQ, pt.Value[0], pt.Value[0])

	for i := 1; i < len(pt.Value); i++ {
		pt.Value[i] = pt.Value[0].CopyNew()
		ringQ.MulByPow2Lvl(levelQ, pt.Value[i], i*logBase2, pt.Value[i])
	}

	return
}
