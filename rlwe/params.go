// Package rlwe implements the generic operations that are common to R-LWE schemes. The other implemented schemes extend this package with their specific operations and structures.
package rlwe

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// MaxLogN is the log2 of the largest supported polynomial modulus degree.
const MaxLogN = 17

// MinLogN is the log2 of the smallest supported polynomial modulus degree (needed to ensure the NTT correctness).
const MinLogN = 4

// MaxModuliCount is the largest supported number of moduli in the RNS representation.
const MaxModuliCount = 34

// MaxModuliSize is the largest bit-length supported for the moduli in the RNS representation.
const MaxModuliSize = 60

// DefaultSigma is the default error distribution standard deviation
const DefaultSigma = 3.2

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = ring.GaloisGen

// ParametersLiteral is a literal representation of BFV parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The NewParametersFromLiteral function is used to generate the actual
// checked parameters from the literal representation.
//
// Users must set the polynomial degree (LogN) and the coefficient modulus, by either setting
// the Q and P fields to the desired moduli chain, or by setting the LogQ and LogP fields to
// the desired moduli sizes.
//
// Optionally, users may specify
// - the base 2 decomposition for the gadget ciphertexts
// - the error variance (Sigma) and secrets' density (H) and the ring
// type (RingType). If left unset, standard default values for these field are substituted at
// parameter creation (see NewParametersFromLiteral).
type ParametersLiteral struct {
	LogN           int
	Q              []uint64
	P              []uint64
	LogQ           []int `json:",omitempty"`
	LogP           []int `json:",omitempty"`
	Pow2Base       int
	Sigma          float64
	H              int
	RingType       ring.Type
	DefaultScale   Scale
	DefaultNTTFlag bool
}

// Parameters represents a set of generic RLWE parameters. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	logN           int
	qi             []uint64
	pi             []uint64
	pow2Base       int
	sigma          float64
	h              int
	ringQ          *ring.Ring
	ringP          *ring.Ring
	ringType       ring.Type
	defaultScale   Scale
	defaultNTTFlag bool
}

// NewParameters returns a new set of generic RLWE parameters from the given ring degree logn, moduli q and p, and
// error distribution parameter sigma. It returns the empty parameters Parameters{} and a non-nil error if the
// specified parameters are invalid.
func NewParameters(logn int, q, p []uint64, pow2Base, h int, sigma float64, ringType ring.Type, defaultScale Scale, defaultNTTFlag bool) (Parameters, error) {

	if pow2Base != 0 && len(p) > 1 {
		return Parameters{}, fmt.Errorf("rlwe.NewParameters: invalid parameters, cannot have pow2Base > 0 if len(P) > 1")
	}

	var lenP int
	if p != nil {
		lenP = len(p)
	}

	var err error
	if err = checkSizeParams(logn, len(q), lenP); err != nil {
		return Parameters{}, err
	}

	params := Parameters{
		logN:           logn,
		qi:             make([]uint64, len(q)),
		pi:             make([]uint64, lenP),
		pow2Base:       pow2Base,
		h:              h,
		sigma:          sigma,
		ringType:       ringType,
		defaultScale:   defaultScale,
		defaultNTTFlag: defaultNTTFlag,
	}

	// pre-check that moduli chain is of valid size and that all factors are prime.
	// note: the Ring instantiation checks that the moduli are valid NTT-friendly primes.
	if err = CheckModuli(q, p); err != nil {
		return Parameters{}, err
	}

	copy(params.qi, q)

	if p != nil {
		copy(params.pi, p)
	}

	return params, params.initRings()
}

// NewParametersFromLiteral instantiate a set of generic RLWE parameters from a ParametersLiteral specification.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
//
// If the moduli chain is specified through the LogQ and LogP fields, the method generates a moduli chain matching
// the specified sizes (see `GenModuli`).
//
// If the secrets' density parameter (H) is left unset, its value is set to 2^(paramDef.LogN-1) to match
// the standard ternary distribution.
//
// If the error variance is left unset, its value is set to `DefaultSigma`.
//
// If the RingType is left unset, the default value is ring.Standard.
func NewParametersFromLiteral(paramDef ParametersLiteral) (Parameters, error) {

	if paramDef.H == 0 {
		paramDef.H = 1 << (paramDef.LogN - 1)
	}

	if paramDef.Sigma == 0 {
		// prevents the zero value of ParameterLiteral to result in a noise-less parameter instance.
		// Users should use the NewParameters method to explicitely create noiseless instances.
		paramDef.Sigma = DefaultSigma
	}

	if paramDef.DefaultScale.Cmp(Scale{}) == 0 {
		paramDef.DefaultScale = NewScale(1)
	}

	switch {
	case paramDef.Q != nil && paramDef.LogQ == nil:
		return NewParameters(paramDef.LogN, paramDef.Q, paramDef.P, paramDef.Pow2Base, paramDef.H, paramDef.Sigma, paramDef.RingType, paramDef.DefaultScale, paramDef.DefaultNTTFlag)
	case paramDef.LogQ != nil && paramDef.Q == nil:
		var q, p []uint64
		var err error
		switch paramDef.RingType {
		case ring.Standard:
			q, p, err = GenModuli(paramDef.LogN, paramDef.LogQ, paramDef.LogP)
		case ring.ConjugateInvariant:
			q, p, err = GenModuli(paramDef.LogN+1, paramDef.LogQ, paramDef.LogP)
		default:
			return Parameters{}, fmt.Errorf("rlwe.NewParametersFromLiteral: invalid ring.Type, must be ring.ConjugateInvariant or ring.Standard")
		}
		if err != nil {
			return Parameters{}, err
		}
		return NewParameters(paramDef.LogN, q, p, paramDef.Pow2Base, paramDef.H, paramDef.Sigma, paramDef.RingType, paramDef.DefaultScale, paramDef.DefaultNTTFlag)
	default:
		return Parameters{}, fmt.Errorf("rlwe.NewParametersFromLiteral: invalid parameter literal")
	}
}

// StandardParameters returns a RLWE parameter set that corresponds to the
// standard dual of a conjugate invariant parameter set. If the receiver is already
// a standard set, then the method returns the receiver.
func (p Parameters) StandardParameters() (pci Parameters, err error) {

	switch p.ringType {
	case ring.Standard:
		return p, nil
	case ring.ConjugateInvariant:
		pci = p
		pci.logN = p.logN + 1
		pci.ringType = ring.Standard
		err = pci.initRings()
	default:
		err = fmt.Errorf("invalid ring type")
	}

	return
}

// ParametersLiteral returns the ParametersLiteral of the target Parameters.
func (p Parameters) ParametersLiteral() ParametersLiteral {

	Q := make([]uint64, len(p.qi))
	copy(Q, p.qi)

	P := make([]uint64, len(p.pi))
	copy(P, p.pi)

	return ParametersLiteral{
		LogN:           p.logN,
		Q:              Q,
		P:              P,
		Pow2Base:       p.pow2Base,
		Sigma:          p.sigma,
		H:              p.h,
		RingType:       p.ringType,
		DefaultScale:   p.defaultScale,
		DefaultNTTFlag: p.defaultNTTFlag,
	}
}

// N returns the ring degree
func (p Parameters) N() int {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p Parameters) LogN() int {
	return p.logN
}

// RingQ returns a pointer to ringQ
func (p Parameters) RingQ() *ring.Ring {
	return p.ringQ
}

// RingP returns a pointer to ringP
func (p Parameters) RingP() *ring.Ring {
	return p.ringP
}

// RingQP returns a pointer to ringQP
func (p Parameters) RingQP() *ringqp.Ring {
	return &ringqp.Ring{RingQ: p.ringQ, RingP: p.ringP}
}

// DefaultScale returns the default scale, if any.
func (p Parameters) DefaultScale() Scale {
	return p.defaultScale
}

// DefaultNTTFlag returns the default NTT flag.
func (p Parameters) DefaultNTTFlag() bool {
	return p.defaultNTTFlag
}

// HammingWeight returns the number of non-zero coefficients in secret-keys.
func (p Parameters) HammingWeight() int {
	return p.h
}

// Sigma returns standard deviation of the noise distribution
func (p Parameters) Sigma() float64 {
	return p.sigma
}

// NoiseBound returns truncation bound for the noise distribution.
func (p Parameters) NoiseBound() uint64 {
	return uint64(math.Floor(p.sigma * 6))
}

// RingType returns the type of the underlying ring.
func (p Parameters) RingType() ring.Type {
	return p.ringType
}

// MaxLevel returns the maximum level of a ciphertext
func (p Parameters) MaxLevel() int {
	return p.QCount() - 1
}

// Q returns a new slice with the factors of the ciphertext modulus q
func (p Parameters) Q() []uint64 {
	qi := make([]uint64, len(p.qi))
	copy(qi, p.qi)
	return qi
}

// QiFloat64 returns the float64 value of the Qi at position level in the modulus chain.
func (p Parameters) QiFloat64(level int) float64 {
	return float64(p.qi[level])
}

// QCount returns the number of factors of the ciphertext modulus Q
func (p Parameters) QCount() int {
	return len(p.qi)
}

// QBigInt return the ciphertext-space modulus Q in big.Integer, reconstructed, representation.
func (p Parameters) QBigInt() *big.Int {
	q := big.NewInt(1)
	for _, qi := range p.qi {
		q.Mul(q, new(big.Int).SetUint64(qi))
	}
	return q
}

// P returns a new slice with the factors of the ciphertext modulus extension P
func (p Parameters) P() []uint64 {
	pi := make([]uint64, len(p.pi))
	copy(pi, p.pi)
	return pi
}

// PCount returns the number of factors of the ciphertext modulus extension P
func (p Parameters) PCount() int {
	return len(p.pi)
}

// PBigInt return the ciphertext-space extention modulus P in big.Integer, reconstructed, representation.
func (p Parameters) PBigInt() *big.Int {
	pInt := big.NewInt(1)
	for _, pi := range p.pi {
		pInt.Mul(pInt, new(big.Int).SetUint64(pi))
	}
	return pInt
}

// QP return the extended ciphertext-space modulus QP in RNS representation.
func (p Parameters) QP() []uint64 {
	qp := make([]uint64, len(p.qi)+len(p.pi))
	copy(qp, p.qi)
	copy(qp[len(p.qi):], p.pi)
	return qp
}

// QPCount returns the number of factors of the ciphertext modulus + the modulus extension P
func (p Parameters) QPCount() int {
	return len(p.qi) + len(p.pi)
}

// QPBigInt return the extended ciphertext-space modulus QP in big.Integer, reconstructed, representation.
func (p Parameters) QPBigInt() *big.Int {
	pqInt := p.QBigInt()
	pqInt.Mul(pqInt, p.PBigInt())
	return pqInt
}

// LogQ returns the size of the extended modulus Q in bits
func (p Parameters) LogQ() int {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	return tmp.BitLen()
}

// LogP returns the size of the extended modulus P in bits
func (p Parameters) LogP() int {
	tmp := ring.NewUint(1)
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return tmp.BitLen()
}

// LogQP returns the size of the extended modulus QP in bits
func (p Parameters) LogQP() int {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return tmp.BitLen()
}

// Pow2Base returns the base 2^x decomposition used for the key-switching keys.
// Returns 0 if no decomposition is used (the case where x = 0).
func (p Parameters) Pow2Base() int {
	return p.pow2Base
}

// MaxBit returns max(max(bitLen(Q[:levelQ+1])), max(bitLen(P[:levelP+1])).
func (p Parameters) MaxBit(levelQ, levelP int) (c int) {
	for _, qi := range p.Q()[:levelQ+1] {
		c = utils.MaxInt(c, bits.Len64(qi))
	}

	for _, pi := range p.P()[:levelP+1] {
		c = utils.MaxInt(c, bits.Len64(pi))
	}
	return
}

// DecompPw2 returns ceil(p.MaxBitQ(levelQ, levelP)/bitDecomp).
func (p Parameters) DecompPw2(levelQ, levelP int) (c int) {
	if p.pow2Base == 0 {
		return 1
	}

	return (p.MaxBit(levelQ, levelP) + p.pow2Base - 1) / p.pow2Base
}

// DecompRNS returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (p Parameters) DecompRNS(levelQ, levelP int) int {

	if levelP == -1 {
		return levelQ + 1
	}

	return (levelQ + levelP + 1) / (levelP + 1)
}

// QiOverflowMargin returns floor(2^64 / max(Qi)), i.e. the number of times elements of Z_max{Qi} can
// be added together before overflowing 2^64.
func (p *Parameters) QiOverflowMargin(level int) int {
	return int(math.Exp2(64) / float64(utils.MaxSliceUint64(p.qi[:level+1])))
}

// PiOverflowMargin returns floor(2^64 / max(Pi)), i.e. the number of times elements of Z_max{Pi} can
// be added together before overflowing 2^64.
func (p *Parameters) PiOverflowMargin(level int) int {
	return int(math.Exp2(64) / float64(utils.MaxSliceUint64(p.pi[:level+1])))
}

// GaloisElementForColumnRotationBy returns the Galois element for plaintext
// column rotations by k position to the left. Providing a negative k is
// equivalent to a right rotation.
func (p Parameters) GaloisElementForColumnRotationBy(k int) uint64 {
	return ring.ModExp(GaloisGen, uint64(k&int(p.ringQ.NthRoot-1)), p.ringQ.NthRoot)
}

// GaloisElementForRowRotation returns the Galois element for generating the row
// rotation automorphism
func (p Parameters) GaloisElementForRowRotation() uint64 {
	if p.ringType == ring.ConjugateInvariant {
		panic("Cannot generate GaloisElementForRowRotation if ringType is ConjugateInvariant")
	}
	return p.ringQ.NthRoot - 1
}

// GaloisElementsForTrace generates the Galois elements for the Trace evaluation.
// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 2^{LogN} <= i < N.
func (p Parameters) GaloisElementsForTrace(logN int) (galEls []uint64) {

	galEls = []uint64{}
	for i, j := logN, 0; i < p.LogN()-1; i, j = i+1, j+1 {
		galEls = append(galEls, p.GaloisElementForColumnRotationBy(1<<i))
	}

	if logN == 0 {
		switch p.ringType {
		case ring.Standard:
			galEls = append(galEls, p.GaloisElementForRowRotation())
		case ring.ConjugateInvariant:
			panic("cannot GaloisElementsForTrace: Galois element 5^-1 is undefined in ConjugateInvariant Ring")
		default:
			panic("cannot GaloisElementsForTrace: invalid ring type")
		}
	}

	return
}

// RotationsForReplicate generates the rotations that will be performed by the
// `Evaluator.Replicate` operation when performed with parameters `batch` and `n`.
func (p Parameters) RotationsForReplicate(batch, n int) (rotations []int) {
	return p.RotationsForInnerSum(-batch, n)
}

// RotationsForInnerSum generates the rotations that will be performed by the
// `Evaluator.RotationsForInnerSum` operation when performed with parameters `batch` and `n`.
func (p Parameters) RotationsForInnerSum(batch, n int) (rotations []int) {

	rotIndex := make(map[int]bool)

	var k int
	for i := 1; i < n; i <<= 1 {

		k = i
		k *= batch
		rotIndex[k] = true

		k = n - (n & ((i << 1) - 1))
		k *= batch
		rotIndex[k] = true
	}

	rotations = make([]int, len(rotIndex))
	var i int
	for j := range rotIndex {
		rotations[i] = j
		i++
	}

	return
}

// GaloisElementsForRowInnerSum returns a list of all Galois elements required to
// perform an InnerSum operation. This corresponds to all the left rotations by
// k-positions where k is a power of two and the row-rotation element.
func (p Parameters) GaloisElementsForRowInnerSum() (galEls []uint64) {
	galEls = make([]uint64, p.logN)
	for i := 0; i < int(p.logN)-1; i++ {
		galEls[i] = p.GaloisElementForColumnRotationBy(1 << i)
	}

	switch p.ringType {
	case ring.Standard:
		galEls[p.logN-1] = p.GaloisElementForRowRotation()
	case ring.ConjugateInvariant:
		galEls[p.logN-1] = p.GaloisElementForColumnRotationBy(1 << (p.logN - 1))
	default:
		panic("cannot GaloisElementsForRowInnerSum: invalid ring type")
	}

	return galEls
}

// GaloisElementForExpand returns the list of Galois elements required
// to perform the Expand operation.
func (p Parameters) GaloisElementForExpand(logN int) (galEls []uint64) {
	galEls = make([]uint64, logN)

	for i := 0; i < logN; i++ {
		galEls[i] = uint64(p.N()/(1<<i) + 1)
	}

	return
}

// GaloisElementsForMerge returns the list of Galois elements required
// to perform the Merge operation.
func (p Parameters) GaloisElementsForMerge() (galEls []uint64) {
	return p.GaloisElementsForRowInnerSum()
}

// InverseGaloisElement takes a Galois element and returns the Galois element
// corresponding to the inverse automorphism
func (p Parameters) InverseGaloisElement(galEl uint64) uint64 {
	return ring.ModExp(galEl, p.ringQ.NthRoot-1, p.ringQ.NthRoot)
}

// RotationFromGaloisElement returns the corresponding rotation
// from the Galois element, i.e. computes k given 5^k = galEl mod NthRoot.
func (p Parameters) RotationFromGaloisElement(galEl uint64) (k uint64) {

	fmt.Println(galEl, p.ringQ.NthRoot)
	N := p.ringQ.NthRoot

	x := N >> 3

	for {

		if ring.ModExpPow2(GaloisGen, k, N) != ring.ModExpPow2(galEl, x, N) {
			k |= N >> 3
		}

		if x == 1 {
			return
		}

		x >>= 1
		k >>= 1
	}
}

// Equals checks two Parameter structs for equality.
func (p Parameters) Equals(other Parameters) bool {
	res := p.logN == other.logN
	res = res && utils.EqualSliceUint64(p.qi, other.qi)
	res = res && utils.EqualSliceUint64(p.pi, other.pi)
	res = res && (p.h == other.h)
	res = res && (p.sigma == other.sigma)
	res = res && (p.ringType == other.ringType)
	res = res && (p.defaultScale.Cmp(other.defaultScale) == 0)
	res = res && (p.defaultNTTFlag == other.defaultNTTFlag)
	return res
}

// CopyNew makes a deep copy of the receiver and returns it.
//
// Deprecated: Parameter is now a read-only struct, except for the UnmarshalBinary method: deep copying should only be
// required to save a Parameter struct before calling its UnmarshalBinary method and it will be deprecated when
// transitioning to a immutable serialization interface.
func (p Parameters) CopyNew() Parameters {
	qi, pi := p.qi, p.pi
	p.qi, p.pi = make([]uint64, len(p.qi)), make([]uint64, len(p.pi))
	copy(p.qi, qi)
	p.ringQ, _ = ring.NewRingFromType(1<<p.logN, p.qi, p.ringType)
	if len(p.pi) > 0 {
		copy(p.pi, pi)
		p.ringP, _ = ring.NewRingFromType(1<<p.logN, p.pi, p.ringType)
	}
	return p
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN() == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	// 1 byte : logN
	// 1 byte : #Q
	// 1 byte : #P
	// 1 byte : pow2Base
	// 8 byte : H
	// 8 byte : sigma
	// 1 byte : ringType
	// 1 byte defaultNTTFlag
	// 40 bytes: defaultScale
	// 8 bytes: scaleMod
	// 8 * (#Q) : Q
	// 8 * (#P) : P
	b := utils.NewBuffer(make([]byte, 0, p.MarshalBinarySize()))
	b.WriteUint8(uint8(p.logN))
	b.WriteUint8(uint8(len(p.qi)))
	b.WriteUint8(uint8(len(p.pi)))
	b.WriteUint8(uint8(p.pow2Base))
	b.WriteUint64(uint64(p.h))
	b.WriteUint64(math.Float64bits(p.sigma))
	b.WriteUint8(uint8(p.ringType))
	if p.defaultNTTFlag {
		b.WriteUint8(1)
	} else {
		b.WriteUint8(0)
	}

	data := make([]byte, p.defaultScale.MarshalBinarySize())
	err := p.defaultScale.Encode(data)
	if err != nil {
		return nil, err
	}
	for i := range data {
		b.WriteUint8(data[i])
	}
	if p.defaultScale.Mod != nil {
		b.WriteUint64(p.defaultScale.Mod.Uint64())
	} else {
		b.WriteUint64(0)
	}

	b.WriteUint64Slice(p.qi)
	b.WriteUint64Slice(p.pi)

	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
func (p *Parameters) UnmarshalBinary(data []byte) error {
	if len(data) < 11 {
		return fmt.Errorf("invalid rlwe.Parameter serialization")
	}
	b := utils.NewBuffer(data)
	logN := int(b.ReadUint8())
	lenQ := int(b.ReadUint8())
	lenP := int(b.ReadUint8())
	logbase2 := int(b.ReadUint8())
	h := int(b.ReadUint64())
	sigma := math.Float64frombits(b.ReadUint64())
	ringType := ring.Type(b.ReadUint8())
	var defaultNTTFlag bool
	if b.ReadUint8() == 1 {
		defaultNTTFlag = true
	}

	var defaultScale Scale
	dataScale := make([]uint8, defaultScale.MarshalBinarySize())
	b.ReadUint8Slice(dataScale)
	defaultScale.Decode(dataScale)
	scaleMod := b.ReadUint64()
	if scaleMod != 0 {
		defaultScale.Mod = big.NewInt(0).SetUint64(scaleMod)
	}

	if err := checkSizeParams(logN, lenQ, lenP); err != nil {
		return err
	}

	qi := make([]uint64, lenQ)
	pi := make([]uint64, lenP)
	b.ReadUint64Slice(qi)
	b.ReadUint64Slice(pi)

	var err error
	*p, err = NewParameters(logN, qi, pi, logbase2, h, sigma, ringType, defaultScale, defaultNTTFlag)
	return err
}

// MarshalBinarySize returns the length of the []byte encoding of the reciever.
func (p Parameters) MarshalBinarySize() int {
	return 22 + 48 + (len(p.qi)+len(p.pi))<<3
}

// MarshalJSON returns a JSON representation of this parameter set. See `Marshal` from the `encoding/json` package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	paramsLit := p.ParametersLiteral()
	return json.Marshal(&paramsLit)
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalJSON(data []byte) (err error) {
	var params ParametersLiteral
	if err = json.Unmarshal(data, &params); err != nil {
		return err
	}
	*p, err = NewParametersFromLiteral(params)
	return
}

// CheckModuli checks that the provided q and p correspond to a valid moduli chain.
func CheckModuli(q, p []uint64) error {

	if len(q) > MaxModuliCount {
		return fmt.Errorf("#Qi is larger than %d", MaxModuliCount)
	}

	for i, qi := range q {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("a Qi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, qi := range q {
		if !ring.IsPrime(qi) {
			return fmt.Errorf("a Qi (i=%d) is not a prime", i)
		}
	}

	if p != nil {
		if len(p) > MaxModuliCount {
			return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
		}

		for i, pi := range p {
			if uint64(bits.Len64(pi)-1) > MaxModuliSize+2 {
				return fmt.Errorf("a Pi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
			}
		}

		for i, pi := range p {
			if !ring.IsPrime(pi) {
				return fmt.Errorf("a Pi (i=%d) is not a prime", i)
			}
		}
	}

	return nil
}

func checkSizeParams(logN int, lenQ, lenP int) error {
	if logN > MaxLogN {
		return fmt.Errorf("logN=%d is larger than MaxLogN=%d", logN, MaxLogN)
	}
	if logN < MinLogN {
		return fmt.Errorf("logN=%d is smaller than MinLogN=%d", logN, MinLogN)
	}
	if lenQ > MaxModuliCount {
		return fmt.Errorf("lenQ=%d is larger than MaxModuliCount=%d", lenQ, MaxModuliCount)
	}
	if lenP > MaxModuliCount {
		return fmt.Errorf("lenP=%d is larger than MaxModuliCount=%d", lenP, MaxModuliCount)
	}
	return nil
}

func checkModuliLogSize(logQ, logP []int) error {

	for i, qi := range logQ {
		if qi <= 0 || qi > MaxModuliSize {
			return fmt.Errorf("logQ[%d]=%d is not in ]0, %d]", i, qi, MaxModuliSize)
		}
	}

	for i, pi := range logP {
		if pi <= 0 || pi > MaxModuliSize+1 {
			return fmt.Errorf("logP[%d]=%d is not in ]0,%d]", i, pi, MaxModuliSize+1)
		}
	}

	return nil
}

// GenModuli generates a valid moduli chain from the provided moduli sizes.
func GenModuli(logN int, logQ, logP []int) (q, p []uint64, err error) {

	if err = checkSizeParams(logN, len(logQ), len(logP)); err != nil {
		return
	}

	if err = checkModuliLogSize(logQ, logP); err != nil {
		return
	}

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[int]int)
	for _, qi := range logQ {
		primesbitlen[qi]++
	}

	for _, pj := range logP {
		primesbitlen[pj]++
	}

	// For each bit-size, finds that many primes
	primes := make(map[int][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(int(key), 2<<logN, int(value))
	}

	// Assigns the primes to the moduli chain
	for _, qi := range logQ {
		q = append(q, primes[qi][0])
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the extended ring
	for _, pj := range logP {
		p = append(p, primes[pj][0])
		primes[pj] = primes[pj][1:]
	}

	return
}

func (p *Parameters) initRings() (err error) {
	if p.ringQ, err = ring.NewRingFromType(1<<p.logN, p.qi, p.ringType); err != nil {
		return err
	}
	if len(p.pi) != 0 {
		p.ringP, err = ring.NewRingFromType(1<<p.logN, p.pi, p.ringType)
	}
	return err
}
