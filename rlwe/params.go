// Package rlwe implements the generic operations that are common to R-LWE schemes. The other implemented schemes extend this package with their specific operations and structures.
package rlwe

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"
	"math/bits"

	"github.com/google/go-cmp/cmp"
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

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = ring.GaloisGen

type DistributionLiteral interface{}

// ParametersLiteral is a literal representation of RLWE parameters. It has public fields and
// is used to express unchecked user-defined parameters literally into Go programs.
// The NewParametersFromLiteral function is used to generate the actual checked parameters
// from the literal representation.
//
// Users must set the polynomial degree (LogN) and the coefficient modulus, by either setting
// the Q and P fields to the desired moduli chain, or by setting the LogQ and LogP fields to
// the desired moduli sizes.
//
// Optionally, users may specify
// - the base 2 decomposition for the gadget ciphertexts
// - the error variance (Sigma) and secrets' density (H) and the ring type (RingType).
// If left unset, standard default values for these field are substituted at
// parameter creation (see NewParametersFromLiteral).
type ParametersLiteral struct {
	LogN           int
	Q              []uint64                    `json:",omitempty"`
	P              []uint64                    `json:",omitempty"`
	LogQ           []int                       `json:",omitempty"`
	LogP           []int                       `json:",omitempty"`
	Xe             ring.DistributionParameters `json:",omitempty"`
	Xs             ring.DistributionParameters `json:",omitempty"`
	RingType       ring.Type                   `json:",omitempty"`
	PlaintextScale Scale                       `json:",omitempty"`
	NTTFlag        bool                        `json:",omitempty"`
}

// Parameters represents a set of generic RLWE parameters. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	logN           int
	qi             []uint64
	pi             []uint64
	xe             distribution
	xs             distribution
	ringQ          *ring.Ring
	ringP          *ring.Ring
	ringType       ring.Type
	plaintextScale Scale
	nttFlag        bool
}

// NewParameters returns a new set of generic RLWE parameters from the given ring degree logn, moduli q and p, and
// error distribution Xs (secret) and Xe (error). It returns the empty parameters Parameters{} and a non-nil error if the
// specified parameters are invalid.
func NewParameters(logn int, q, p []uint64, xs, xe DistributionLiteral, ringType ring.Type, plaintextScale Scale, NTTFlag bool) (params Parameters, err error) {

	var lenP int
	if p != nil {
		lenP = len(p)
	}

	if err = checkSizeParams(logn, len(q), lenP); err != nil {
		return Parameters{}, err
	}

	params = Parameters{
		logN:           logn,
		qi:             make([]uint64, len(q)),
		pi:             make([]uint64, lenP),
		ringType:       ringType,
		plaintextScale: plaintextScale,
		nttFlag:        NTTFlag,
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

	if err = params.initRings(); err != nil {
		return Parameters{}, fmt.Errorf("cannot NewParameters: %w", err)
	}

	logQP := params.LogQP()

	switch xs := xs.(type) {
	case ring.Ternary, ring.DiscreteGaussian:
		params.xs = newDistribution(xs.(ring.DistributionParameters), logn, logQP)
	default:
		return Parameters{}, fmt.Errorf("secret distribution type must be Ternary or DiscretGaussian but is %T", xs)
	}
	if err != nil {
		return Parameters{}, err
	}

	switch xe := xe.(type) {
	case ring.Ternary, ring.DiscreteGaussian:
		params.xe = newDistribution(xe.(ring.DistributionParameters), logn, logQP)
	default:
		return Parameters{}, fmt.Errorf("error distribution type must be Ternary or DiscretGaussian but is %T", xe)
	}
	if err != nil {
		return Parameters{}, err
	}

	var warning error
	if params.XsHammingWeight() == 0 {
		warning = fmt.Errorf("warning secret standard HammingWeight is 0")
	}

	if params.xe.std <= 0 {
		if warning != nil {
			warning = fmt.Errorf("%w; warning error standard deviation 0", warning)
		} else {
			warning = fmt.Errorf("warning error standard deviation 0")
		}
	}

	return params, warning
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
// If the error variance is left unset, its value is set to `DefaultError`.
//
// If the RingType is left unset, the default value is ring.Standard.
func NewParametersFromLiteral(paramDef ParametersLiteral) (params Parameters, err error) {

	if paramDef.Xs == nil {
		paramDef.Xs = DefaultXs
	}

	if paramDef.Xe == nil {
		// prevents the zero value of ParameterLiteral to result in a noise-less parameter instance.
		// Users should use the NewParameters method to explicitely create noiseless instances.
		paramDef.Xe = DefaultXe
	}

	if paramDef.PlaintextScale.Cmp(Scale{}) == 0 {
		s := NewScale(1)
		paramDef.PlaintextScale = s
	}

	switch {
	case paramDef.Q != nil && paramDef.LogQ == nil:
		return NewParameters(paramDef.LogN, paramDef.Q, paramDef.P, paramDef.Xs, paramDef.Xe, paramDef.RingType, paramDef.PlaintextScale, paramDef.NTTFlag)
	case paramDef.LogQ != nil && paramDef.Q == nil:
		var q, p []uint64
		switch paramDef.RingType {
		case ring.Standard:
			q, p, err = GenModuli(paramDef.LogN+1, paramDef.LogQ, paramDef.LogP) //2NthRoot
		case ring.ConjugateInvariant:
			q, p, err = GenModuli(paramDef.LogN+2, paramDef.LogQ, paramDef.LogP) //4NthRoot
		default:
			return Parameters{}, fmt.Errorf("rlwe.NewParametersFromLiteral: invalid ring.Type, must be ring.ConjugateInvariant or ring.Standard")
		}
		if err != nil {
			return Parameters{}, err
		}
		return NewParameters(paramDef.LogN, q, p, paramDef.Xs, paramDef.Xe, paramDef.RingType, paramDef.PlaintextScale, paramDef.NTTFlag)
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
		Xe:             p.xe.params,
		Xs:             p.xs.params,
		RingType:       p.ringType,
		PlaintextScale: p.plaintextScale,
		NTTFlag:        p.nttFlag,
	}
}

// NewScale creates a new scale using the stored default scale as template.
func (p Parameters) NewScale(scale interface{}) Scale {
	newScale := NewScale(scale)
	newScale.Mod = p.plaintextScale.Mod
	return newScale
}

// N returns the ring degree
func (p Parameters) N() int {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p Parameters) LogN() int {
	return p.logN
}

// PlaintextDimensions returns the dimensions of the matrix that can be SIMD packed in a single plaintext polynomial.
// Returns [0, 0] by default.
func (p Parameters) PlaintextDimensions() [2]int {
	return [2]int{0, 0}
}

// PlaintextLogDimensions returns the log dimensions of the matrix that can be SIMD packed in a single plaintext polynomial.
// Returns [-1, -1] by default.
func (p Parameters) PlaintextLogDimensions() [2]int {
	return [2]int{-1, -1}
}

// PlaintextSlots returns the total number of entries (`slots`) that a plaintext can store.
// This value is obtained by multiplying all dimensions from PlaintextDimensions.
func (p Parameters) PlaintextSlots() int {
	dims := p.PlaintextDimensions()
	return dims[0] * dims[1]
}

// PlaintextLogSlots returns the total number of entries (`slots`) that a plaintext can store.
// This value is obtained by summing all log dimensions from PlaintextLogDimensions.
func (p Parameters) PlaintextLogSlots() int {
	dims := p.PlaintextLogDimensions()
	return dims[0] + dims[1]
}

// PlaintextScale returns the default scaling factor of the plaintext, if any.
func (p Parameters) PlaintextScale() Scale {
	return p.plaintextScale
}

// PlaintextModulus returns the plaintext modulus, if any. Else returns 0.
func (p Parameters) PlaintextModulus() uint64 {
	if p.plaintextScale.Mod != nil {
		return p.plaintextScale.Mod.Uint64()
	}

	return 0
}

// PlaintextPrecision returns the default precision in bits of the plaintext values which
// is max(53, log2(PlaintextScale)).
func (p Parameters) PlaintextPrecision() (prec uint) {
	if log2scale := math.Log2(p.PlaintextScale().Float64()); log2scale <= 53 {
		prec = 53
	} else {
		prec = uint(log2scale)
	}

	return
}

// PlaintextScaleToModuliRatio returns the default ratio between the scaling factor and moduli.
// This default ratio is computed as ceil(PlaintextScale/2^{60}).
// Returns 0 if the scaling factor is 0 (e.g. scale invariant scheme such as BFV).
func (p Parameters) PlaintextScaleToModuliRatio() int {

	if p.PlaintextScale().Mod != nil {
		return 1
	}

	scale := p.PlaintextScale().Float64()
	nbModuli := 0
	for scale > 1 {
		scale /= 0xfffffffffffffff
		nbModuli++
	}
	return nbModuli
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

// MaxDepth returns MaxLevel / PlaintextScaleToModuliRatio which is the maximum number of multiplicaitons
// followed by a rescaling that can be carried out with on a ciphertext with the plaintextScale.
// Returns 0 if the scaling factor is zero (e.g. scale invariant scheme such as BFV).
func (p Parameters) MaxDepth() int {
	if ratio := p.PlaintextScaleToModuliRatio(); ratio > 0 {
		return p.MaxLevel() / ratio
	}
	return 0
}

// NTTFlag returns a boolean indicating if elements are stored by default in the NTT domain.
func (p Parameters) NTTFlag() bool {
	return p.nttFlag
}

// Xs returns the Distribution of the secret
func (p Parameters) Xs() ring.DistributionParameters {
	return p.xs.params
}

// XsHammingWeight returns the expected Hamming weight of the secret.
func (p Parameters) XsHammingWeight() int {
	switch xs := p.xs.params.(type) {
	case ring.Ternary:
		if xs.H != 0 {
			return xs.H
		} else {
			return int(math.Ceil(float64(p.N()) * xs.P))
		}
	case ring.DiscreteGaussian:
		return int(math.Ceil(float64(p.N()) * float64(xs.Sigma) * math.Sqrt(2.0/math.Pi)))
	default:
		panic(fmt.Sprintf("invalid error distribution: must be DiscretGaussian, Ternary but is %T", xs))
	}
}

// Xe returns Distribution of the error
func (p Parameters) Xe() ring.DistributionParameters {
	return p.xe.params
}

// NoiseBound returns truncation bound for the error distribution.
func (p Parameters) NoiseBound() float64 {
	return p.xe.absBound
}

// NoiseFreshPK returns the standard deviation
// of a fresh encryption with the public key.
func (p Parameters) NoiseFreshPK() (std float64) {

	std = float64(p.XsHammingWeight() + 1)

	if p.RingP() != nil {
		std *= 1 / 12.0
	} else {
		sigma := float64(p.xe.std)
		std *= sigma * sigma
	}

	return math.Sqrt(std)
}

// NoiseFreshSK returns the standard deviation
// of a fresh encryption with the secret key.
func (p Parameters) NoiseFreshSK() (std float64) {
	return float64(p.xe.std)
}

// RingType returns the type of the underlying ring.
func (p Parameters) RingType() ring.Type {
	return p.ringType
}

// MaxLevel returns the maximum level of a ciphertext.
func (p Parameters) MaxLevel() int {
	return p.MaxLevelQ()
}

// MaxLevelQ returns the maximum level of the modulus Q.
func (p Parameters) MaxLevelQ() int {
	return p.QCount() - 1
}

// MaxLevelP returns the maximum level of the modulus P.
func (p Parameters) MaxLevelP() int {
	return p.PCount() - 1
}

// Q returns a new slice with the factors of the ciphertext modulus q
func (p Parameters) Q() []uint64 {
	qi := make([]uint64, len(p.qi))
	copy(qi, p.qi)
	return qi
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

// PBigInt return the ciphertext-space extension modulus P in big.Integer, reconstructed, representation.
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
func (p Parameters) LogQ() (logq float64) {
	return p.ringQ.LogModuli()
}

// LogP returns the size of the extended modulus P in bits
func (p Parameters) LogP() (logp float64) {
	if p.ringP == nil {
		return 0
	}
	return p.ringP.LogModuli()
}

// LogQP returns the size of the extended modulus QP in bits
func (p Parameters) LogQP() (logqp float64) {
	return p.LogQ() + p.LogP()
}

// MaxBit returns max(max(bitLen(Q[:levelQ+1])), max(bitLen(P[:levelP+1])).
func (p Parameters) MaxBit(levelQ, levelP int) (c int) {
	for _, qi := range p.Q()[:levelQ+1] {
		c = utils.Max(c, bits.Len64(qi))
	}

	for _, pi := range p.P()[:levelP+1] {
		c = utils.Max(c, bits.Len64(pi))
	}
	return
}

// DecompPw2 returns ceil(p.MaxBitQ(levelQ, levelP)/Base2Decomposition).
func (p Parameters) DecompPw2(levelQ, levelP, Base2Decomposition int) (c int) {
	if Base2Decomposition == 0 || levelP > 0 {
		return 1
	}

	return (p.MaxBit(levelQ, levelP) + Base2Decomposition - 1) / Base2Decomposition
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
func (p Parameters) QiOverflowMargin(level int) int {
	return int(math.Exp2(64) / float64(utils.MaxSlice(p.qi[:level+1])))
}

// PiOverflowMargin returns floor(2^64 / max(Pi)), i.e. the number of times elements of Z_max{Pi} can
// be added together before overflowing 2^64.
func (p Parameters) PiOverflowMargin(level int) int {
	return int(math.Exp2(64) / float64(utils.MaxSlice(p.pi[:level+1])))
}

// GaloisElements takes a list of integers k and returns the list [GaloisGen^{k[i]} mod NthRoot, ...].
func (p Parameters) GaloisElements(k []int) (galEls []uint64) {
	galEls = make([]uint64, len(k))
	for i, ki := range k {
		galEls[i] = p.GaloisElement(ki)
	}
	return
}

// GaloisElement takes an integer k and returns GaloisGen^{k} mod NthRoot.
func (p Parameters) GaloisElement(k int) uint64 {
	return ring.ModExp(GaloisGen, uint64(k)&(p.ringQ.NthRoot()-1), p.ringQ.NthRoot())
}

// ModInvGaloisElement takes a Galois element of the form GaloisGen^{k} mod NthRoot
// and returns GaloisGen^{-k} mod NthRoot.
func (p Parameters) ModInvGaloisElement(galEl uint64) uint64 {
	return ring.ModExp(galEl, p.ringQ.NthRoot()-1, p.ringQ.NthRoot())
}

// GaloisElementInverse returns GaloisGen^{-1} mod NthRoot
func (p Parameters) GaloisElementInverse() uint64 {
	if p.ringType == ring.ConjugateInvariant {
		panic("Cannot generate GaloisElementInverse if ringType is ConjugateInvariant")
	}
	return p.ringQ.NthRoot() - 1
}

// GaloisElementsForTrace returns the list of Galois elements requored for the for the `Trace` operation.
// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 2^{LogN} <= i < N.
func (p Parameters) GaloisElementsForTrace(logN int) (galEls []uint64) {

	galEls = []uint64{}
	for i, j := logN, 0; i < p.LogN()-1; i, j = i+1, j+1 {
		galEls = append(galEls, p.GaloisElement(1<<i))
	}

	if logN == 0 {
		switch p.ringType {
		case ring.Standard:
			galEls = append(galEls, p.GaloisElementInverse())
		case ring.ConjugateInvariant:
			panic("cannot GaloisElementsForTrace: Galois element GaloisGen^-1 is undefined in ConjugateInvariant Ring")
		default:
			panic("cannot GaloisElementsForTrace: invalid ring type")
		}
	}

	return
}

// GaloisElementsForReplicate returns the list of Galois elements necessary to perform the
// `Replicate` operation with parameters `batch` and `n`.
func (p Parameters) GaloisElementsForReplicate(batch, n int) (galEls []uint64) {
	return p.GaloisElementsForInnerSum(-batch, n)
}

// GaloisElementsForInnerSum returns the list of Galois elements necessary to apply the method
// `InnerSum` operation with parameters `batch` and `n`.
func (p Parameters) GaloisElementsForInnerSum(batch, n int) (galEls []uint64) {

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

	rotations := make([]int, len(rotIndex))
	var i int
	for j := range rotIndex {
		rotations[i] = j
		i++
	}

	return p.GaloisElements(rotations)
}

// GaloisElementsForExpand returns the list of Galois elements required
// to perform the `Expand` operation with parameter `logN`.
func (p Parameters) GaloisElementsForExpand(logN int) (galEls []uint64) {
	galEls = make([]uint64, logN)

	NthRoot := p.RingQ().NthRoot()

	for i := 0; i < logN; i++ {
		galEls[i] = uint64(NthRoot/(2<<i) + 1)
	}

	return
}

// GaloisElementsForPack returns the list of Galois elements required
// to perform the `Merge` operation.
func (p Parameters) GaloisElementsForPack(logGap int) (galEls []uint64) {

	if logGap > p.logN || logGap < 0 {
		panic("cannot GaloisElementsForPack: logGap > logN || logGap < 0")
	}

	galEls = make([]uint64, 0, logGap)
	for i := 0; i < logGap; i++ {
		galEls = append(galEls, p.GaloisElement(1<<i))
	}

	switch p.ringType {
	case ring.Standard:
		if logGap == p.logN {
			galEls = append(galEls, p.GaloisElementInverse())
		}
	default:
		panic("cannot GaloisElementsForPack: invalid ring type")
	}
	return
}

// GaloisElementsForLinearTransform returns the list of Galois elements required to perform a linear transform
// with the provided non-zero diagonales.
// Set LogBSGSRatio < 0 to return the Galois elements for a naive evaluation of the linear transform.
func (p Parameters) GaloisElementsForLinearTransform(nonZeroDiagonals []int, LogSlots, LogBSGSRatio int) (galEls []uint64) {

	slots := 1 << LogSlots

	if LogBSGSRatio < 0 {

		_, _, rotN2 := BSGSIndex(nonZeroDiagonals, slots, slots)

		galEls = make([]uint64, len(rotN2))

		for i := range rotN2 {
			galEls[i] = p.GaloisElement(rotN2[i])
		}

		return
	}

	N1 := FindBestBSGSRatio(nonZeroDiagonals, slots, LogBSGSRatio)

	_, rotN1, rotN2 := BSGSIndex(nonZeroDiagonals, slots, N1)

	return p.GaloisElements(utils.GetDistincts(append(rotN1, rotN2...)))
}

// SolveDiscreteLogGaloisElement takes a Galois element of the form GaloisGen^{k} mod NthRoot and returns k.
func (p Parameters) SolveDiscreteLogGaloisElement(galEl uint64) (k int) {

	N := p.ringQ.NthRoot()

	var kuint uint64

	x := N >> 3

	for {

		if ring.ModExpPow2(GaloisGen, kuint, N) != ring.ModExpPow2(galEl, x, N) {
			kuint |= N >> 3
		}

		if x == 1 {
			return int(kuint)
		}

		x >>= 1
		kuint >>= 1
	}
}

// Equal checks two Parameter structs for equality.
func (p Parameters) Equal(other ParametersInterface) (res bool) {

	switch other := other.(type) {
	case Parameters:
		res = p.logN == other.logN
		res = res && (p.xs.params == other.xs.params)
		res = res && (p.xe.params == other.xe.params)
		res = res && cmp.Equal(p.qi, other.qi)
		res = res && cmp.Equal(p.pi, other.pi)
		res = res && (p.ringType == other.ringType)
		res = res && (p.plaintextScale.Equal(other.plaintextScale))
		res = res && (p.nttFlag == other.nttFlag)
		return
	}

	return false
}

// MarshalBinary returns a []byte representation of the parameter set.
// This representation corresponds to the MarshalJSON representation.
func (p Parameters) MarshalBinary() ([]byte, error) {
	return p.MarshalJSON()
}

// UnmarshalBinary decodes a slice of bytes on the target Parameters.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	return p.UnmarshalJSON(data)
}

// MarshalJSON returns a JSON representation of this parameter set. See `Marshal` from the `encoding/json` package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	return json.Marshal(p.ParametersLiteral())
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
func GenModuli(LogNthRoot int, logQ, logP []int) (q, p []uint64, err error) {

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
	for bitsize, value := range primesbitlen {

		g := ring.NewNTTFriendlyPrimesGenerator(uint64(bitsize), uint64(1<<LogNthRoot))

		if bitsize == 61 {
			if primes[bitsize], err = g.NextDownstreamPrimes(value); err != nil {
				return q, p, fmt.Errorf("cannot GenModuli: failed to generate %d primes of bit-size=61 for LogNthRoot=%d: %w", value, LogNthRoot, err)
			}
		} else {
			if primes[bitsize], err = g.NextAlternatingPrimes(value); err != nil {
				return q, p, fmt.Errorf("cannot GenModuli: failed to generate %d primes of bit-size=%d for LogNthRoot=%d: %w", value, bitsize, LogNthRoot, err)
			}
		}
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
		return fmt.Errorf("initRings/ringQ: %w", err)
	}
	if len(p.pi) != 0 {
		if p.ringP, err = ring.NewRingFromType(1<<p.logN, p.pi, p.ringType); err != nil {
			return fmt.Errorf("initRings/ringP: %w", err)
		}
	}
	return
}

func (p *ParametersLiteral) UnmarshalJSON(b []byte) (err error) {
	var pl struct {
		LogN           int
		Q              []uint64
		P              []uint64
		LogQ           []int
		LogP           []int
		Xe             map[string]interface{}
		Xs             map[string]interface{}
		RingType       ring.Type
		PlaintextScale Scale
		NTTFlag        bool
	}

	err = json.Unmarshal(b, &pl)
	if err != nil {
		return err
	}

	p.LogN = pl.LogN
	p.Q, p.P, p.LogQ, p.LogP = pl.Q, pl.P, pl.LogQ, pl.LogP
	if pl.Xs != nil {
		p.Xs, err = ring.ParametersFromMap(pl.Xs)
		if err != nil {
			return err
		}
	}
	if pl.Xe != nil {
		p.Xe, err = ring.ParametersFromMap(pl.Xe)
		if err != nil {
			return err
		}
	}
	p.RingType = pl.RingType
	p.PlaintextScale = pl.PlaintextScale
	p.NTTFlag = pl.NTTFlag

	return err
}

// LatticeEstimatorSageMathCell returns a string formated SageMath cell of the code
// to run using the Lattice estimator (https://github.com/malb/lattice-estimator)
// to estimate the security of the target Parameters.
func LatticeEstimatorSageMathCell(p Parameters) string {

	LogN := p.LogN()
	LogQP := p.LogQP()
	Xs := p.xs
	Xe := p.xe

	return fmt.Sprintf(`# 1) Clone https://github.com/malb/lattice-estimator
# 2) Create a new SageMath notebook in the folder 
# 3) Copy-past the following code in a new cell
# ================================================================
from estimator import *
from estimator.nd import NoiseDistribution
from estimator import LWE

n = 1<<%d
q = 1<<%d
Xs = NoiseDistribution.(stddev=%f, mean=0, n=n, bounds=(%f, %f), density=%f, tag=%s)
Xe = NoiseDistribution.(stddev=%f, mean=0, n=n, bounds=(%f, %f), density=%f, tag=%s)

params = LWE.Parameters(n=n, q=q, Xs=Xs, Xe=Xe)

print(params)

LWE.estimate(params)
`,
		LogN,
		int(math.Round(LogQP)),
		Xs.std, Xs.bounds[0], Xs.bounds[1], Xs.density, Xs.params.Type(),
		Xe.std, Xe.bounds[0], Xe.bounds[1], Xe.density, Xe.params.Type())
}
