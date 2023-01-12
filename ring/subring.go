package ring

import (
	"encoding/binary"
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// SubRing is a struct storing precomputation
// for fast modular reduction and NTT for
// a given modulus.
type SubRing struct {
	ntt NumberTheoreticTransformer

	// Polynomial nb.Coefficients
	N int

	// Nthroot used for the NTT
	NthRoot uint64

	// Modulus
	Modulus uint64

	// Unique factors of Modulus-1
	Factors []uint64

	// 2^bit_length(Modulus) - 1
	Mask uint64

	// Fast reduction constants
	BRedConstant []uint64 // Barrett Reduction
	MRedConstant uint64   // Montgomery Reduction

	AllowsNTT bool // Indicates whether NTT can be used with the current ring.

	PrimitiveRoot uint64 // 2N-th primitive root

	RootsForward  []uint64 //powers of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	RootsBackward []uint64 //powers of the inverse of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	NInv          uint64   //[N^-1] mod Modulus in Montgomery form
}

// NewSubRing creates a new SubRing with the standard NTT.
// NTT constants still need to be generated using .GenNTTConstants(NthRoot uint64).
func NewSubRing(N int, Modulus uint64) (s *SubRing, err error) {
	return NewSubRingWithCustomNTT(N, Modulus, NewNumberTheoreticTransformerStandard, 2*N)
}

// NewSubRingWithCustomNTT creates a new SubRing with degree N and modulus Modulus with user-defined NTT transform and primitive Nth root of unity.
// Modulus should be equal to 1 modulo the root of unity.
// N must be a power of two larger than 8. An error is returned with a nil *SubRing in the case of non NTT-enabling parameters.
func NewSubRingWithCustomNTT(N int, Modulus uint64, ntt func(*SubRing, int) NumberTheoreticTransformer, NthRoot int) (s *SubRing, err error) {

	// Checks if N is a power of 2
	if (N < 16) || (N&(N-1)) != 0 && N != 0 {
		return nil, fmt.Errorf("invalid degree (must be a power of 2 >= 8)")
	}

	s = &SubRing{}

	s.N = N

	s.NthRoot = uint64(NthRoot)

	s.AllowsNTT = false

	s.Modulus = Modulus
	s.Mask = (1 << uint64(bits.Len64(Modulus-1))) - 1

	// Computes the fast modular reduction constants for the Ring
	s.BRedConstant = BRedConstant(Modulus)

	// If qi is not a power of 2, we can compute the MRed (otherwise, it
	// would return an error as there is no valid Montgomery form mod a power of 2)
	if (Modulus&(Modulus-1)) != 0 && Modulus != 0 {
		s.MRedConstant = MRedConstant(Modulus)
	}

	s.RootsForward = make([]uint64, NthRoot>>1)
	s.RootsBackward = make([]uint64, NthRoot>>1)

	s.ntt = ntt(s, N)

	return
}

// Type returns the Type of subring which might be either `Standard` or `ConjugateInvariant`.
func (s *SubRing) Type() Type {
	switch s.ntt.(type) {
	case NumberTheoreticTransformerStandard:
		return Standard
	case NumberTheoreticTransformerConjugateInvariant:
		return ConjugateInvariant
	default:
		panic(fmt.Errorf("invalid NumberTheoreticTransformer type: %T", s.ntt))
	}
}

// generateNTTConstants generates the NTT constant for the target SubRing.
// The fields `PrimitiveRoot` and `Factors` can be set manually to
// bypasse the search for the primitive root (which requires to
// factor Modulus-1) and speedup the generation of the constants.
func (s *SubRing) generateNTTConstants() (err error) {

	if s.N == 0 || s.Modulus == 0 {
		return fmt.Errorf("invalid t parameters (missing)")
	}

	Modulus := s.Modulus
	NthRoot := s.NthRoot

	// Checks if each qi is prime and equal to 1 mod NthRoot
	if !IsPrime(Modulus) {
		return fmt.Errorf("invalid modulus: %d is not prime)", Modulus)
	}

	if Modulus&(NthRoot-1) != 1 {
		return fmt.Errorf("invalid modulus: %d != 1 mod NthRoot)", Modulus)
	}

	// It is possible to manually set the primitive root along with the factors of q-1.
	// This is notably useful when marhsalling the SubRing, to avoid re-factoring q-1.
	// If both are set, then checks that that the root is indeed primitive.
	// Else, factorize q-1 and finds a primitive roos.
	if s.PrimitiveRoot != 0 && s.Factors != nil {
		if err = CheckPrimitiveRoot(s.PrimitiveRoot, s.Modulus, s.Factors); err != nil {
			return
		}
	} else {
		if s.PrimitiveRoot, s.Factors, err = PrimitiveRoot(Modulus, s.Factors); err != nil {
			return
		}
	}

	logNthRoot := uint64(bits.Len64(NthRoot>>1) - 1)

	// 1.1 Computes N^(-1) mod Q in Montgomery form
	s.NInv = MForm(ModExp(NthRoot>>1, Modulus-2, Modulus), Modulus, s.BRedConstant)

	// 1.2 Computes Psi and PsiInv in Montgomery form

	// Computes Psi and PsiInv in Montgomery form
	PsiMont := MForm(ModExp(s.PrimitiveRoot, (Modulus-1)/NthRoot, Modulus), Modulus, s.BRedConstant)
	PsiInvMont := MForm(ModExp(s.PrimitiveRoot, Modulus-((Modulus-1)/NthRoot)-1, Modulus), Modulus, s.BRedConstant)

	s.RootsForward[0] = MForm(1, Modulus, s.BRedConstant)
	s.RootsBackward[0] = MForm(1, Modulus, s.BRedConstant)

	// Computes nttPsi[j] = nttPsi[j-1]*Psi and RootsBackward[j] = RootsBackward[j-1]*PsiInv
	for j := uint64(1); j < NthRoot>>1; j++ {

		indexReversePrev := utils.BitReverse64(uint64(j-1), logNthRoot)
		indexReverseNext := utils.BitReverse64(uint64(j), logNthRoot)

		s.RootsForward[indexReverseNext] = MRed(s.RootsForward[indexReversePrev], PsiMont, Modulus, s.MRedConstant)
		s.RootsBackward[indexReverseNext] = MRed(s.RootsBackward[indexReversePrev], PsiInvMont, Modulus, s.MRedConstant)
	}

	s.AllowsNTT = true

	return
}

// PrimitiveRoot computes the smallest primitive root of the given prime q
// The unique factors of q-1 can be given to speed up the search for the roos.
func PrimitiveRoot(q uint64, factors []uint64) (uint64, []uint64, error) {

	if factors != nil {
		if err := CheckFactors(q-1, factors); err != nil {
			return 0, factors, err
		}
	} else {

		factorsBig := utils.GetFactors(new(big.Int).SetUint64(q - 1)) //Factor q-1, might be slow

		factors = make([]uint64, len(factorsBig))
		for i := range factors {
			factors[i] = factorsBig[i].Uint64()
		}
	}

	notFoundPrimitiveRoot := true

	var g uint64 = 2

	for notFoundPrimitiveRoot {
		g++
		for _, factor := range factors {
			// if for any factor of q-1, g^(q-1)/factor = 1 mod q, g is not a primitive root
			if ModExp(g, (q-1)/factor, q) == 1 {
				notFoundPrimitiveRoot = true
				break
			}
			notFoundPrimitiveRoot = false
		}
	}

	return g, factors, nil
}

// CheckFactors checks that the given list of factors contains
// all the unique primes of m.
func CheckFactors(m uint64, factors []uint64) (err error) {

	for _, factor := range factors {

		if !IsPrime(factor) {
			return fmt.Errorf("composite factor")
		}

		for m%factor == 0 {
			m /= factor
		}
	}

	if m != 1 {
		return fmt.Errorf("incomplete factor list")
	}

	return
}

// CheckPrimitiveRoot checks that g is a valid primtive root mod q,
// given the factors of q-1.
func CheckPrimitiveRoot(g, q uint64, factors []uint64) (err error) {

	if err = CheckFactors(q-1, factors); err != nil {
		return
	}

	for _, factor := range factors {
		if ModExp(g, (q-1)/factor, q) == 1 {
			return fmt.Errorf("invalid primitive root")
		}
	}

	return
}

// MarshalBinarySize returns the length in bytes of the target SubRing.
func (s *SubRing) MarshalBinarySize() (dataLen int) {
	dataLen++                     // RingType
	dataLen++                     // LogN
	dataLen++                     // NthRoot
	dataLen += 8                  // Modulus
	dataLen++                     // #Factors
	dataLen += len(s.Factors) * 8 // Factors
	dataLen += 8                  // PrimitiveRoot
	return
}

// Encode encodes the target SubRing on a slice of bytes and returns
// the number of bytes written.
func (s *SubRing) Encode(data []byte) (ptr int, err error) {
	data[ptr] = uint8(s.Type())
	ptr++
	data[ptr] = uint8(bits.Len64(uint64(s.N - 1)))
	ptr++
	data[ptr] = uint8(int(s.NthRoot) / s.N)
	ptr++
	binary.LittleEndian.PutUint64(data[ptr:], s.Modulus)
	ptr += 8
	data[ptr] = uint8(len(s.Factors))
	ptr++
	for i := range s.Factors {
		binary.LittleEndian.PutUint64(data[ptr:], s.Factors[i])
		ptr += 8
	}
	binary.LittleEndian.PutUint64(data[ptr:], s.PrimitiveRoot)
	ptr += 8
	return
}

// Decode decodes the input slice of bytes on the target SubRing and
// returns the number of bytes read.
func (s *SubRing) Decode(data []byte) (ptr int, err error) {
	ringType := Type(data[ptr])
	ptr++
	s.N = 1 << int(data[ptr])
	ptr++
	s.NthRoot = uint64(s.N) * uint64(data[ptr])
	ptr++
	s.Modulus = binary.LittleEndian.Uint64(data[ptr:])
	ptr += 8
	s.Factors = make([]uint64, data[ptr])
	ptr++
	for i := range s.Factors {
		s.Factors[i] = binary.LittleEndian.Uint64(data[ptr:])
		ptr += 8
	}
	s.PrimitiveRoot = binary.LittleEndian.Uint64(data[ptr:])
	ptr += 8

	s.Mask = (1 << uint64(bits.Len64(s.Modulus-1))) - 1

	// Computes the fast modular reduction parameters for the Ring
	s.BRedConstant = BRedConstant(s.Modulus)

	// If qi is not a power of 2, we can compute the MRed (otherwise, it
	// would return an error as there is no valid Montgomery form mod a power of 2)
	if (s.Modulus&(s.Modulus-1)) != 0 && s.Modulus != 0 {
		s.MRedConstant = MRedConstant(s.Modulus)
	}

	s.RootsForward = make([]uint64, s.NthRoot>>1)
	s.RootsBackward = make([]uint64, s.NthRoot>>1)

	switch ringType {
	case Standard:

		s.ntt = NewNumberTheoreticTransformerStandard(s, s.N)

		if int(s.NthRoot) < s.N<<1 {
			return ptr, fmt.Errorf("invalid ring type: NthRoot must be at least 2N but is %dN", int(s.NthRoot)/s.N)
		}

	case ConjugateInvariant:

		s.ntt = NewNumberTheoreticTransformerConjugateInvariant(s, s.N)

		if int(s.NthRoot) < s.N<<2 {
			return ptr, fmt.Errorf("invalid ring type: NthRoot must be at least 4N but is %dN", int(s.NthRoot)/s.N)
		}

	default:
		return ptr, fmt.Errorf("invalid ring type")
	}

	if err = s.generateNTTConstants(); err != nil {
		return
	}

	return
}
