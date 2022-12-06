package ring

import (
	"encoding/binary"
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// Table is a struct storing precomputation
// for fast modular reduction and NTT for
// a given modulus.
type Table struct {

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

	// Fast reduction parameters
	BRedParams []uint64
	MRedParams uint64

	AllowsNTT bool // Indicates whether NTT can be used with the current ring.

	PrimitiveRoot uint64 // 2N-th primitive root

	RootsForward  []uint64 //powers of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	RootsBackward []uint64 //powers of the inverse of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	NInv          uint64   //[N^-1] mod Modulus in Montgomery form
}

// NewTable generates a new table from a ring degree N and a modulus.
func NewTable(N int, Modulus uint64) (t *Table) {

	t = &Table{}

	t.N = N

	t.AllowsNTT = false

	t.Modulus = Modulus
	t.Mask = (1 << uint64(bits.Len64(Modulus-1))) - 1

	// Computes the fast modular reduction parameters for the Ring
	t.BRedParams = BRedParams(Modulus)

	// If qi is not a power of 2, we can compute the MRedParams (otherwise, it
	// would return an error as there is no valid Montgomery form mod a power of 2)
	if (Modulus&(Modulus-1)) != 0 && Modulus != 0 {
		t.MRedParams = MRedParams(Modulus)
	}

	return
}

// GenNTTParams generates the NTT tables.
// The fields `PrimitiveRoot` and `Factors` can be set manually to
// bypasse the search for the primitive root (which requires to
// factor Modulus-1) and speedup the generation of the Table.
func (t *Table) GenNTTParams(NthRoot uint64) (err error) {

	if t.N == 0 || t.Modulus == 0 || NthRoot < 1 {
		return fmt.Errorf("invalid t parameters (missing)")
	}

	Modulus := t.Modulus

	// Checks if each qi is prime and equal to 1 mod NthRoot
	if !IsPrime(Modulus) {
		return fmt.Errorf("invalid modulus: %d is not prime)", Modulus)
	}

	if Modulus&(NthRoot-1) != 1 {
		return fmt.Errorf("invalid modulus: %d != 1 mod NthRoot)", Modulus)
	}

	t.NthRoot = NthRoot

	// It is possible to manually set the primitive root along with the factors of q-1.
	// This is notably useful when marhsalling the table, to avoid re-factoring q-1.
	// If both are set, then checks that that the root is indeed primitive.
	// Else, factorize q-1 and finds a primitive root.
	if t.PrimitiveRoot != 0 && t.Factors != nil {
		if err = CheckPrimitiveRoot(t.PrimitiveRoot, t.Modulus, t.Factors); err != nil {
			return
		}
	} else {
		if t.PrimitiveRoot, t.Factors, err = PrimitiveRoot(Modulus, t.Factors); err != nil {
			return
		}
	}

	logNthRoot := uint64(bits.Len64(NthRoot>>1) - 1)

	// 1.1 Computes N^(-1) mod Q in Montgomery form
	t.NInv = MForm(ModExp(NthRoot>>1, Modulus-2, Modulus), Modulus, t.BRedParams)

	// 1.2 Computes Psi and PsiInv in Montgomery form
	t.RootsForward = make([]uint64, NthRoot>>1)
	t.RootsBackward = make([]uint64, NthRoot>>1)

	// Computes Psi and PsiInv in Montgomery form
	PsiMont := MForm(ModExp(t.PrimitiveRoot, (Modulus-1)/NthRoot, Modulus), Modulus, t.BRedParams)
	PsiInvMont := MForm(ModExp(t.PrimitiveRoot, Modulus-((Modulus-1)/NthRoot)-1, Modulus), Modulus, t.BRedParams)

	t.RootsForward[0] = MForm(1, Modulus, t.BRedParams)
	t.RootsBackward[0] = MForm(1, Modulus, t.BRedParams)

	// Computes nttPsi[j] = nttPsi[j-1]*Psi and RootsBackward[j] = RootsBackward[j-1]*PsiInv
	for j := uint64(1); j < NthRoot>>1; j++ {

		indexReversePrev := utils.BitReverse64(uint64(j-1), logNthRoot)
		indexReverseNext := utils.BitReverse64(uint64(j), logNthRoot)

		t.RootsForward[indexReverseNext] = MRed(t.RootsForward[indexReversePrev], PsiMont, Modulus, t.MRedParams)
		t.RootsBackward[indexReverseNext] = MRed(t.RootsBackward[indexReversePrev], PsiInvMont, Modulus, t.MRedParams)
	}

	t.AllowsNTT = true

	return
}

// PrimitiveRoot computes the smallest primitive root of the given prime q
// The unique factors of q-1 can be given to speed up the search for the root.
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

// MarshalBinarySize returns the length in bytes of the target Table.
func (t *Table) MarshalBinarySize() int {
	return 1 + 1 + 8 + 1 + len(t.Factors)*8 + 8
}

// Encode encodes the target table on a slice of bytes and returns
// the number of bytes written.
func (t *Table) Encode(data []byte) (ptr int, err error) {
	data[ptr] = uint8(bits.Len64(uint64(t.N - 1)))
	ptr++
	data[ptr] = uint8(int(t.NthRoot) / t.N)
	ptr++
	binary.LittleEndian.PutUint64(data[ptr:], t.Modulus)
	ptr += 8
	data[ptr] = uint8(len(t.Factors))
	ptr++
	for i := range t.Factors {
		binary.LittleEndian.PutUint64(data[ptr:], t.Factors[i])
		ptr += 8
	}
	binary.LittleEndian.PutUint64(data[ptr:], t.PrimitiveRoot)
	ptr += 8
	return
}

// Decode decodes the input slice of bytes on the target Table and
// returns the number of bytes read.
func (t *Table) Decode(data []byte) (ptr int, err error) {
	t.N = 1 << int(data[ptr])
	ptr++
	t.NthRoot = uint64(t.N) * uint64(data[ptr])
	ptr++
	t.Modulus = binary.LittleEndian.Uint64(data[ptr:])
	ptr += 8
	t.Factors = make([]uint64, data[ptr])
	ptr++
	for i := range t.Factors {
		t.Factors[i] = binary.LittleEndian.Uint64(data[ptr:])
		ptr += 8
	}
	t.PrimitiveRoot = binary.LittleEndian.Uint64(data[ptr:])
	ptr += 8

	t.Mask = (1 << uint64(bits.Len64(t.Modulus-1))) - 1

	// Computes the fast modular reduction parameters for the Ring
	t.BRedParams = BRedParams(t.Modulus)

	// If qi is not a power of 2, we can compute the MRedParams (otherwise, it
	// would return an error as there is no valid Montgomery form mod a power of 2)
	if (t.Modulus&(t.Modulus-1)) != 0 && t.Modulus != 0 {
		t.MRedParams = MRedParams(t.Modulus)
	}

	if err = t.GenNTTParams(t.NthRoot); err != nil {
		return
	}

	return
}
