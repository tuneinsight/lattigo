package rlwe

import (
	"math"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
)

// MaxLogN is the log2 of the largest supported polynomial modulus degree.
const MaxLogN = 16

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
const GaloisGen uint64 = 5

type Parameters struct {
	logN  int
	qi    []uint64
	pi    []uint64
	sigma float64
}

func NewRLWEParameters(logn int, q, p []uint64, sigma float64) *Parameters { // TEMPORARY constructor
	return &Parameters{
		logN:  logn,
		pi:    p,
		qi:    q,
		sigma: sigma,
	}
}

// N returns the ring degree
func (p *Parameters) N() int {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) LogN() int {
	return p.logN
}

// Sigma returns standard deviation of the noise distribution
func (p *Parameters) Sigma() float64 {
	return p.sigma
}

// Qi returns a new slice with the factors of the ciphertext modulus q
func (p *Parameters) Q() []uint64 {
	qi := make([]uint64, len(p.qi))
	copy(qi, p.qi)
	return qi
}

// QiCount returns the number of factors of the ciphertext modulus Q
func (p *Parameters) QiCount() int {
	return len(p.qi)
}

func (p *Parameters) QBigInt() *big.Int {
	q := big.NewInt(1)
	for _, qi := range p.qi {
		q.Mul(q, new(big.Int).SetUint64(qi))
	}
	return q
}

// Pi returns a new slice with the factors of the ciphertext modulus extension P
func (p *Parameters) P() []uint64 {
	pi := make([]uint64, len(p.pi))
	copy(pi, p.pi)
	return pi
}

// PiCount returns the number of factors of the ciphertext modulus extension P
func (p *Parameters) PiCount() int {
	return len(p.pi)
}

func (p *Parameters) PBigInt() *big.Int {
	pInt := big.NewInt(1)
	for _, pi := range p.pi {
		pInt.Mul(pInt, new(big.Int).SetUint64(pi))
	}
	return pInt
}

// Pi returns a new slice with the factors of the ciphertext modulus extension P
func (p *Parameters) QP() []uint64 {
	qp := make([]uint64, len(p.qi)+len(p.pi))
	copy(qp, p.qi)
	copy(qp[len(p.qi):], p.pi)
	return qp
}

// QPiCount returns the number of factors of the ciphertext modulus + the modulus extension P
func (p *Parameters) QPCount() int {
	return len(p.qi) + len(p.pi)
}

func (p *Parameters) QPBigInt() *big.Int {
	pqInt := p.QBigInt()
	pqInt.Mul(pqInt, p.PBigInt())
	return pqInt
}

// LogQP returns the size of the extended modulus QP in bits
func (p *Parameters) LogQP() int {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return tmp.BitLen()
}

// Alpha returns the number of moduli in in P
func (p *Parameters) Alpha() int {
	return p.PiCount()
}

// Beta returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (p *Parameters) Beta() int {
	if p.Alpha() != 0 {
		return int(math.Ceil(float64(p.QiCount()) / float64(p.Alpha())))
	}

	return 1
}

func (p *Parameters) RingQ() *ring.Ring {
	ringQ, err := ring.NewRing(p.N(), p.qi)
	if err != nil {
		panic(err) // Parameter type invariant
	}
	return ringQ
}

func (p *Parameters) RingP() *ring.Ring {
	if len(p.pi) == 0 {
		return nil
	}
	ringP, err := ring.NewRing(p.N(), p.pi)
	if err != nil {
		panic(err) // Parameter type invariant
	}
	return ringP
}

func (p *Parameters) RingQP() *ring.Ring {
	ringQP, err := ring.NewRing(p.N(), append(p.qi, p.pi...))
	if err != nil {
		panic(err) // Parameter type invariant
	}
	return ringQP
}

// GaloisElementForColumnRotationBy returns the galois element for plaintext
// column rotations by k position to the left. Providing a negative k is
// equivalent to a right rotation.
func (p *Parameters) GaloisElementForColumnRotationBy(k int) uint64 {
	twoN := 1 << (p.logN + 1)
	mask := twoN - 1
	kRed := k & mask
	return ring.ModExp(GaloisGen, kRed, uint64(twoN))
}

// GaloisElementForRowRotation returns the galois element for generating the row
// rotation automorphism
func (p *Parameters) GaloisElementForRowRotation() uint64 {
	return (1 << (p.logN + 1)) - 1
}

// GaloisElementsForRowInnerSum returns a list of all galois elements required to
// perform an InnerSum operation. This corresponds to all the left rotations by
// k-positions where k is a power of two and the row-rotation element.
func (p *Parameters) GaloisElementsForRowInnerSum() (galEls []uint64) {
	galEls = make([]uint64, p.logN+1, p.logN+1)
	galEls[0] = p.GaloisElementForRowRotation()
	for i := 0; i < int(p.logN)-1; i++ {
		galEls[i+1] = p.GaloisElementForColumnRotationBy(1 << i)
	}
	return galEls
}

// InverseGaloisElement takes a galois element and returns the galois element
//  corresponding to the inverse automorphism
func (p *Parameters) InverseGaloisElement(galEl uint64) uint64 {
	twoN := 1 << (p.logN + 1)
	return ring.ModExp(galEl, twoN-1, uint64(twoN))
}
