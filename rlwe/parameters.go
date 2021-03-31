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

type Parameters struct {
	logN  uint64
	qi    []uint64
	pi    []uint64
	sigma float64
}

func NewRLWEParameters(logn uint64, q, p []uint64, sigma float64) *Parameters { // TEMPORARY constructor
	return &Parameters{
		logN:  logn,
		pi:    p,
		qi:    q,
		sigma: sigma,
	}
}

// N returns the ring degree
func (p *Parameters) N() uint64 {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) LogN() uint64 {
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
func (p *Parameters) QiCount() uint64 {
	return uint64(len(p.qi))
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
func (p *Parameters) PiCount() uint64 {
	return uint64(len(p.pi))
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
func (p *Parameters) QPCount() uint64 {
	return uint64(len(p.qi) + len(p.pi))
}

func (p *Parameters) QPBigInt() *big.Int {
	pqInt := p.QBigInt()
	pqInt.Mul(pqInt, p.PBigInt())
	return pqInt
}

// LogQP returns the size of the extended modulus QP in bits
func (p *Parameters) LogQP() uint64 {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return uint64(tmp.BitLen())
}

// Alpha returns the number of moduli in in P
func (p *Parameters) Alpha() uint64 {
	return p.PiCount()
}

// Beta returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (p *Parameters) Beta() uint64 {
	if p.Alpha() != 0 {
		return uint64(math.Ceil(float64(p.QiCount()) / float64(p.Alpha())))
	}

	return 0
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
