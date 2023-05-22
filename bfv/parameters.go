package bfv

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

var (
	// PN11QP54 is a set of default parameters with logN=11 and logQP=54
	PN11QP54 = ParametersLiteral{
		LogN:     11,
		Q:        []uint64{0x3001, 0x15400000001}, // 13.5 + 40.4 bits
		Pow2Base: 6,
		T:        0x3001,
	}

	// PN12QP109 is a set of default parameters with logN=12 and logQP=109
	PN12QP109 = ParametersLiteral{
		LogN: 12,
		Q:    []uint64{0x7ffffec001, 0x8000016001}, // 39 + 39 bits
		P:    []uint64{0x40002001},                 // 30 bits
		T:    65537,
	}
	// PN13QP218 is a set of default parameters with logN=13 and logQP=218
	PN13QP218 = ParametersLiteral{
		LogN: 13,
		Q:    []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001}, // 54 + 54 + 54 bits
		P:    []uint64{0x7ffffffffb4001},                                     // 55 bits
		T:    65537,
	}

	// PN14QP438 is a set of default parameters with logN=14 and logQP=438
	PN14QP438 = ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
			0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
		P: []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
		T: 65537,
	}

	// PN15QP880 is a set of default parameters with logN=15 and logQP=880
	PN15QP880 = ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
			0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
			0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
			0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
		P: []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
		T: 65537,
	}

	// PN12QP101pq is a set of default (post quantum) parameters with logN=12 and logQP=101
	PN12QP101pq = ParametersLiteral{ // LogQP = 101.00005709794536
		LogN: 12,
		Q:    []uint64{0x800004001, 0x800008001}, // 2*35
		P:    []uint64{0x80014001},               // 1*31
		T:    65537,
	}

	// PN13QP202pq is a set of default (post quantum) parameters with logN=13 and logQP=202
	PN13QP202pq = ParametersLiteral{ // LogQP = 201.99999999994753
		LogN: 13,
		Q:    []uint64{0x7fffffffe0001, 0x7fffffffcc001, 0x3ffffffffc001}, // 2*51 + 50
		P:    []uint64{0x4000000024001},                                   // 50,
		T:    65537,
	}

	// PN14QP411pq is a set of default (post quantum) parameters with logN=14 and logQP=411
	PN14QP411pq = ParametersLiteral{ // LogQP = 410.9999999999886
		LogN: 14,
		Q:    []uint64{0x7fffffffff18001, 0x8000000000f8001, 0x7ffffffffeb8001, 0x800000000158001, 0x7ffffffffe70001}, // 5*59
		P:    []uint64{0x7ffffffffe10001, 0x400000000068001},                                                          // 59+58
		T:    65537,
	}

	// PN15QP827pq is a set of default (post quantum) parameters with logN=15 and logQP=827
	PN15QP827pq = ParametersLiteral{ // LogQP = 826.9999999999509
		LogN: 15,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, 0x7ffffffffba0001, 0x8000000004a0001,
			0x7ffffffffb00001, 0x800000000890001, 0x8000000009d0001, 0x7ffffffff630001, 0x800000000a70001,
			0x7ffffffff510001}, // 11*59
		P: []uint64{0x800000000b80001, 0x800000000bb0001, 0xffffffffffc0001}, // 2*59+60
		T: 65537,
	}
)

func NewParameters(rlweParams rlwe.Parameters, t uint64) (p Parameters, err error) {
	var pbgv bgv.Parameters
	pbgv, err = bgv.NewParameters(rlweParams, t)
	return Parameters(pbgv), err
}

func NewParametersFromLiteral(pl ParametersLiteral) (p Parameters, err error) {
	var pbgv bgv.Parameters
	pbgv, err = bgv.NewParametersFromLiteral(bgv.ParametersLiteral(pl))
	return Parameters(pbgv), err
}

type ParametersLiteral bgv.ParametersLiteral

func (p ParametersLiteral) RLWEParametersLiteral() rlwe.ParametersLiteral {
	return bgv.ParametersLiteral(p).RLWEParametersLiteral()
}

type Parameters bgv.Parameters

func (p Parameters) ParametersLiteral() ParametersLiteral {
	return ParametersLiteral(bgv.Parameters(p).ParametersLiteral())
}

func (p Parameters) RingQMul() *ring.Ring {
	return bgv.Parameters(p).RingQMul()
}

func (p Parameters) T() uint64 {
	return bgv.Parameters(p).T()
}

func (p Parameters) LogT() float64 {
	return bgv.Parameters(p).LogT()
}

func (p Parameters) RingT() *ring.Ring {
	return bgv.Parameters(p).RingT()
}

func (p Parameters) Equal(other Parameters) bool {
	return bgv.Parameters(p).Equal(bgv.Parameters(other))
}

func (p Parameters) CopyNew() Parameters {
	return Parameters(bgv.Parameters(p))
}

func (p Parameters) MarshalBinary() (data []byte, err error) {
	return p.MarshalJSON()
}

func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	return p.UnmarshalJSON(data)
}

// MarshalJSON returns a JSON representation of this parameter set. See `Marshal` from the `encoding/json` package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	return bgv.Parameters(p).MarshalJSON()
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalJSON(data []byte) (err error) {

	pp := bgv.Parameters(*p)

	if err = pp.UnmarshalJSON(data); err != nil {
		return
	}

	*p = Parameters(pp)

	return
}
