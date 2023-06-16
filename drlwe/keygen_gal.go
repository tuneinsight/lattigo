package drlwe

import (
	"bufio"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// GaloisKeyGenProtocol is the structure storing the parameters for the collective GaloisKeys generation.
type GaloisKeyGenProtocol struct {
	params           rlwe.Parameters
	buff             [2]*ringqp.Poly
	gaussianSamplerQ ring.Sampler
}

// GaloisKeyGenShare is represent a Party's share in the GaloisKey Generation protocol.
type GaloisKeyGenShare struct {
	GaloisElement uint64
	Value         structs.Matrix[ringqp.Poly]
}

// GaloisKeyGenCRP is a type for common reference polynomials in the GaloisKey Generation protocol.
type GaloisKeyGenCRP struct {
	Value structs.Matrix[ringqp.Poly]
}

// ShallowCopy creates a shallow copy of GaloisKeyGenProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// GaloisKeyGenProtocol can be used concurrently.
func (gkg *GaloisKeyGenProtocol) ShallowCopy() GaloisKeyGenProtocol {
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := gkg.params

	return GaloisKeyGenProtocol{
		params:           gkg.params,
		buff:             [2]*ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
		gaussianSamplerQ: ring.NewSampler(prng, gkg.params.RingQ(), gkg.params.Xe(), false),
	}
}

// NewGaloisKeyGenProtocol creates a GaloisKeyGenProtocol instance.
func NewGaloisKeyGenProtocol(params rlwe.Parameters) (gkg GaloisKeyGenProtocol) {
	gkg = GaloisKeyGenProtocol{}
	gkg.params = params

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}
	gkg.buff = [2]*ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()}
	gkg.gaussianSamplerQ = ring.NewSampler(prng, params.RingQ(), params.Xe(), false)
	return
}

// AllocateShare allocates a party's share in the GaloisKey Generation.
func (gkg GaloisKeyGenProtocol) AllocateShare() (gkgShare GaloisKeyGenShare) {
	params := gkg.params
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	p := make([][]ringqp.Poly, decompRNS)
	for i := range p {
		vec := make([]ringqp.Poly, decompPw2)
		for j := range vec {
			vec[j] = *ringqp.NewPoly(params.N(), params.MaxLevelQ(), params.MaxLevelP())
		}
		p[i] = vec
	}

	return GaloisKeyGenShare{Value: structs.Matrix[ringqp.Poly](p)}
}

// SampleCRP samples a common random polynomial to be used in the GaloisKey Generation from the provided
// common reference string.
func (gkg GaloisKeyGenProtocol) SampleCRP(crs CRS) GaloisKeyGenCRP {

	params := gkg.params
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	m := make([][]ringqp.Poly, decompRNS)
	for i := range m {
		vec := make([]ringqp.Poly, decompPw2)
		for j := range vec {
			vec[j] = *ringqp.NewPoly(params.N(), params.MaxLevelQ(), params.MaxLevelP())
		}
		m[i] = vec
	}

	us := ringqp.NewUniformSampler(crs, *params.RingQP())

	for _, v := range m {
		for _, p := range v {
			us.Read(&p)
		}
	}

	return GaloisKeyGenCRP{Value: structs.Matrix[ringqp.Poly](m)}
}

// GenShare generates a party's share in the GaloisKey Generation.
func (gkg GaloisKeyGenProtocol) GenShare(sk *rlwe.SecretKey, galEl uint64, crp GaloisKeyGenCRP, shareOut *GaloisKeyGenShare) {

	ringQ := gkg.params.RingQ()
	ringQP := gkg.params.RingQP()

	levelQ := sk.LevelQ()
	levelP := sk.LevelP()

	galElInv := ring.ModExp(galEl, ringQ.NthRoot()-1, ringQ.NthRoot())

	// Important
	shareOut.GaloisElement = galEl

	ringQ.AutomorphismNTT(sk.Value.Q, galElInv, gkg.buff[1].Q)

	var hasModulusP bool

	if levelP > -1 {
		hasModulusP = true
		gkg.params.RingP().AutomorphismNTT(sk.Value.P, galElInv, gkg.buff[1].P)
		ringQ.MulScalarBigint(sk.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], gkg.buff[0].Q)
	} else {
		levelP = 0
		ring.CopyLvl(levelQ, sk.Value.Q, gkg.buff[0].Q)
	}

	m := shareOut.Value
	c := crp.Value

	RNSDecomp := len(m)
	BITDecomp := len(m[0])

	N := ringQ.N()

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e
			gkg.gaussianSamplerQ.Read(m[i][j].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(m[i][j].Q, levelP, nil, m[i][j].P)
			}

			ringQP.NTTLazy(&m[i][j], &m[i][j])
			ringQP.MForm(&m[i][j], &m[i][j])

			// a is the CRP

			// e + sk_in * (qiBarre*qiStar) * 2^w
			// (qiBarre*qiStar)%qi = 1, else 0
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handles the case where nb pj does not divides nb qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.SubRings[index].Modulus
				tmp0 := gkg.buff[0].Q.Coeffs[index]
				tmp1 := m[i][j].Q.Coeffs[index]

				for w := 0; w < N; w++ {
					tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
				}
			}

			// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
			ringQP.MulCoeffsMontgomeryThenSub(&c[i][j], gkg.buff[1], &m[i][j])
		}

		ringQ.MulScalar(gkg.buff[0].Q, 1<<gkg.params.Pow2Base(), gkg.buff[0].Q)
	}
}

// AggregateShares computes share3 = share1 + share2.
func (gkg GaloisKeyGenProtocol) AggregateShares(share1, share2, share3 *GaloisKeyGenShare) {

	if share1.GaloisElement != share2.GaloisElement {
		panic(fmt.Sprintf("cannot aggregate: GaloisKeyGenShares do not share the same GaloisElement: %d != %d", share1.GaloisElement, share2.GaloisElement))
	}

	share3.GaloisElement = share1.GaloisElement

	m1 := share1.Value
	m2 := share2.Value
	m3 := share3.Value

	levelQ := m1[0][0].Q.Level()

	var levelP int
	if m1[0][0].P != nil {
		levelP = m1[0][0].P.Level()
	}

	ringQP := gkg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(m1)
	BITDecomp := len(m1[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(&m1[i][j], &m2[i][j], &m3[i][j])
		}
	}
}

// GenGaloisKey finalizes the GaloisKey Generation and populates the input GaloisKey with the computed collective GaloisKey.
func (gkg GaloisKeyGenProtocol) GenGaloisKey(share GaloisKeyGenShare, crp GaloisKeyGenCRP, gk *rlwe.GaloisKey) {

	m := share.Value
	p := crp.Value

	RNSDecomp := len(m)
	BITDecomp := len(m[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			gk.Value[i][j][0].Copy(&m[i][j])
			gk.Value[i][j][1].Copy(&p[i][j])
		}
	}

	gk.GaloisElement = share.GaloisElement
	gk.NthRoot = gkg.params.RingQ().NthRoot()
}

// BinarySize returns the serialized size of the object in bytes.
func (share GaloisKeyGenShare) BinarySize() int {
	return 8 + share.Value.BinarySize()
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (share GaloisKeyGenShare) WriteTo(w io.Writer) (n int64, err error) {
	switch w := w.(type) {
	case buffer.Writer:
		var inc int

		if inc, err = buffer.WriteUint64(w, share.GaloisElement); err != nil {
			return n + int64(inc), err
		}

		n += int64(inc)

		var inc2 int64
		if inc2, err = share.Value.WriteTo(w); err != nil {
			return n + inc2, err
		}

		n += inc2

		return n, err

	default:
		return share.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// Unless r implements the buffer.Reader interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (share *GaloisKeyGenShare) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:

		var inc int
		if inc, err = buffer.ReadUint64(r, &share.GaloisElement); err != nil {
			return n + int64(inc), err
		}
		n += int64(inc)

		var inc64 int64
		if inc64, err = share.Value.ReadFrom(r); err != nil {
			return n + inc64, err
		}

		return n + inc64, nil
	default:
		return share.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share GaloisKeyGenShare) MarshalBinary() (p []byte, err error) {
	buf := buffer.NewBufferSize(share.BinarySize())
	_, err = share.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (share *GaloisKeyGenShare) UnmarshalBinary(p []byte) (err error) {
	_, err = share.ReadFrom(buffer.NewBuffer(p))
	return
}
