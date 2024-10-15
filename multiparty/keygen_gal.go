package multiparty

import (
	"bufio"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils/buffer"
)

// GaloisKeyGenProtocol is the structure storing the parameters for the collective GaloisKeys generation.
type GaloisKeyGenProtocol struct {
	skOut ringqp.Poly
	EvaluationKeyGenProtocol
}

// GaloisKeyGenShare is represent a Party's share in the GaloisKey Generation protocol.
type GaloisKeyGenShare struct {
	GaloisElement uint64
	EvaluationKeyGenShare
}

// GaloisKeyGenCRP is a type for common reference polynomials in the GaloisKey Generation protocol.
type GaloisKeyGenCRP struct {
	EvaluationKeyGenCRP
}

// ShallowCopy creates a shallow copy of [GaloisKeyGenProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [GaloisKeyGenProtocol] can be used concurrently.
func (gkg GaloisKeyGenProtocol) ShallowCopy() GaloisKeyGenProtocol {
	return GaloisKeyGenProtocol{EvaluationKeyGenProtocol: gkg.EvaluationKeyGenProtocol.ShallowCopy(), skOut: gkg.params.RingQP().NewPoly()}
}

// NewGaloisKeyGenProtocol creates a [GaloisKeyGenProtocol] instance.
func NewGaloisKeyGenProtocol(params rlwe.ParameterProvider) (gkg GaloisKeyGenProtocol) {
	return GaloisKeyGenProtocol{EvaluationKeyGenProtocol: NewEvaluationKeyGenProtocol(params), skOut: params.GetRLWEParameters().RingQP().NewPoly()}
}

// AllocateShare allocates a party's share in the GaloisKey Generation.
func (gkg GaloisKeyGenProtocol) AllocateShare(evkParams ...rlwe.EvaluationKeyParameters) (gkgShare GaloisKeyGenShare) {
	levelQ, levelP, BaseTwoDecomposition, _ := rlwe.ResolveEvaluationKeyParameters(gkg.params, evkParams)
	return GaloisKeyGenShare{EvaluationKeyGenShare: gkg.EvaluationKeyGenProtocol.allocateShare(levelQ, levelP, BaseTwoDecomposition)}
}

// SampleCRP samples a common random polynomial to be used in the GaloisKey Generation from the provided
// common reference string.
func (gkg GaloisKeyGenProtocol) SampleCRP(crs CRS, evkParams ...rlwe.EvaluationKeyParameters) GaloisKeyGenCRP {
	levelQ, levelP, BaseTwoDecomposition, _ := rlwe.ResolveEvaluationKeyParameters(gkg.params, evkParams)
	return GaloisKeyGenCRP{gkg.EvaluationKeyGenProtocol.sampleCRP(crs, levelQ, levelP, BaseTwoDecomposition)}
}

// GenShare generates a party's share in the GaloisKey Generation.
func (gkg GaloisKeyGenProtocol) GenShare(sk *rlwe.SecretKey, galEl uint64, crp GaloisKeyGenCRP, shareOut *GaloisKeyGenShare) (err error) {

	levelQ := shareOut.LevelQ()
	levelP := shareOut.LevelP()

	ringQ := gkg.params.RingQ().AtLevel(levelQ)
	ringP := gkg.params.RingP().AtLevel(levelP)

	galElInv := ring.ModExp(galEl, ringQ.NthRoot()-1, ringQ.NthRoot())

	// Important
	shareOut.GaloisElement = galEl

	ringQ.AutomorphismNTT(sk.Value.Q, galElInv, gkg.skOut.Q)

	if levelP > -1 {
		ringP.AutomorphismNTT(sk.Value.P, galElInv, gkg.skOut.P)
	}

	return gkg.EvaluationKeyGenProtocol.GenShare(sk, &rlwe.SecretKey{Value: gkg.skOut}, crp.EvaluationKeyGenCRP, &shareOut.EvaluationKeyGenShare)

}

// AggregateShares computes share3 = share1 + share2.
func (gkg GaloisKeyGenProtocol) AggregateShares(share1, share2 GaloisKeyGenShare, share3 *GaloisKeyGenShare) (err error) {

	if share1.GaloisElement != share2.GaloisElement {
		return fmt.Errorf("cannot aggregate: GaloisKeyGenShares do not share the same GaloisElement: %d != %d", share1.GaloisElement, share2.GaloisElement)
	}

	share3.GaloisElement = share1.GaloisElement

	return gkg.EvaluationKeyGenProtocol.AggregateShares(share1.EvaluationKeyGenShare, share2.EvaluationKeyGenShare, &share3.EvaluationKeyGenShare)
}

// GenGaloisKey finalizes the GaloisKey Generation and populates the input GaloisKey with the computed collective GaloisKey.
func (gkg GaloisKeyGenProtocol) GenGaloisKey(share GaloisKeyGenShare, crp GaloisKeyGenCRP, gk *rlwe.GaloisKey) (err error) {
	if err = gkg.EvaluationKeyGenProtocol.GenEvaluationKey(share.EvaluationKeyGenShare, crp.EvaluationKeyGenCRP, &gk.EvaluationKey); err != nil {
		return
	}
	gk.GaloisElement = share.GaloisElement
	gk.NthRoot = gkg.params.RingQ().NthRoot()
	return
}

// BinarySize returns the serialized size of the object in bytes.
func (share GaloisKeyGenShare) BinarySize() int {
	return 8 + share.EvaluationKeyGenShare.BinarySize()
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the [buffer.Writer] interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a [bufio.Writer]. Since this requires allocations, it
// is preferable to pass a [buffer.Writer] directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated [bufio.Writer].
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (share GaloisKeyGenShare) WriteTo(w io.Writer) (n int64, err error) {
	switch w := w.(type) {
	case buffer.Writer:
		var inc int64

		if inc, err = buffer.WriteUint64(w, share.GaloisElement); err != nil {
			return n + inc, err
		}

		n += inc

		if inc, err = share.EvaluationKeyGenShare.WriteTo(w); err != nil {
			return n + inc, err
		}

		n += inc

		return n, err

	default:
		return share.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap [io.Reader] in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (share *GaloisKeyGenShare) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:

		var inc int64
		if inc, err = buffer.ReadUint64(r, &share.GaloisElement); err != nil {
			return n + inc, err
		}
		n += inc

		if inc, err = share.EvaluationKeyGenShare.ReadFrom(r); err != nil {
			return n + inc, err
		}

		return n + inc, nil
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
// [GaloisKeyGenShare.MarshalBinary] or [GaloisKeyGenShare.WriteTo] on the object.
func (share *GaloisKeyGenShare) UnmarshalBinary(p []byte) (err error) {
	_, err = share.ReadFrom(buffer.NewBuffer(p))
	return
}
