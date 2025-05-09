package bootstrapping

import (
	"bufio"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/buffer"
)

// EvaluationKeys is a struct storing the different
// evaluation keys required by the bootstrapper.
type EvaluationKeys struct {

	// EvkN1ToN2 is the evaluation key to switch from the residual parameters'
	// ring degree (N1) to the bootstrapping parameters' ring degree (N2)
	EvkN1ToN2 *rlwe.EvaluationKey

	// EvkN2ToN1 is the evaluation key to switch from the bootstrapping parameters'
	// ring degree (N2) to the residual parameters' ring degree (N1)
	EvkN2ToN1 *rlwe.EvaluationKey

	// EvkRealToCmplx is the evaluation key to switch from the standard ring to the
	// conjugate invariant ring.
	EvkRealToCmplx *rlwe.EvaluationKey

	// EvkCmplxToReal is the evaluation key to switch from the conjugate invariant
	// ring to the standard ring.
	EvkCmplxToReal *rlwe.EvaluationKey

	// EvkDenseToSparse is the evaluation key to switch
	// from the dense secret to the sparse secret.
	// https://eprint.iacr.org/2022/024
	EvkDenseToSparse *rlwe.EvaluationKey

	// EvkSparseToDense is the evaluation key to switch
	// from the sparse secret to the dense secret.
	// https://eprint.iacr.org/2022/024
	EvkSparseToDense *rlwe.EvaluationKey

	// MemEvaluationKeySet is the evaluation key set storing the relinearization
	// key and the Galois keys necessary for the bootstrapping circuit.
	*rlwe.MemEvaluationKeySet
}

// GenEvaluationKeys generates the bootstrapping evaluation keys, which include:
//
// If the bootstrapping parameters' ring degree > residual parameters' ring degree:
//   - An evaluation key to switch from the residual parameters' ring to the bootstrapping parameters' ring
//   - An evaluation key to switch from the bootstrapping parameters' ring to the residual parameters' ring
//
// If the residual parameters use the Conjugate Invariant ring:
//   - An evaluation key to switch from the conjugate invariant ring to the standard ring
//   - An evaluation key to switch from the standard ring to the conjugate invariant ring
//
// The core bootstrapping circuit evaluation keys:
//   - Relinearization key
//   - Galois keys
//   - The encapsulation evaluation keys (https://eprint.iacr.org/2022/024)
//
// Note:
//   - These evaluation keys are generated under an ephemeral secret key skN2 using the distribution
//     specified in the bootstrapping parameters.
//   - The ephemeral key used to generate the bootstrapping keys is returned by this method for debugging purposes.
//   - !WARNING! The bootstrapping parameters use their own and independent cryptographic parameters (i.e. float.Parameters)
//     and it is the user's responsibility to ensure that these parameters meet the target security and tweak them if necessary.
func (p Parameters) GenEvaluationKeys(skN1 *rlwe.SecretKey) (btpkeys *EvaluationKeys, skN2 *rlwe.SecretKey, err error) {

	var EvkN1ToN2, EvkN2ToN1 *rlwe.EvaluationKey
	var EvkRealToCmplx *rlwe.EvaluationKey
	var EvkCmplxToReal *rlwe.EvaluationKey
	paramsN2 := p.BootstrappingParameters

	kgen := rlwe.NewKeyGenerator(paramsN2)

	if p.ResidualParameters.N() != paramsN2.N() {
		// If the ring degree do not match
		// (if the residual parameters are Conjugate Invariant, N1 = N2/2)
		skN2 = kgen.GenSecretKeyNew()

		if p.ResidualParameters.RingType() == ring.ConjugateInvariant {
			EvkCmplxToReal, EvkRealToCmplx = kgen.GenEvaluationKeysForRingSwapNew(skN2, skN1)
		} else {
			EvkN1ToN2 = kgen.GenEvaluationKeyNew(skN1, skN2)
			EvkN2ToN1 = kgen.GenEvaluationKeyNew(skN2, skN1)
		}

	} else {

		ringQ := paramsN2.RingQ()
		ringP := paramsN2.RingP()

		// Else, keeps the same secret, but extends to the full modulus of the bootstrapping parameters.
		skN2 = rlwe.NewSecretKey(paramsN2)
		buff := ringQ.NewPoly()

		// Extends basis Q0 -> QL
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringQ, skN1.Value.Q, buff, skN2.Value.Q)

		// Extends basis Q0 -> P
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringP, skN1.Value.Q, buff, skN2.Value.P)
	}

	EvkDenseToSparse, EvkSparseToDense := p.genEncapsulationEvaluationKeysNew(skN2)

	rlk := kgen.GenRelinearizationKeyNew(skN2)
	gks := kgen.GenGaloisKeysNew(append(p.GaloisElements(paramsN2), paramsN2.GaloisElementForComplexConjugation()), skN2)

	return &EvaluationKeys{
		EvkN1ToN2:           EvkN1ToN2,
		EvkN2ToN1:           EvkN2ToN1,
		EvkRealToCmplx:      EvkRealToCmplx,
		EvkCmplxToReal:      EvkCmplxToReal,
		MemEvaluationKeySet: rlwe.NewMemEvaluationKeySet(rlk, gks...),
		EvkDenseToSparse:    EvkDenseToSparse,
		EvkSparseToDense:    EvkSparseToDense,
	}, skN2, nil
}

// GenEncapsulationEvaluationKeysNew generates the low level encapsulation EvaluationKeys for the bootstrapping.
func (p Parameters) genEncapsulationEvaluationKeysNew(skDense *rlwe.SecretKey) (EvkDenseToSparse, EvkSparseToDense *rlwe.EvaluationKey) {

	params := p.BootstrappingParameters

	if p.EphemeralSecretWeight == 0 {
		return
	}

	paramsSparse, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN: params.LogN(),
		Q:    params.Q()[:1],
		P:    params.P()[:1],
	})

	kgenSparse := rlwe.NewKeyGenerator(paramsSparse)
	kgenDense := rlwe.NewKeyGenerator(params)
	skSparse := kgenSparse.GenSecretKeyWithHammingWeightNew(p.EphemeralSecretWeight)

	EvkDenseToSparse = kgenSparse.GenEvaluationKeyNew(skDense, skSparse)
	EvkSparseToDense = kgenDense.GenEvaluationKeyNew(skSparse, skDense)
	return
}

// BinarySize returns the total binary size of the bootstrapper's keys.
func (b EvaluationKeys) BinarySize() (dLen int) {
	if b.EvkN1ToN2 != nil {
		dLen += b.EvkN1ToN2.BinarySize()
	}
	dLen++

	if b.EvkN2ToN1 != nil {
		dLen += b.EvkN2ToN1.BinarySize()
	}
	dLen++

	if b.EvkRealToCmplx != nil {
		dLen += b.EvkRealToCmplx.BinarySize()
	}
	dLen++

	if b.EvkCmplxToReal != nil {
		dLen += b.EvkCmplxToReal.BinarySize()
	}
	dLen++

	if b.EvkDenseToSparse != nil {
		dLen += b.EvkDenseToSparse.BinarySize()
	}
	dLen++

	if b.EvkSparseToDense != nil {
		dLen += b.EvkSparseToDense.BinarySize()
	}
	dLen++

	if b.MemEvaluationKeySet != nil {
		dLen += b.MemEvaluationKeySet.BinarySize()
	}
	dLen++

	return
}

// MarshalBinary encodes the object into a binary form on a newly allocated slices of bytes.
func (b EvaluationKeys) MarshalBinary() (p []byte, err error) {
	buf := buffer.NewBufferSize(b.BinarySize())
	_, err = b.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary encodes a slices of bytes generated by MarshalBinary or WriteTo on the object.
func (b *EvaluationKeys) UnmarshalBinary(p []byte) (err error) {
	_, err = b.ReadFrom(buffer.NewBuffer(p))
	return
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
func (b EvaluationKeys) WriteTo(w io.Writer) (n int64, err error) {
	switch w := w.(type) {
	case buffer.Writer:

		var inc int64

		inc, err = writeEvkKey(b.EvkN1ToN2, w)
		if err != nil {
			return inc, fmt.Errorf("cannot write EvkN1ToN2 evaluation key: %w", err)
		}
		n += inc

		inc, err = writeEvkKey(b.EvkN2ToN1, w)
		if err != nil {
			return inc, fmt.Errorf("cannot write EvkN2ToN1 evaluation key: %w", err)
		}
		n += inc

		inc, err = writeEvkKey(b.EvkRealToCmplx, w)
		if err != nil {
			return inc, fmt.Errorf("cannot write EvkRealToCmplx evaluation key: %w", err)
		}
		n += inc

		inc, err = writeEvkKey(b.EvkCmplxToReal, w)
		if err != nil {
			return inc, fmt.Errorf("cannot write EvkCmplxToReal evaluation key: %w", err)
		}
		n += inc

		inc, err = writeEvkKey(b.EvkDenseToSparse, w)
		if err != nil {
			return inc, fmt.Errorf("cannot write EvkDenseToSparse evaluation key: %w", err)
		}
		n += inc

		inc, err = writeEvkKey(b.EvkSparseToDense, w)
		if err != nil {
			return inc, fmt.Errorf("cannot write EvkSparseToDense evaluation key: %w", err)
		}
		n += inc

		if b.MemEvaluationKeySet != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return inc, err
			}
			n += inc

			if inc, err = b.MemEvaluationKeySet.WriteTo(w); err != nil {
				return n + inc, err
			}
			n += inc

		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return inc, err
			}
			n += inc
		}

		return n, w.Flush()
	default:
		return b.WriteTo(bufio.NewWriter(w))
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
func (b *EvaluationKeys) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:

		var inc int64

		b.EvkN1ToN2, inc, err = readEvkKey(r)
		if err != nil {
			return inc, fmt.Errorf("unable to read EvkN1ToN2 evaluation key: %w", err)
		}
		n += inc

		b.EvkN2ToN1, inc, err = readEvkKey(r)
		if err != nil {
			return inc, fmt.Errorf("unable to read EvkN2ToN1 evaluation key: %w", err)
		}
		n += inc

		b.EvkRealToCmplx, inc, err = readEvkKey(r)
		if err != nil {
			return inc, fmt.Errorf("unable to read EvkRealToCmplx evaluation key: %w", err)
		}
		n += inc

		b.EvkCmplxToReal, inc, err = readEvkKey(r)
		if err != nil {
			return inc, fmt.Errorf("unable to read EvkCmplxToReal evaluation key: %w", err)
		}
		n += inc

		b.EvkDenseToSparse, inc, err = readEvkKey(r)
		if err != nil {
			return inc, fmt.Errorf("unable to read EvkDenseToSparse evaluation key: %w", err)
		}
		n += inc

		b.EvkSparseToDense, inc, err = readEvkKey(r)
		if err != nil {
			return inc, fmt.Errorf("unable to read EvkSparseToDense evaluation key: %w", err)
		}
		n += inc

		var hasKey uint8

		if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
			return inc, err
		}
		n += inc

		if hasKey == 1 {

			b.MemEvaluationKeySet = new(rlwe.MemEvaluationKeySet)

			if inc, err = b.MemEvaluationKeySet.ReadFrom(r); err != nil {
				return inc, err
			}

			n += inc
		}

		return n, nil

	default:
		return b.ReadFrom(bufio.NewReader(r))
	}
}

// writeEvkKey is a helper function for marshalling a bootstrapping evaluation key.
func writeEvkKey(key *rlwe.EvaluationKey, w buffer.Writer) (n int64, err error) {
	var inc int64

	if key != nil {
		if inc, err = buffer.WriteUint8(w, 1); err != nil {
			return inc, err
		}
		n += inc

		if inc, err = key.WriteTo(w); err != nil {
			return inc, err
		}
		n += inc
	} else {
		if inc, err = buffer.WriteUint8(w, 0); err != nil {
			return inc, err
		}
		n += inc
	}
	return
}

// readEvkKey is a helper function for unmarshalling a bootstrapping evaluation key.
func readEvkKey(r buffer.Reader) (key *rlwe.EvaluationKey, n int64, err error) {
	var (
		inc    int64
		hasKey uint8
	)

	if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
		return nil, inc, err
	}
	n += inc

	if hasKey == 1 {

		key = new(rlwe.EvaluationKey)

		if inc, err = key.ReadFrom(r); err != nil {
			return nil, n + inc, err
		}

		n += inc
	}
	return
}
