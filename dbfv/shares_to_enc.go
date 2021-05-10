package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

type S2EProtocol struct {
	CKSProtocol

	ringQ   *ring.Ring
	encoder bfv.Encoder

	zero         *rlwe.SecretKey
	tmpPlaintext *bfv.Plaintext
}

func NewS2EProtocol(params bfv.Parameters, sigmaSmudging float64) *S2EProtocol {
	e2s := new(S2EProtocol)
	e2s.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	e2s.ringQ = params.RingQ()
	e2s.encoder = bfv.NewEncoder(params)
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.tmpPlaintext = bfv.NewPlaintext(params)
	return e2s
}

func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, c1 *ring.Poly, secretShare AdditiveShare, c0ShareOut *drlwe.CKSShare) {
	s2e.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: &secretShare.Value}}, s2e.tmpPlaintext)
	s2e.CKSProtocol.GenShare(s2e.zero, sk, &rlwe.Element{Value: []*ring.Poly{c0ShareOut.Value, c1}}, c0ShareOut)
	s2e.ringQ.Add(c0ShareOut.Value, s2e.tmpPlaintext.Value, c0ShareOut.Value)
}

func (s2e *S2EProtocol) Finalize(c0Agg *drlwe.CKSShare, c1 *ring.Poly, ctOut bfv.Ciphertext) {
	if ctOut.Degree() != 1 {
		panic("ctOut must have degree 1.")
	}
	ctOut.Value[0].Copy(c0Agg.Value)
	ctOut.Value[1].Copy(c1)
}
