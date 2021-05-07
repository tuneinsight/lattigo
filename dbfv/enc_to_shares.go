package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

type E2SProtocol struct {
	CKSProtocol

	ringQ       *ring.Ring
	ringT       *ring.Ring
	maskSampler *ring.UniformSampler
	encoder     bfv.Encoder

	zero              *rlwe.SecretKey
	tmpPlaintextRingT *bfv.PlaintextRingT
	tmpPlaintext      *bfv.Plaintext
}

type AdditiveShare struct {
	Value ring.Poly
}

func NewAdditiveShare(params bfv.Parameters) AdditiveShare {
	return AdditiveShare{Value: *ring.NewPoly(params.N(), 1)}
}

func NewE2SProtocol(params bfv.Parameters, sigmaSmudging float64) *E2SProtocol {
	e2s := new(E2SProtocol)
	e2s.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	e2s.ringQ = params.RingQ()
	e2s.ringT = params.RingT()
	e2s.encoder = bfv.NewEncoder(params)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	e2s.maskSampler = ring.NewUniformSampler(prng, e2s.ringT)
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.tmpPlaintext = bfv.NewPlaintext(params)
	e2s.tmpPlaintextRingT = bfv.NewPlaintextRingT(params)
	return e2s
}

func (e2s *E2SProtocol) GenShare(sk *rlwe.SecretKey, ct *bfv.Ciphertext, secretShareOut AdditiveShare, publicShareOut *drlwe.CKSShare) {
	e2s.CKSProtocol.GenShare(sk, e2s.zero, ct, publicShareOut)
	e2s.maskSampler.Read(&secretShareOut.Value)
	//e2s.encoder.EncodeUint(secretShareOut.Value.Coeffs[0], e2s.tmpPlaintext)
	e2s.tmpPlaintextRingT.Value[0].Copy(&secretShareOut.Value)
	e2s.encoder.ScaleUp(e2s.tmpPlaintextRingT, e2s.tmpPlaintext)
	e2s.ringQ.Sub(publicShareOut.Value, e2s.tmpPlaintext.Value[0], publicShareOut.Value)
}

func (e2s *E2SProtocol) Finalize(secretShare *AdditiveShare, aggregatePublicShare *drlwe.CKSShare, ct *bfv.Ciphertext, secretShareOut *AdditiveShare) {
	e2s.ringQ.Add(aggregatePublicShare.Value, ct.Value[0], e2s.tmpPlaintext.Value[0])
	//e2s.encoder.DecodeUint(e2s.tmpPlaintext, e2s.tmpPlaintextRingT.Value[0].Coeffs[0])
	e2s.encoder.ScaleDown(e2s.tmpPlaintext, e2s.tmpPlaintextRingT)
	if secretShare != nil {
		e2s.ringT.Add(&secretShare.Value, e2s.tmpPlaintextRingT.Value[0], &secretShareOut.Value)
	} else {
		secretShareOut.Value.Copy(e2s.tmpPlaintextRingT.Value[0])
	}
}
