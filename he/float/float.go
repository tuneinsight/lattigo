// Package float implements Homomorphic Encryption for fixed-point approximate arithmetic over the reals/complexes.
package float

import (
	"testing"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type Float interface {
	ckks.Float
}

type ParametersLiteral ckks.ParametersLiteral

func NewParametersFromLiteral(paramsLit ParametersLiteral) (Parameters, error) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral(paramsLit))
	return Parameters{Parameters: params}, err
}

type Parameters struct {
	ckks.Parameters
}

func (p Parameters) MarshalJSON() (d []byte, err error) {
	return p.Parameters.MarshalJSON()
}

func (p *Parameters) UnmarshalJSON(d []byte) (err error) {
	return p.Parameters.UnmarshalJSON(d)
}

func (p Parameters) MarshalBinary() (d []byte, err error) {
	return p.Parameters.MarshalBinary()
}

func (p *Parameters) UnmarshalBinary(d []byte) (err error) {
	return p.Parameters.UnmarshalBinary(d)
}

func (p Parameters) Equal(other *Parameters) bool {
	return p.Parameters.Equal(&other.Parameters)
}

func NewPlaintext(params Parameters, level int) *rlwe.Plaintext {
	return ckks.NewPlaintext(params.Parameters, level)
}

func NewCiphertext(params Parameters, degree, level int) *rlwe.Ciphertext {
	return ckks.NewCiphertext(params.Parameters, degree, level)
}

type Encoder struct {
	ckks.Encoder
}

func NewEncoder(params Parameters, prec ...uint) *Encoder {

	var ecd *ckks.Encoder
	if len(prec) == 0 {
		ecd = ckks.NewEncoder(params.Parameters)
	} else {
		ecd = ckks.NewEncoder(params.Parameters, prec[0])
	}

	return &Encoder{Encoder: *ecd}
}

func (ecd Encoder) ShallowCopy() *Encoder {
	return &Encoder{Encoder: *ecd.Encoder.ShallowCopy()}
}

type Evaluator struct {
	ckks.Evaluator
}

func NewEvaluator(params Parameters, evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{Evaluator: *ckks.NewEvaluator(params.Parameters, evk)}
}

func (eval Evaluator) GetParameters() *Parameters {
	return &Parameters{*eval.Evaluator.GetParameters()}
}

func (eval Evaluator) WithKey(evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{Evaluator: *eval.Evaluator.WithKey(evk)}
}

func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{Evaluator: *eval.Evaluator.ShallowCopy()}
}

func GetPrecisionStats(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, want, have interface{}, logprec float64, computeDCF bool) (prec ckks.PrecisionStats) {
	return ckks.GetPrecisionStats(params.Parameters, &encoder.Encoder, decryptor, want, have, logprec, computeDCF)
}

func VerifyTestVectors(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, valuesWant, valuesHave interface{}, log2MinPrec int, logprec float64, printPrecisionStats bool, t *testing.T) {
	ckks.VerifyTestVectors(params.Parameters, &encoder.Encoder, decryptor, valuesWant, valuesHave, log2MinPrec, logprec, printPrecisionStats, t)
}
