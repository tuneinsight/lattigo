// Package heint implements Homomorphic Encryption for encrypted modular arithmetic over the integers.
package heint

import (
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
)

type Integer interface {
	bgv.Integer
}

type ParametersLiteral bgv.ParametersLiteral

func NewParametersFromLiteral(paramsLit ParametersLiteral) (Parameters, error) {
	params, err := bgv.NewParametersFromLiteral(bgv.ParametersLiteral(paramsLit))
	return Parameters{Parameters: params}, err
}

type Parameters struct {
	bgv.Parameters
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
	return bgv.NewPlaintext(params.Parameters, level)
}

func NewCiphertext(params Parameters, degree, level int) *rlwe.Ciphertext {
	return bgv.NewCiphertext(params.Parameters, degree, level)
}

func NewEncryptor(params Parameters, key rlwe.EncryptionKey) *rlwe.Encryptor {
	return rlwe.NewEncryptor(params, key)
}

func NewDecryptor(params Parameters, key *rlwe.SecretKey) *rlwe.Decryptor {
	return rlwe.NewDecryptor(params, key)
}

func NewKeyGenerator(params Parameters) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params)
}

type Encoder struct {
	bgv.Encoder
}

func NewEncoder(params Parameters) *Encoder {
	return &Encoder{Encoder: *bgv.NewEncoder(params.Parameters)}
}

func (ecd Encoder) ShallowCopy() *Encoder {
	return &Encoder{Encoder: *ecd.Encoder.ShallowCopy()}
}

type Evaluator struct {
	bgv.Evaluator
}

func NewEvaluator(params Parameters, evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{Evaluator: *bgv.NewEvaluator(params.Parameters, evk)}
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
