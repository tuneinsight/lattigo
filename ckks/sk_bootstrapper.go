package ckks

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// SecretKeyBootstrapper is an implementation of the rlwe.Bootstrapping interface that
// uses the secret-key to decrypt and re-encrypt the bootstrapped ciphertext.
type SecretKeyBootstrapper struct {
	Parameters
	*Encoder
	*rlwe.Decryptor
	rlwe.EncryptorInterface
	Values            []*bignum.Complex
	Counter           int // records the number of bootstrapping
	minimumInputLevel int
	outputLevel       int
}

func NewSecretKeyBootstrapper(params Parameters, sk *rlwe.SecretKey, MinimumInputLevel, OutputLevel int) (rlwe.Bootstrapper, error) {

	dec, err := NewDecryptor(params, sk)

	if err != nil {
		return nil, err
	}

	enc, err := NewEncryptor(params, sk)

	if err != nil {
		return nil, err
	}

	return &SecretKeyBootstrapper{
		Parameters:         params,
		Encoder:            NewEncoder(params),
		EncryptorInterface: enc,
		Decryptor:          dec,
		Values:             make([]*bignum.Complex, params.N()),
		Counter:            0,
		minimumInputLevel:  MinimumInputLevel,
		outputLevel:        OutputLevel,
	}, nil
}

func (d SecretKeyBootstrapper) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	if ct.Level() < d.MinimumInputLevel() {
		return nil, fmt.Errorf("cannot Bootstrap: input ciphertext.Level()=%d < MinimumInputLevel=%d", ct.Level(), d.MinimumInputLevel())
	}

	values := d.Values[:1<<ct.PlaintextLogDimensions[1]]
	if err := d.Decode(d.DecryptNew(ct), values); err != nil {
		return nil, err
	}

	pt := NewPlaintext(d.Parameters, d.MaxLevel())
	pt.MetaData = ct.MetaData
	pt.PlaintextScale = d.parameters.PlaintextScale()

	if err := d.Encode(values, pt); err != nil {
		return nil, err
	}

	ct.Resize(1, d.MaxLevel())
	if err := d.Encrypt(pt, ct); err != nil {
		return nil, err
	}

	d.Counter++
	return ct, nil
}

func (d SecretKeyBootstrapper) BootstrapMany(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error) {
	var err error
	for i := range cts {

		var ct *rlwe.Ciphertext
		if ct, err = d.Bootstrap(&cts[i]); err != nil {
			return nil, err
		}

		cts[i] = *ct
	}
	return cts, nil
}

func (d SecretKeyBootstrapper) MinimumInputLevel() int {
	return d.minimumInputLevel
}

func (d SecretKeyBootstrapper) OutputLevel() int {
	return d.outputLevel
}
