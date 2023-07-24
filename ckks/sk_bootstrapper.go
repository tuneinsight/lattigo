package ckks

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// SecretKeyBootstrapper is an implementation of the rlwe.Bootstrapping interface that
// uses the secret-key to decrypt and re-encrypt the bootstrapped ciphertext.
type SecretKeyBootstrapper struct {
	Parameters
	*Encoder
	*rlwe.Decryptor
	*rlwe.Encryptor
	sk      *rlwe.SecretKey
	Values  []*bignum.Complex
	Counter int // records the number of bootstrapping
}

func NewSecretKeyBootstrapper(params Parameters, sk *rlwe.SecretKey) (rlwe.Bootstrapper, error) {

	dec, err := NewDecryptor(params, sk)

	if err != nil {
		return nil, err
	}

	enc, err := NewEncryptor(params, sk)

	if err != nil {
		return nil, err
	}

	return &SecretKeyBootstrapper{
		params,
		NewEncoder(params),
		dec,
		enc,
		sk,
		make([]*bignum.Complex, params.N()),
		0}, nil
}

func (d *SecretKeyBootstrapper) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	values := d.Values[:1<<ct.LogDimensions.Cols]
	if err := d.Decode(d.DecryptNew(ct), values); err != nil {
		return nil, err
	}
	pt := NewPlaintext(d.Parameters, d.MaxLevel())
	pt.MetaData = ct.MetaData
	pt.Scale = d.Parameters.DefaultScale()
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

func (d SecretKeyBootstrapper) BootstrapMany(cts []*rlwe.Ciphertext) ([]*rlwe.Ciphertext, error) {
	for i := range cts {
		cts[i], _ = d.Bootstrap(cts[i])
	}
	return cts, nil
}

func (d SecretKeyBootstrapper) Depth() int {
	return 0
}

func (d SecretKeyBootstrapper) MinimumInputLevel() int {
	return 0
}

func (d SecretKeyBootstrapper) OutputLevel() int {
	return d.MaxLevel()
}
