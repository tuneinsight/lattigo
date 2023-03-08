package bootstrapping

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// SecretKeyBootstrapper is an implementation of the rlwe.Bootstrapping interface that
// uses the secret-key to decrypt and re-encrypt the bootstrapped ciphertext.
type SecretKeyBootstrapper struct {
	ckks.Parameters
	ckks.Encoder
	rlwe.Decryptor
	rlwe.Encryptor
	Counter int // records the number of bootstrapping
}

func NewSecretKeyBootstrapper(params ckks.Parameters, sk *rlwe.SecretKey) rlwe.Bootstrapper {
	return &SecretKeyBootstrapper{params, ckks.NewEncoder(params), ckks.NewDecryptor(params, sk), ckks.NewEncryptor(params, sk), 0}
}

func (d *SecretKeyBootstrapper) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	pt := d.EncodeNew(d.Decode(d.DecryptNew(ct), d.LogSlots()), d.MaxLevel(), d.DefaultScale(), d.LogSlots())
	ct.Resize(1, d.MaxLevel())
	d.Encrypt(pt, ct)
	d.Counter++
	return ct, nil
}

func (d *SecretKeyBootstrapper) BootstrapMany(cts []*rlwe.Ciphertext) ([]*rlwe.Ciphertext, error) {
	for i := range cts {
		cts[i], _ = d.Bootstrap(cts[i])
	}
	return cts, nil
}

func (d *SecretKeyBootstrapper) Depth() int {
	return 0
}

func (d *SecretKeyBootstrapper) MinimumInputLevel() int {
	return 0
}

func (d *SecretKeyBootstrapper) OuputLevel() int {
	return d.MaxLevel()
}
