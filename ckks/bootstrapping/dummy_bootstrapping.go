package bootstrapping

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// DummyBootstrapper is a dummy bootstrapper that complies to
// the rlwe.Bootstrapping interface
type DummyBootstrapper struct {
	ckks.Parameters
	ckks.Encoder
	rlwe.Decryptor
	rlwe.Encryptor
	Counter int // records the number of bootstrapping
}

func NewDummyBootstrapper(params ckks.Parameters, sk *rlwe.SecretKey) *DummyBootstrapper {
	return &DummyBootstrapper{params, ckks.NewEncoder(params), ckks.NewDecryptor(params, sk), ckks.NewEncryptor(params, sk), 0}
}

func (dbtp *DummyBootstrapper) Bootstrap(ct interface{}) (err error) {

	switch ct := ct.(type) {
	case *rlwe.Ciphertext:
		pt := dbtp.EncodeNew(dbtp.Decode(dbtp.DecryptNew(ct), dbtp.LogSlots()), dbtp.MaxLevel(), dbtp.DefaultScale(), dbtp.LogSlots())
		ct.Resize(1, dbtp.MaxLevel())
		dbtp.Encrypt(pt, ct)
		dbtp.Counter++
	default:
		return fmt.Errorf("DummyBootstrapper: invalid ciphertext type, must be *ckks.Ciphertext")
	}

	return
}

func (dbtp *DummyBootstrapper) Depth() int {
	return 0
}

func (dbtp *DummyBootstrapper) StartingLevel() int {
	return 0
}

func (dbtp *DummyBootstrapper) EndLevel() int {
	return dbtp.MaxLevel()
}
