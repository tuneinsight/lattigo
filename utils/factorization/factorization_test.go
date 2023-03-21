package factorization_test

import (
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/utils/factorization"
)

func TestGetFactors(t *testing.T) {

	m := new(big.Int).SetUint64(35184372088631)

	t.Run("ECM", func(t *testing.T) {

		factor := factorization.GetFactorECM(m)

		if factor.Cmp(new(big.Int).SetUint64(6292343)) != 0 && factor.Cmp(new(big.Int).SetUint64(5591617)) != 0 {
			t.Fail()
		}
	})

	t.Run("PollardRho", func(t *testing.T) {
		factor := factorization.GetFactorPollardRho(m)

		if factor.Cmp(new(big.Int).SetUint64(6292343)) != 0 && factor.Cmp(new(big.Int).SetUint64(5591617)) != 0 {
			t.Fail()
		}
	})
}
