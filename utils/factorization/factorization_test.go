package factorization_test

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/tuneinsight/lattigo/v4/utils/factorization"
)

func TestIsPrime(t *testing.T) {
	// 2^64 - 59 is prime
	assert.True(t, factorization.IsPrime(new(big.Int).SetUint64(0xffffffffffffffc5)))
	// 2^64 + 13 is prime
	bigPrime, _ := new(big.Int).SetString("18446744073709551629", 10)
	assert.True(t, factorization.IsPrime(bigPrime))
	// 2^64 is not prime
	assert.False(t, factorization.IsPrime(new(big.Int).SetUint64(0xffffffffffffffff)))
}

func TestGetFactors(t *testing.T) {

	m := new(big.Int).SetUint64(35184372088631)

	t.Run("GetFactors", func(t *testing.T) {
		factors := factorization.GetFactors(m)
		if factors[0].Cmp(new(big.Int).SetUint64(5591617)) != 0 && factors[0].Cmp(new(big.Int).SetUint64(6292343)) != 0 {
			t.Fail()
		}
	})

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
