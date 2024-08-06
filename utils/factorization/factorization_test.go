package factorization_test

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/utils/factorization"
)

const (
	prime uint64 = 0x1fffffffffe00001
)

func TestIsPrime(t *testing.T) {
	// 2^64 - 59 is prime
	require.True(t, factorization.IsPrime(new(big.Int).SetUint64(0xffffffffffffffc5)))
	// 2^64 + 13 is prime
	bigPrime, _ := new(big.Int).SetString("18446744073709551629", 10)
	require.True(t, factorization.IsPrime(bigPrime))
	// 2^64 is not prime
	require.False(t, factorization.IsPrime(new(big.Int).SetUint64(0xffffffffffffffff)))
}

func TestGetFactors(t *testing.T) {

	t.Run("GetFactors", func(t *testing.T) {
		m := new(big.Int).SetUint64(prime - 1)
		require.True(t, checkFactorization(new(big.Int).Set(m), factorization.GetFactors(m)))
	})

	t.Run("ECM", func(t *testing.T) {
		m := new(big.Int).SetUint64(prime - 1)
		require.True(t, m.Mod(m, factorization.GetFactorECM(m)).Cmp(new(big.Int)) == 0)
	})

	t.Run("PollardRho", func(t *testing.T) {
		m := new(big.Int).SetUint64(prime - 1)
		require.True(t, m.Mod(m, factorization.GetFactorPollardRho(m)).Cmp(new(big.Int)) == 0)
	})
}

func checkFactorization(p *big.Int, factors []*big.Int) bool {
	zero := new(big.Int)
	for _, factor := range factors {
		for new(big.Int).Mod(p, factor).Cmp(zero) == 0 {
			p.Quo(p, factor)
		}
	}

	return p.Cmp(new(big.Int).SetUint64(1)) == 0
}
