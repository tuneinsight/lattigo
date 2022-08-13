package utils

import (
	"testing"

	"github.com/stretchr/testify/require"
)

func Test_PRNG(t *testing.T) {

	t.Run("PRNG", func(t *testing.T) {

		key := []byte{0x49, 0x0a, 0x42, 0x3d, 0x97, 0x9d, 0xc1, 0x07, 0xa1, 0xd7, 0xe9, 0x7b, 0x3b, 0xce, 0xa1, 0xdb,
			0x42, 0xf3, 0xa6, 0xd5, 0x75, 0xd2, 0x0c, 0x92, 0xb7, 0x35, 0xce, 0x0c, 0xee, 0x09, 0x7c, 0x98}

		Ha, _ := NewKeyedPRNG(key)
		Hb, _ := NewKeyedPRNG(key)

		sum0 := make([]byte, 512)
		sum1 := make([]byte, 512)

		Ha.SetClock(sum0, 256)
		Hb.SetClock(sum1, 128)

		for i := 0; i < 128; i++ {
			Hb.Clock(sum1)
		}

		Ha.Clock(sum0)
		Hb.Clock(sum1)

		require.Equal(t, sum0, sum1)
	})

}
