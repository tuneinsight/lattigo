package ckks

import (
	"testing"

	"github.com/stretchr/testify/require"
)

func TestParams_BinaryMarshaller(t *testing.T) {
	t.Run("ZeroValue", func(t *testing.T) {
		bytes, err := (&Parameters{}).MarshalBinary()

		require.Nil(t, err)
		require.Equal(t, []byte{}, bytes)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		require.NotNil(t, err)
	})
	t.Run("SupportedParams", func(t *testing.T) {
		for _, params := range DefaultParams {
			bytes, err := params.MarshalBinary()
			require.Nil(t, err)
			var p Parameters
			err = p.UnmarshalBinary(bytes)
			require.Nil(t, err)
			require.Equal(t, params, &p)
		}
	})
}
