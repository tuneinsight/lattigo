package bfv

import (
	"github.com/stretchr/testify/assert"
	"testing"
)

func TestParams_BinaryMarshaller(t *testing.T) {
	t.Run("ZeroValue", func(t *testing.T) {
		bytes, err := (&Parameters{}).MarshalBinary()
		assert.Nil(t, err)
		assert.Equal(t, []byte{}, bytes)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		assert.NotNil(t, err)
	})
	t.Run("SupportedParams", func(t *testing.T) {
		for _, params := range DefaultParams {
			bytes, err := params.MarshalBinary()
			assert.Nil(t, err)
			var p Parameters
			err = p.UnmarshalBinary(bytes)
			assert.Nil(t, err)
			assert.Equal(t, params, p)
		}
	})
}
