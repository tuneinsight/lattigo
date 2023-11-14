package examples

import (
	"testing"

	"github.com/tuneinsight/lattigo/v4/he/hefloat"
	"github.com/tuneinsight/lattigo/v4/he/heint"
)

func TestExampleParams(t *testing.T) {
	for _, pl := range HEIntParams {
		p, err := heint.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
	}
	for _, pl := range HEFloatParams {
		p, err := hefloat.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
	}
}
