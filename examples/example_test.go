package examples

import (
	"testing"

	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/he/heint"
)

func TestExampleParams(t *testing.T) {
	for _, pl := range HEIntParams {
		p, err := heint.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
		t.Logf("HEIntParams: LogN: %d - LogQP: %12.7f - LogSlots: %d", p.LogN(), p.LogQP(), p.LogMaxSlots())
	}

	for _, pl := range HEIntScaleInvariantParams {
		p, err := heint.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
		t.Logf("HEIntScaleInvariantParams: LogN: %d - LogQP: %12.7f - LogSlots: %d", p.LogN(), p.LogQP(), p.LogMaxSlots())
	}

	for _, pl := range HEFloatComplexParams {
		p, err := hefloat.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
		t.Logf("HEFloatComplex: LogN: %d - LogQP: %12.7f - LogSlots: %d", p.LogN(), p.LogQP(), p.LogMaxSlots())
	}

	for _, pl := range HEFloatRealParams {
		p, err := hefloat.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
		t.Logf("HEFloatReal: LogN: %d - LogQP: %12.7f - LogSlots: %d", p.LogN(), p.LogQP(), p.LogMaxSlots())
	}
}
