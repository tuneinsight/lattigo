package examples

import (
	"testing"

	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestExampleParams(t *testing.T) {
	for _, pl := range BGVParams {
		p, err := bgv.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
		t.Logf("BGVParams: LogN: %d - LogQP: %12.7f - LogSlots: %d", p.LogN(), p.LogQP(), p.LogMaxSlots())
	}

	for _, pl := range BGVScaleInvariantParams {
		p, err := bgv.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
		t.Logf("BGVScaleInvariantParams: LogN: %d - LogQP: %12.7f - LogSlots: %d", p.LogN(), p.LogQP(), p.LogMaxSlots())
	}

	for _, pl := range CKKSComplexParams {
		p, err := ckks.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
		t.Logf("CKKSComplex: LogN: %d - LogQP: %12.7f - LogSlots: %d", p.LogN(), p.LogQP(), p.LogMaxSlots())
	}

	for _, pl := range CKKSRealParams {
		p, err := ckks.NewParametersFromLiteral(pl)
		if err != nil {
			t.Fatal(err)
		}
		p.RingQ()
		t.Logf("CKKSReal: LogN: %d - LogQP: %12.7f - LogSlots: %d", p.LogN(), p.LogQP(), p.LogMaxSlots())
	}
}
