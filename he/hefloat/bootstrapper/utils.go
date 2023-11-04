package bootstrapper

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/core/rlwe"
	"github.com/tuneinsight/lattigo/v4/he/hefloat"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func (b Bootstrapper) SwitchRingDegreeN1ToN2New(ctN1 *rlwe.Ciphertext) (ctN2 *rlwe.Ciphertext) {
	ctN2 = hefloat.NewCiphertext(b.Parameters.Parameters.Parameters, 1, ctN1.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := b.bootstrapper.ApplyEvaluationKey(ctN1, b.evk.EvkN1ToN2, ctN2); err != nil {
		panic(err)
	}
	return
}

func (b Bootstrapper) SwitchRingDegreeN2ToN1New(ctN2 *rlwe.Ciphertext) (ctN1 *rlwe.Ciphertext) {
	ctN1 = hefloat.NewCiphertext(b.ResidualParameters, 1, ctN2.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := b.bootstrapper.ApplyEvaluationKey(ctN2, b.evk.EvkN2ToN1, ctN1); err != nil {
		panic(err)
	}
	return
}

func (b Bootstrapper) ComplexToRealNew(ctCmplx *rlwe.Ciphertext) (ctReal *rlwe.Ciphertext) {
	ctReal = hefloat.NewCiphertext(b.ResidualParameters, 1, ctCmplx.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := b.bridge.ComplexToReal(&b.bootstrapper.Evaluator.Evaluator, ctCmplx, ctReal); err != nil {
		panic(err)
	}
	return
}

func (b Bootstrapper) RealToComplexNew(ctReal *rlwe.Ciphertext) (ctCmplx *rlwe.Ciphertext) {
	ctCmplx = hefloat.NewCiphertext(b.Parameters.Parameters.Parameters, 1, ctReal.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := b.bridge.RealToComplex(&b.bootstrapper.Evaluator.Evaluator, ctReal, ctCmplx); err != nil {
		panic(err)
	}
	return
}

func (b Bootstrapper) PackAndSwitchN1ToN2(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error) {

	var err error

	if b.ResidualParameters.N() != b.Parameters.Parameters.Parameters.N() {
		if cts, err = b.Pack(cts, b.ResidualParameters, b.xPow2N1); err != nil {
			return nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN1: %w", err)
		}
	}

	for i := range cts {
		cts[i] = *b.SwitchRingDegreeN1ToN2New(&cts[i])
	}

	if cts, err = b.Pack(cts, b.Parameters.Parameters.Parameters, b.xPow2N2); err != nil {
		return nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN2: %w", err)
	}

	return cts, nil
}

func (b Bootstrapper) UnpackAndSwitchN2Tn1(cts []rlwe.Ciphertext, LogSlots, Nb int) ([]rlwe.Ciphertext, error) {

	var err error

	if cts, err = b.UnPack(cts, b.Parameters.Parameters.Parameters, LogSlots, Nb, b.xPow2InvN2); err != nil {
		return nil, fmt.Errorf("cannot UnpackAndSwitchN2Tn1: UnpackN2: %w", err)
	}

	for i := range cts {
		cts[i] = *b.SwitchRingDegreeN2ToN1New(&cts[i])
	}

	for i := range cts {
		cts[i].LogDimensions.Cols = LogSlots
	}

	return cts, nil
}

func (b Bootstrapper) UnPack(cts []rlwe.Ciphertext, params hefloat.Parameters, LogSlots, Nb int, xPow2Inv []ring.Poly) ([]rlwe.Ciphertext, error) {
	LogGap := params.LogMaxSlots() - LogSlots

	if LogGap == 0 {
		return cts, nil
	}

	cts = append(cts, make([]rlwe.Ciphertext, Nb-1)...)

	for i := 1; i < len(cts); i++ {
		cts[i] = *cts[0].CopyNew()
	}

	r := params.RingQ().AtLevel(cts[0].Level())

	N := len(cts)

	for i := 0; i < utils.Min(bits.Len64(uint64(N-1)), LogGap); i++ {

		step := 1 << (i + 1)

		for j := 0; j < N; j += step {

			for k := step >> 1; k < step; k++ {

				if (j + k) >= N {
					break
				}

				r.MulCoeffsMontgomery(cts[j+k].Value[0], xPow2Inv[i], cts[j+k].Value[0])
				r.MulCoeffsMontgomery(cts[j+k].Value[1], xPow2Inv[i], cts[j+k].Value[1])
			}
		}
	}

	return cts, nil
}

func (b Bootstrapper) Pack(cts []rlwe.Ciphertext, params hefloat.Parameters, xPow2 []ring.Poly) ([]rlwe.Ciphertext, error) {

	var LogSlots = cts[0].LogSlots()
	RingDegree := params.N()

	for i, ct := range cts {
		if N := ct.LogSlots(); N != LogSlots {
			return nil, fmt.Errorf("cannot Pack: cts[%d].PlaintextLogSlots()=%d != cts[0].PlaintextLogSlots=%d", i, N, LogSlots)
		}

		if N := ct.Value[0].N(); N != RingDegree {
			return nil, fmt.Errorf("cannot Pack: cts[%d].Value[0].N()=%d != params.N()=%d", i, N, RingDegree)
		}
	}

	LogGap := params.LogMaxSlots() - LogSlots

	if LogGap == 0 {
		return cts, nil
	}

	for i := 0; i < LogGap; i++ {

		for j := 0; j < len(cts)>>1; j++ {

			eve := cts[j*2+0]
			odd := cts[j*2+1]

			level := utils.Min(eve.Level(), odd.Level())

			r := params.RingQ().AtLevel(level)

			r.MulCoeffsMontgomeryThenAdd(odd.Value[0], xPow2[i], eve.Value[0])
			r.MulCoeffsMontgomeryThenAdd(odd.Value[1], xPow2[i], eve.Value[1])

			cts[j] = eve
		}

		if len(cts)&1 == 1 {
			cts[len(cts)>>1] = cts[len(cts)-1]
			cts = cts[:len(cts)>>1+1]
		} else {
			cts = cts[:len(cts)>>1]
		}
	}

	LogMaxDimensions := params.LogMaxDimensions()
	for i := range cts {
		cts[i].LogDimensions = LogMaxDimensions
	}

	return cts, nil
}
