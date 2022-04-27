package bootstrapping

import (
	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ckks/advanced"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Parameters is a struct for the default bootstrapping parameters
type Parameters struct {
	SlotsToCoeffsParameters advanced.EncodingMatrixLiteral
	EvalModParameters       advanced.EvalModLiteral
	CoeffsToSlotsParameters advanced.EncodingMatrixLiteral
	EphemeralSecretWeight   int // Hamming weight of the ephemeral secret. If 0, no ephemeral secret is used during the bootstrapping.
}

// MarshalBinary encode the target Parameters on a slice of bytes.
func (p *Parameters) MarshalBinary() (data []byte, err error) {
	data = []byte{}

	var tmp []byte
	if tmp, err = p.SlotsToCoeffsParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	if tmp, err = p.EvalModParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	if tmp, err = p.CoeffsToSlotsParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	tmp = make([]byte, 4)
	tmp[0] = uint8(p.EphemeralSecretWeight >> 24)
	tmp[1] = uint8(p.EphemeralSecretWeight >> 16)
	tmp[2] = uint8(p.EphemeralSecretWeight >> 8)
	tmp[3] = uint8(p.EphemeralSecretWeight >> 0)
	data = append(data, tmp...)
	return
}

// UnmarshalBinary decodes a slice of bytes on the target Parameters.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {

	pt := 0
	dLen := int(data[pt])

	if err := p.SlotsToCoeffsParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++
	dLen = int(data[pt])

	if err := p.EvalModParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++
	dLen = int(data[pt])

	if err := p.CoeffsToSlotsParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++

	p.EphemeralSecretWeight = int(data[pt])<<24 | int(data[pt+1])<<16 | int(data[pt+2])<<8 | int(data[pt+3])

	return
}

// RotationsForBootstrapping returns the list of rotations performed during the Bootstrapping operation.
func (p *Parameters) RotationsForBootstrapping(params ckks.Parameters) (rotations []int) {

	logN := params.LogN()
	logSlots := params.LogSlots()

	// List of the rotation key values to needed for the bootstrapp
	rotations = []int{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < logN-1; i++ {
		if !utils.IsInSliceInt(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	p.CoeffsToSlotsParameters.LogN = logN
	p.SlotsToCoeffsParameters.LogN = logN

	p.CoeffsToSlotsParameters.LogSlots = logSlots
	p.SlotsToCoeffsParameters.LogSlots = logSlots

	rotations = append(rotations, p.CoeffsToSlotsParameters.Rotations()...)
	rotations = append(rotations, p.SlotsToCoeffsParameters.Rotations()...)

	return
}
