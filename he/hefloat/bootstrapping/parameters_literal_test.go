package bootstrapping

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils"
)

func TestParameterLiteralMarshal(t *testing.T) {
	schemaParamsLit := testPrec45
	btpParamsLit := &ParametersLiteral{}

	schemaParamsLit.LogN = 16
	btpParamsLit.EphemeralSecretWeight = utils.Pointy(DefaultEphemeralSecretWeight)

	params, err := hefloat.NewParametersFromLiteral(schemaParamsLit)
	require.NoError(t, err)

	btpParams, err := NewParametersFromLiteral(params, *btpParamsLit)
	require.NoError(t, err)

	binParams, err := btpParams.MarshalBinary()
	require.NoError(t, err)

	paramsLit2 := &ParametersLiteral{}

	err = paramsLit2.UnmarshalBinary(binParams)
	require.NoError(t, err)

	require.Equal(t, btpParamsLit.LogN, paramsLit2.LogN)
	require.Equal(t, btpParamsLit.LogP, paramsLit2.LogP)
	require.Equal(t, btpParamsLit.LogSlots, paramsLit2.LogSlots)
	require.Equal(t, btpParamsLit.K, paramsLit2.K)
	require.Equal(t, btpParamsLit.Xs, paramsLit2.Xs)
	require.Equal(t, btpParamsLit.Xe, paramsLit2.Xe)
	require.Equal(t, btpParamsLit.CoeffsToSlotsFactorizationDepthAndLogScales, paramsLit2.CoeffsToSlotsFactorizationDepthAndLogScales)
	require.Equal(t, btpParamsLit.SlotsToCoeffsFactorizationDepthAndLogScales, paramsLit2.SlotsToCoeffsFactorizationDepthAndLogScales)
	require.Equal(t, btpParamsLit.EvalModLogScale, paramsLit2.EvalModLogScale)
	require.Equal(t, btpParamsLit.EphemeralSecretWeight, paramsLit2.EphemeralSecretWeight)
	require.Equal(t, btpParamsLit.IterationsParameters, paramsLit2.IterationsParameters)
	require.Equal(t, btpParamsLit.Mod1Type, paramsLit2.Mod1Type)
	require.Equal(t, btpParamsLit.LogMessageRatio, paramsLit2.LogMessageRatio)
	require.Equal(t, btpParamsLit.Mod1Degree, paramsLit2.Mod1Degree)
	require.Equal(t, btpParamsLit.DoubleAngle, paramsLit2.DoubleAngle)
	require.Equal(t, btpParamsLit.Mod1InvDegree, paramsLit2.Mod1InvDegree)
}

func TestParameterMarshal(t *testing.T) {

	schemaParamsLit := testPrec45
	btpParamsLit := &ParametersLiteral{}

	schemaParamsLit.LogN = 16

	params, err := hefloat.NewParametersFromLiteral(schemaParamsLit)
	require.NoError(t, err)

	btpParams, err := NewParametersFromLiteral(params, *btpParamsLit)
	require.NoError(t, err)

	binParams, err := btpParams.MarshalBinary()
	require.NoError(t, err)

	params2 := &Parameters{}

	err = params2.UnmarshalBinary(binParams)
	require.NoError(t, err)

	require.True(t, btpParams.Equal(params2))
}
