package bootstrapper

import (
	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ckks"
)

// ParametersLiteral is a wrapper of bootstrapping.ParametersLiteral.
// See bootstrapping.ParametersLiteral for additional information.
type ParametersLiteral bootstrapping.ParametersLiteral

// Parameters is a wrapper of the bootstrapping.Parameters.
// See bootstrapping.Parameters for additional information.
type Parameters struct {
	bootstrapping.Parameters
	ResidualParameters ckks.Parameters
}

// NewParametersFromLiteral is a wrapper of bootstrapping.NewParametersFromLiteral.
// See bootstrapping.NewParametersFromLiteral for additional information.
func NewParametersFromLiteral(paramsResidual ckks.Parameters, paramsBootstrapping ParametersLiteral) (Parameters, error) {
	params, err := bootstrapping.NewParametersFromLiteral(paramsResidual, bootstrapping.ParametersLiteral(paramsBootstrapping))
	return Parameters{
		Parameters:         params,
		ResidualParameters: paramsResidual,
	}, err
}
