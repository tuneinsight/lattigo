package bootstrapper

import (
	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ckks"
)

type ParametersLiteral bootstrapping.ParametersLiteral

type Parameters struct {
	bootstrapping.Parameters
}

func NewParametersFromLiteral(paramsResidual ckks.Parameters, paramsBootstrapping ParametersLiteral) (Parameters, error) {
	params, err := bootstrapping.NewParametersFromLiteral(paramsResidual, bootstrapping.ParametersLiteral(paramsBootstrapping))
	return Parameters{
		Parameters: params,
	}, err
}
