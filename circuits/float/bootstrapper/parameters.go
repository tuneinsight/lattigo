package bootstrapper

import (
	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type ParametersLiteral bootstrapping.ParametersLiteral

type Parameters struct {
	bootstrapping.Parameters
	RingType ring.Type
}

func NewParametersFromLiteral(paramsResidual ckks.Parameters, paramsBootstrapping ParametersLiteral) (Parameters, error) {
	params, err := bootstrapping.NewParametersFromLiteral(paramsResidual, bootstrapping.ParametersLiteral(paramsBootstrapping))
	return Parameters{
		Parameters: params,
		RingType:   paramsResidual.RingType(),
	}, err
}
