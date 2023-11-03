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
//
// >>>>>>>!WARNING!<<<<<<<
// The bootstrapping parameters use their own and independent cryptographic parameters (i.e. ckks.Parameters)
// which are instantiated based on the option specified in `paramsBootstrapping` (and the default values of
// bootstrapping.Parameters).
// It is the user's responsibility to ensure that these scheme parameters meet the target security and to tweak them
// if necessary.
// It is possible to access information about these cryptographic parameters directly through the
// instantiated bootstrapper.Parameters struct which supports and API an identical to the ckks.Parameters.
func NewParametersFromLiteral(paramsResidual ckks.Parameters, paramsBootstrapping ParametersLiteral) (Parameters, error) {
	params, err := bootstrapping.NewParametersFromLiteral(paramsResidual, bootstrapping.ParametersLiteral(paramsBootstrapping))
	return Parameters{
		Parameters:         params,
		ResidualParameters: paramsResidual,
	}, err
}
