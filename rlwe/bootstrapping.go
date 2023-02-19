package rlwe

// Bootstrapper is a scheme independant generic interface to handle bootstrapping.
type Bootstrapper interface {

	// Bootstrapp defines a method that takes an empty interface as input and applies
	// an in place scheme specific bootstrapping.
	Bootstrap(ct interface{}) (err error)

	// The depth of the bootstrapping circuit, can be zero.
	Depth() int

	// StartingLevel defines the level at which the bootstrapping starts.
	// For the centralized bootstrapping this value is usually zero.
	// For the collective bootstrapping it is given by the user-defined
	// security parameters
	StartingLevel() int

	// EndLevel defines the level at which the bootstrapping ends.
	// For the centralized bootstrapping this value is the maximum
	// level minus the depth of the bootstrapping circuit.
	// For the collective bootstrapping this value is usually the
	// maximum level.
	EndLevel() int
}
