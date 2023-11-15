package he

// Bootstrapper is a scheme-independent generic interface to handle bootstrapping.
type Bootstrapper[CiphertextType any] interface {

	// Bootstrap defines a method that takes a single Ciphertext as input and applies
	// an in place scheme-specific bootstrapping. The result is also returned.
	// An error should notably be returned if ct.Level() < Bootstrapper.MinimumInputLevel().
	Bootstrap(ct *CiphertextType) (*CiphertextType, error)

	// BootstrapMany defines a method that takes a slice of Ciphertexts as input and applies an
	// in place scheme-specific bootstrapping to each Ciphertext. The result is also returned.
	// An error should notably be returned if cts[i].Level() < Bootstrapper.MinimumInputLevel().
	BootstrapMany(cts []CiphertextType) ([]CiphertextType, error)

	// Depth is the number of levels consumed by the bootstrapping circuit.
	// This value is equivalent to params.MaxLevel() - OutputLevel().
	Depth() int

	// MinimumInputLevel defines the minimum level that the ciphertext
	// must be at when given to the bootstrapper.
	// For the centralized bootstrapping this value is usually zero.
	// For the collective bootstrapping it is given by the user-defined
	// security parameters
	MinimumInputLevel() int

	// OutputLevel defines the level that the ciphertext will be at
	// after the bootstrapping.
	// For the centralized bootstrapping this value is the maximum
	// level minus the depth of the bootstrapping circuit.
	// For the collective bootstrapping this value is usually the
	// maximum level.
	OutputLevel() int
}
