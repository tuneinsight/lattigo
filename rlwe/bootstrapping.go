package rlwe

// Bootstrapper is a scheme independent generic interface to handle bootstrapping.
type Bootstrapper interface {

	// Bootstrap defines a method that takes a single Ciphertext as input and applies
	// an in place scheme-specific bootstrapping. The result is also returned.
	// An error should notably be returned if ct.Level() < MinimumInputLevel().
	Bootstrap(ct *Ciphertext) (*Ciphertext, error)

	// BootstrapMany defines a method that takes a slice of Ciphertexts as input and applies an
	// in place scheme-specific bootstrapping to each Ciphertext. The result is also returned.
	// An error should notably be returned if ct.Level() < MinimumInputLevel().
	BootstrapMany(cts []Ciphertext) ([]Ciphertext, error)

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
