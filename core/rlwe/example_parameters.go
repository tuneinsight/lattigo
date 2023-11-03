package rlwe

var (
	// ExampleParameterLogN14LogQP438 is an example parameters set with logN=14 and logQP=438
	// offering 128-bit of security.
	ExampleParametersLogN14LogQP438 = ParametersLiteral{
		LogN:    14,
		Q:       []uint64{0x200000440001, 0x7fff80001, 0x800280001, 0x7ffd80001, 0x7ffc80001},
		P:       []uint64{0x3ffffffb80001, 0x4000000800001},
		NTTFlag: true,
	}
)
