package bfv

var (
	// ExampleParameters128BitLogN14LogQP438 is an example parameters set with logN=14, logQP=438
	// and a 16-bit plaintext modulus, offering 128-bit of security.
	ExampleParameters128BitLogN14LogQP438 = ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
			0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
		P:                []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
		PlaintextModulus: 0x10001,
	}
)
