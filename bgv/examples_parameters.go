package bgv

var (
	// ExampleParameters128BitLogN14LogQP438 is an example parameters set with logN=14, logQP=438
	// and a 16-bit plaintext modulus, offering 128-bit of security.
	ExampleParameters128BitLogN14LogQP438 = ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x10000048001, 0x20008001, 0x1ffc8001,
			0x20040001, 0x1ffc0001, 0x1ffb0001,
			0x20068001, 0x1ff60001, 0x200b0001,
			0x200d0001, 0x1ff18001, 0x200f8001}, // 40 + 11*29 bits
		P:                []uint64{0x10000140001, 0x7ffffb0001}, // 40 + 39 bits
		PlaintextModulus: 0x10001,                               // 16 bits
	}
)
