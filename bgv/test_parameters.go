package bgv

var (
	// TESTN13QP218 is a of 128-bit secure test parameters set with a 32-bit plaintext and depth 4.
	testInsecure = ParametersLiteral{
		LogN: 10,
		Q:    []uint64{0x3fffffa8001, 0x1000090001, 0x10000c8001, 0x10000f0001, 0xffff00001},
		P:    []uint64{0x7fffffd8001},
	}

	testPlaintextModulus = []uint64{0x101, 0xffc001}

	// TestParams is a set of test parameters for BGV ensuring 128 bit security in the classic setting.
	testParams = []ParametersLiteral{testInsecure}
)
