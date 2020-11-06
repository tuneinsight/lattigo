package ring

// Parameters is a struct storing test parameters for the package Ring.
type Parameters struct {
	logN uint64
	qi   []uint64
	pi   []uint64
}

// DefaultParams is a struct storing default test parameters of the Qi and Pi moduli for the package Ring.
var DefaultParams = []*Parameters{
	{12, Qi60[len(Qi60)-2:], Pi60[len(Pi60)-2:]},
	{13, Qi60[len(Qi60)-4:], Pi60[len(Pi60)-4:]},
	{14, Qi60[len(Qi60)-7:], Pi60[len(Pi60)-7:]},
	{15, Qi60[len(Qi60)-14:], Pi60[len(Pi60)-14:]},
	{16, Qi60[len(Qi60)-29:], Pi60[len(Pi60)-29:]},
}

// Qi60 are the first [0:32] 61-bit close to 2^{62} NTT-friendly primes for N up to 2^{17}
var Qi60 = []uint64{0x1fffffffffe00001, 0x1fffffffffc80001, 0x1fffffffffb40001, 0x1fffffffff500001,
	0x1fffffffff380001, 0x1fffffffff000001, 0x1ffffffffef00001, 0x1ffffffffee80001,
	0x1ffffffffeb40001, 0x1ffffffffe780001, 0x1ffffffffe600001, 0x1ffffffffe4c0001,
	0x1ffffffffdf40001, 0x1ffffffffdac0001, 0x1ffffffffda40001, 0x1ffffffffc680001,
	0x1ffffffffc000001, 0x1ffffffffb880001, 0x1ffffffffb7c0001, 0x1ffffffffb300001,
	0x1ffffffffb1c0001, 0x1ffffffffadc0001, 0x1ffffffffa400001, 0x1ffffffffa140001,
	0x1ffffffff9d80001, 0x1ffffffff9140001, 0x1ffffffff8ac0001, 0x1ffffffff8a80001,
	0x1ffffffff81c0001, 0x1ffffffff7800001, 0x1ffffffff7680001, 0x1ffffffff7080001}

// Pi60 are the next [32:64] 61-bit close to 2^{62} NTT-friendly primes for N up to 2^{17}
var Pi60 = []uint64{0x1ffffffff6c80001, 0x1ffffffff6140001, 0x1ffffffff5f40001, 0x1ffffffff5700001,
	0x1ffffffff4bc0001, 0x1ffffffff4380001, 0x1ffffffff3240001, 0x1ffffffff2dc0001,
	0x1ffffffff1a40001, 0x1ffffffff11c0001, 0x1ffffffff0fc0001, 0x1ffffffff0d80001,
	0x1ffffffff0c80001, 0x1ffffffff08c0001, 0x1fffffffefd00001, 0x1fffffffef9c0001,
	0x1fffffffef600001, 0x1fffffffeef40001, 0x1fffffffeed40001, 0x1fffffffeed00001,
	0x1fffffffeebc0001, 0x1fffffffed540001, 0x1fffffffed440001, 0x1fffffffed2c0001,
	0x1fffffffed200001, 0x1fffffffec940001, 0x1fffffffec6c0001, 0x1fffffffebe80001,
	0x1fffffffebac0001, 0x1fffffffeba40001, 0x1fffffffeb4c0001, 0x1fffffffeb280001}
