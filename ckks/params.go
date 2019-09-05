package ckks

type Parameters struct {
	logN        uint64
	modulichain []uint64
	logScale    uint64
	sigma       float64
}

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security.
var DefaultParams = map[uint64]*Parameters{
	11: {11, []uint64{54}, 45, 3.2},                                                                         //logQ = 54
	12: {12, []uint64{44, 32, 32}, 32, 3.2},                                                                 //logQ = 109
	13: {13, []uint64{49, 42, 42, 42, 42}, 42, 3.2},                                                         //logQ = 218
	14: {14, []uint64{50, 43, 43, 43, 43, 43, 43, 43, 43, 43}, 43, 3.2},                                     //logQ = 438
	15: {15, []uint64{53, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46}, 46, 3.2}, //logQ = 881
}
