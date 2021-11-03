package lwe


type Plaintext struct{
	Value []uint64
}

type Ciphertext struct{
	Value [][]uint64
}

func (pt *Plaintext) Level() int {
	return len(pt.Value)-1
}

func (ct *Ciphertext) Level() int {
	return len(ct.Value)-1
}

func NewCiphertext(N, level int) (ct *Ciphertext){
	ct = new(Ciphertext)
	ct.Value = make([][]uint64, level+1)
	for i := 0; i < level+1; i++{
		ct.Value[i] = make([]uint64, N+1)
	}
	return ct
}