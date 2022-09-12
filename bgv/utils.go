package bgv

import (
	"math/big"
)

// Norm returns the log2 of the standard deviation, minimum and maximum absolute norm of
// the decrypted ciphertext, before the decoding (i.e. including the error).
func Norm(ct *Ciphertext, dec Decryptor) (std, min, max float64) {

	params := dec.(*decryptor).params

	coeffsBigint := make([]*big.Int, params.N())
	for i := range coeffsBigint {
		coeffsBigint[i] = new(big.Int)
	}

	buffQ := params.RingQ().NewPoly()
	pt := NewPlaintextAtLevelFromPoly(ct.Level(), buffQ).Plaintext
	dec.(*decryptor).Decryptor.Decrypt(ct.El(), pt)
	params.RingQ().InvNTTLvl(ct.Level(), buffQ, buffQ)
	params.RingQ().PolyToBigintCenteredLvl(ct.Level(), buffQ, 1, coeffsBigint)

	return errorStats(coeffsBigint)
}
