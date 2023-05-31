package rlwe

// PublicKey is a type for generic RLWE public keys.
// The Value field stores the polynomials in NTT and Montgomery form.
type PublicKey struct {
	OperandQP
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params ParametersInterface) (pk *PublicKey) {
	pk = &PublicKey{*NewOperandQP(params, 1, params.MaxLevelQ(), params.MaxLevelP())}
	pk.IsNTT = true
	pk.IsMontgomery = true
	return
}

func (p *PublicKey) CopyNew() *PublicKey {
	return &PublicKey{*p.OperandQP.CopyNew()}
}

func (p *PublicKey) Equal(other *PublicKey) bool {
	return p.OperandQP.Equal(&other.OperandQP)
}
