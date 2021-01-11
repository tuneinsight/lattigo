package rlwe

import "github.com/ldsec/lattigo/v2/ring"

// PublicKey is a generic public key type for RLWE schemes
type PublicKey struct {
	pk [2]*ring.Poly
}

// Get returns the polynomials of the PublicKey.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the polynomial of the PublicKey as the input polynomials.
func (pk *PublicKey) Set(p [2]*ring.Poly) {
	pk.pk[0] = p[0].CopyNew()
	pk.pk[1] = p[1].CopyNew()
}

type SwitchingKey struct {
	Value [][2]*ring.Poly
}

type RotationKey struct {
}
