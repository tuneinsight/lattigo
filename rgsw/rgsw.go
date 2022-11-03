// Package rgsw implements an RLWE-based RGSW encryption scheme. In RSGW, ciphertexts are tuples of two gadget ciphertexts
// where the first gadget ciphertext encrypts the message and the second gadget ciphertext encrypts the message times the
// secret. This package only implements a subset of the RGSW scheme that is necessary for bridging between RLWE and LWE-based
// schemes and for supporting look-up table evaluation.
package rgsw
