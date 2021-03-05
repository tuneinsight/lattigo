package mkbfv

// Convert creates a switching key K_ij from the evaluation key of the i-th peer and the public key of the j-th peer.
func Convert(D *MKEvaluationKey, publicKey *MKPublicKey) *MKSwitchingKey {
	res := new(MKSwitchingKey)

	return res
}

// CreateSharedRelinearizationKey generates a shared relinearization key containing the switching key for all pair of participants.
func CreateSharedRelinearizationKey() *MKRelinearizationKey {
	res := new(MKRelinearizationKey)

	return res
}

func RelinearizationWithSharedRelinKey() {
	return
}
