package mkbfv

// RelinearizationWithSharedRelinKey implements the algorithm 1 in section 3.3.1 of the Chen paper
// It does relin using a precomputed shared key
func RelinearizationWithSharedRelinKey(relinKey *MKRelinearizationKey, ciphertext *MKCiphertext) {

}

// RelinearizationOnTheFly implements the algorithm 2 in section 3.3.1 of the Chen paper
// It does relin directly by linearizing each entry of the extended ciphertext
func RelinearizationOnTheFly(switchingKeys []*MKSwitchingKey, ciphertext *MKCiphertext) {

}
