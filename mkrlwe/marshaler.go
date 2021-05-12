package mkrlwe

/*
// MarshalBinary encodes a Ciphertext in a byte slice.
func (pk *MKPublicKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(len(ciphertext.Value))

	var pointer, inc uint64

	pointer = 1

	for _, el := range ciphertext.Value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext in the target Ciphertext.
func (pk *MKPublicKey) UnmarshalBinary(data [2][]byte) (err error) {

	bfvElement := new(rlwe.Element)

	ciphertext.Value = make([]*ring.Poly, uint8(data[0][0]))

	var pointer, inc uint64
	pointer = 1

	for i := range ciphertext.Value {

		ciphertext.Value[i] = new(ring.Poly)

		if inc, err = ciphertext.Value[i].DecodePolyNew(data[0][pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	// retrieve peer IDs

	ciphertext.PeerIDs = peers

	return nil
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (pk *MKPublicKey) GetDataLen(WithMetaData bool) (dataLen uint64) {
	if WithMetaData {
		dataLen++
	}

	for _, el := range ciphertext.Value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}

// MarshalBinary encodes a Ciphertext in a byte slice.
func (evalGal *MKEvalGalKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(len(ciphertext.Value))

	var pointer, inc uint64

	pointer = 1

	for _, el := range ciphertext.Value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext in the target Ciphertext.
func (evalGal *MKEvalGalKey) UnmarshalBinary(data [2][]byte) (err error) {

	bfvElement := new(rlwe.Element)

	ciphertext.Value = make([]*ring.Poly, uint8(data[0][0]))

	var pointer, inc uint64
	pointer = 1

	for i := range ciphertext.Value {

		ciphertext.Value[i] = new(ring.Poly)

		if inc, err = ciphertext.Value[i].DecodePolyNew(data[0][pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	// retrieve peer IDs

	ciphertext.PeerIDs = peers

	return nil
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (evalGal *MKEvalGalKey) GetDataLen(WithMetaData bool) (dataLen uint64) {
	if WithMetaData {
		dataLen++
	}

	for _, el := range ciphertext.Value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}

// MarshalBinary encodes a Ciphertext in a byte slice.
func (evalKey *MKEvaluationKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(len(ciphertext.Value))

	var pointer, inc uint64

	pointer = 1

	for _, el := range ciphertext.Value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext in the target Ciphertext.
func (evalKey *MKEvaluationKey) UnmarshalBinary(data [2][]byte) (err error) {

	bfvElement := new(rlwe.Element)

	ciphertext.Value = make([]*ring.Poly, uint8(data[0][0]))

	var pointer, inc uint64
	pointer = 1

	for i := range ciphertext.Value {

		ciphertext.Value[i] = new(ring.Poly)

		if inc, err = ciphertext.Value[i].DecodePolyNew(data[0][pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	// retrieve peer IDs

	ciphertext.PeerIDs = peers

	return nil
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (evalKey *MKEvaluationKey) GetDataLen(WithMetaData bool) (dataLen uint64) {
	if WithMetaData {
		dataLen++
	}

	for _, el := range ciphertext.Value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}
*/
