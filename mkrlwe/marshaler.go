package mkrlwe

/*
// MarshalBinary encodes a Galois evaluation key in a byte slice.
func (pubKey *MKPublicKey) MarshalBinary() (data [3][]byte, err error) {

	var err1, err2 error

	data[0], err1 = pubKey.Key[0].MarshalBinary()
	data[1], err2 = pubKey.Key[1].MarshalBinary()
	data[2] = make([]byte, 8)
	binary.LittleEndian.PutUint64(data[2], pubKey.PeerID)

	if err1 != nil || err2 != nil {
		panic("Couldn't Marshall evaluation key")
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Galois Evaluation Key in the target MKEvalGalKey.
func (pubKey *MKPublicKey) UnmarshalBinary(data [3][]byte) (err error) {

	pubKey.Key[0] = new(MKDecomposedPoly)
	pubKey.Key[1] = new(MKDecomposedPoly)

	err1 := pubKey.Key[0].UnmarshalBinary(data[0])
	err2 := pubKey.Key[1].UnmarshalBinary(data[1])
	pubKey.PeerID = binary.LittleEndian.Uint64(data[2])

	if err1 != nil || err2 != nil {
		panic("Couldn't UnMarshall evaluation key")
	}

	return nil
}

// MarshalBinary encodes a Galois evaluation key in a byte slice.
func (evalGal *MKEvalGalKey) MarshalBinary() (data [3][]byte, err error) {

	var err1, err2 error

	data[0], err1 = evalGal.Key[0].MarshalBinary()
	data[1], err2 = evalGal.Key[1].MarshalBinary()
	data[2] = make([]byte, 8)
	binary.LittleEndian.PutUint64(data[2], evalGal.PeerID)

	if err1 != nil || err2 != nil {
		panic("Couldn't Marshall evaluation key")
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Galois Evaluation Key in the target MKEvalGalKey.
func (evalGal *MKEvalGalKey) UnmarshalBinary(data [3][]byte) (err error) {

	evalGal.Key[0] = new(MKDecomposedPoly)
	evalGal.Key[1] = new(MKDecomposedPoly)

	err1 := evalGal.Key[0].UnmarshalBinary(data[0])
	err2 := evalGal.Key[1].UnmarshalBinary(data[1])
	evalGal.PeerID = binary.LittleEndian.Uint64(data[2])

	if err1 != nil || err2 != nil {
		panic("Couldn't UnMarshall evaluation key")
	}

	return nil
}

// MarshalBinary encodes an evaluation key in a byte slice.
func (evalKey *MKEvaluationKey) MarshalBinary() (data [4][]byte, err error) {

	var err1, err2, err3 error

	data[0], err1 = evalKey.Key[0].MarshalBinary()
	data[1], err2 = evalKey.Key[0].MarshalBinary()
	data[2], err3 = evalKey.Key[0].MarshalBinary()
	data[3] = make([]byte, 8)
	binary.LittleEndian.PutUint64(data[3], evalKey.PeerID)

	if err1 != nil || err2 != nil || err3 != nil {
		panic("Couldn't Marshall evaluation key")
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled evaluation key in the target MKEvaluationKey.
func (evalKey *MKEvaluationKey) UnmarshalBinary(data [4][]byte) (err error) {

	evalKey.Key[0] = new(MKDecomposedPoly)
	evalKey.Key[1] = new(MKDecomposedPoly)
	evalKey.Key[2] = new(MKDecomposedPoly)

	err1 := evalKey.Key[0].UnmarshalBinary(data[0])
	err2 := evalKey.Key[1].UnmarshalBinary(data[1])
	err3 := evalKey.Key[2].UnmarshalBinary(data[2])
	evalKey.PeerID = binary.LittleEndian.Uint64(data[3])

	if err1 != nil || err2 != nil || err3 != nil {
		panic("Couldn't UnMarshall evaluation key")
	}

	return nil
}

// MarshalBinary encodes a decomposed poly in a byte slice.
func (decPol *MKDecomposedPoly) MarshalBinary() (data []byte, err error) {

	data = make([]byte, decPol.GetDataLen(true))

	data[0] = uint8(len(decPol.Poly))

	var pointer, inc uint64

	pointer = 1

	for _, p := range decPol.Poly {

		if inc, err = p.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled MKDecomposedPoly in the target MKDecomposedPoly.
func (decPol *MKDecomposedPoly) UnmarshalBinary(data []byte) (err error) {

	decPol.Poly = make([]*ring.Poly, uint8(data[0]))

	var pointer, inc uint64
	pointer = 1

	for i := range decPol.Poly {

		decPol.Poly[i] = new(ring.Poly)

		if inc, err = decPol.Poly[i].DecodePolyNew(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	return nil
}

// GetDataLen returns the length in bytes of the target MKDecomposedPoly
func (decPol *MKDecomposedPoly) GetDataLen(WithMetaData bool) (dataLen uint64) {
	if WithMetaData {
		dataLen++
	}

	for _, p := range decPol.Poly {
		dataLen += p.GetDataLen(WithMetaData)
	}

	return dataLen
}

// MarshallPeerID transforms a list of uint64 in an array of bytes
func MarshallPeerID(peerID []uint64) (data []byte) {

	data = make([]byte, 8*len(peerID))

	for i, p := range peerID {
		binary.LittleEndian.PutUint64(data[8*i:8*(i+1)], p)
	}
	return
}

// UnMarshallPeerID retrieve an array of uint representing the peers from an array of bytes
func UnMarshallPeerID(data []byte) []uint64 {

	res := make([]uint64, len(data)/8)

	for i := range res {
		res[i] = binary.LittleEndian.Uint64(data[8*i : 8*(i+1)])
	}

	return res
}
*/
