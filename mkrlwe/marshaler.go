package mkrlwe

import (
	"encoding/binary"

	"github.com/ldsec/lattigo/v2/rlwe"

	"github.com/ldsec/lattigo/v2/ring"
)

// MarshalBinary encodes a MKPublic key in a byte slice.
func (pubKey *MKPublicKey) MarshalBinary() (data []byte, err error) {

	d1, err1 := pubKey.Key[0].MarshalBinary()
	d2, err2 := pubKey.Key[1].MarshalBinary()
	d3 := make([]byte, 8)
	binary.LittleEndian.PutUint64(d3, pubKey.PeerID)

	if err1 != nil || err2 != nil {
		panic("Couldn't Marshall public key")
	}

	// Merge the byte slices
	d12 := MergeTwoByteSlices(d1, d2)
	data = make([]byte, len(d12)+8)
	for i := 0; i < len(d12); i++ {
		data[i] = d12[i]
	}
	for i := 0; i < len(d3); i++ {
		data[len(d12)+i] = d3[i]
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Galois Evaluation Key in the target MKEvalGalKey.
func (pubKey *MKPublicKey) UnmarshalBinary(data []byte) (err error) {

	pubKey.Key[0] = new(MKDecomposedPoly)
	pubKey.Key[1] = new(MKDecomposedPoly)

	separator := binary.LittleEndian.Uint64(data[:8])
	err1 := pubKey.Key[0].UnmarshalBinary(data[8 : 8+separator])
	err2 := pubKey.Key[1].UnmarshalBinary(data[8+separator : len(data)-8])

	pubKey.PeerID = binary.LittleEndian.Uint64(data[len(data)-8:])

	if err1 != nil || err2 != nil {
		panic("Couldn't UnMarshall evaluation key")
	}

	return nil
}

// MarshalBinary encodes a Galois evaluation key in a byte slice.
func (evalGal *MKRotationKey) MarshalBinary() (data []byte, err error) {

	d1, err1 := MarshalSwitchingKey(evalGal.Key)
	d2 := make([]byte, 8)
	binary.LittleEndian.PutUint64(d2, evalGal.PeerID)

	if err1 != nil {
		panic("Couldn't Marshall evaluation key")
	}
	data = make([]byte, len(d1)+8)
	for i := 0; i < len(d1); i++ {
		data[i] = d1[i]
	}
	for i := 0; i < len(d2); i++ {
		data[len(d1)+i] = d2[i]
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Galois Evaluation Key in the target MKEvalGalKey.
func (evalGal *MKRotationKey) UnmarshalBinary(data []byte) (err error) {

	evalKey, err := UnmarshalSwitchingKey(data[:len(data)-8])
	pid := binary.LittleEndian.Uint64(data[len(data)-8:])

	if err != nil {
		panic("Couldn't UnMarshall evaluation key")
	}
	evalGal.Key = &evalKey
	evalGal.PeerID = pid

	return nil
}

// MarshalBinary encodes an evaluation key in a byte slice.
func (evalKey *MKRelinearizationKey) MarshalBinary() (data []byte, err error) {

	// marshal each component of evalKey
	d1, err1 := MarshalSwitchingKey(evalKey.Key01)
	d2, err2 := evalKey.Key2.MarshalBinary()
	d3 := make([]byte, 8)
	binary.LittleEndian.PutUint64(d3, evalKey.PeerID)
	if err1 != nil || err2 != nil {
		panic("Couldn't Marshall evaluation key")
	}

	// Merge the byte slices
	d12 := MergeTwoByteSlices(d1, d2)
	data = make([]byte, len(d12)+8)
	for i := 0; i < len(d12); i++ {
		data[i] = d12[i]
	}
	for i := 0; i < len(d3); i++ {
		data[len(d12)+i] = d3[i]
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled evaluation key in the target MKEvaluationKey.
func (evalKey *MKRelinearizationKey) UnmarshalBinary(data []byte) (err error) {

	// Unmarshal the switchingKey first
	separator := binary.LittleEndian.Uint64(data[:8])
	swkSize := binary.LittleEndian.Uint64(data[8:16])
	value := make([][2]*ring.Poly, swkSize/2)

	swk := new(rlwe.SwitchingKey)
	var pointer, inc uint64
	pointer = 16
	dpStartIndex := separator + 8

	for j := uint64(0); j < swkSize/2; j++ {
		for k := 0; k < 2; k++ {
			value[j][k] = new(ring.Poly)
			if inc, err = value[j][k].DecodePolyNew(data[pointer:dpStartIndex]); err != nil {
				return err
			}
			pointer += inc

		}
	}
	swk.Value = value
	evalKey.Key01 = swk

	// Unmarshal the decomposedPoly
	k2 := new(MKDecomposedPoly)
	dpSize := binary.LittleEndian.Uint64(data[dpStartIndex : dpStartIndex+8])
	dpoly := make([]*ring.Poly, dpSize)
	pointer = dpStartIndex

	for i := range dpoly {

		dpoly[i] = new(ring.Poly)

		if inc, err = dpoly[i].DecodePolyNew(data[pointer : len(data)-8]); err != nil {
			return err
		}

		pointer += inc
	}
	k2.Poly = dpoly
	evalKey.Key2 = k2

	// Unmarshal the peerId
	pid := binary.LittleEndian.Uint64(data[len(data)-8:])
	evalKey.PeerID = pid

	return nil

}

// MarshalSwitchingKey encodes a switching in a byte slice.
func MarshalSwitchingKey(swk *rlwe.SwitchingKey) (data []byte, err error) {
	dataLen, numEl := GetSwkDataLen(swk, true)
	data = make([]byte, dataLen+7)
	binary.LittleEndian.PutUint64(data, uint64(numEl))
	var pointer, inc uint64

	pointer = 8
	for _, p := range swk.Value {
		for _, pol := range p {
			if inc, err = pol.WriteTo(data[pointer:]); err != nil {
				return nil, err
			}

			pointer += inc
		}
	}

	return data, nil
}

// UnmarshalSwitchingKey encodes a switching in a byte slice.
func UnmarshalSwitchingKey(data []byte) (rlwe.SwitchingKey, error) {

	size := binary.LittleEndian.Uint64(data[:8])
	value := make([][2]*ring.Poly, size/2)

	swk := new(rlwe.SwitchingKey)
	var pointer, inc uint64
	var err error
	pointer = 8
	for j := uint64(0); j < size/2; j++ {
		for k := 0; k < 2; k++ {
			value[j][k] = new(ring.Poly)
			if inc, err = value[j][k].DecodePolyNew(data[pointer:]); err != nil {
				return *swk, err
			}
			pointer += inc

		}
	}
	swk.Value = value
	return *swk, err

}

// GetSwkDataLen returns the length in bytes of the switching key
func GetSwkDataLen(swk *rlwe.SwitchingKey, WithMetaData bool) (dataLen uint64, numEl int) {
	if WithMetaData {
		dataLen++
	}

	for _, p := range swk.Value {
		for _, pol := range p {
			numEl++
			dataLen += pol.GetDataLen(WithMetaData)
		}
	}

	return dataLen, numEl
}

// MergeTwoByteSlices takes two slices of bytes resulting from marshalling, and places the separator in the first bytes
func MergeTwoByteSlices(d1 []byte, d2 []byte) (data []byte) {
	data = make([]byte, len(d1)+len(d2)+8)
	binary.LittleEndian.PutUint64(data, uint64(len(d1)))

	for i := 0; i < len(d1); i++ {
		data[i+8] = d1[i]
	}
	offset := len(d1) + 8
	for i := 0; i < len(d2); i++ {
		data[offset+i] = d2[i]
	}
	return data
}

// MergeByteSlicesWithoutSum takes two slices of bytes resulting from marshalling, and merges the two while retaining only the length of the first one
func MergeByteSlicesWithoutSum(d1 []byte, d2 []byte) (data []byte) {
	data = make([]byte, len(d1)+len(d2)-1)

	for i := 0; i < len(d1); i++ {
		data[i] = d1[i]
	}
	for i := len(d1); i < len(d1)+len(d2); i++ {
		data[i] = d2[i-len(d1)]
	}
	return data
}

// MarshalBinary encodes a decomposed poly in a byte slice.
func (decPol *MKDecomposedPoly) MarshalBinary() (data []byte, err error) {

	data = make([]byte, decPol.GetDataLen(true)+7)
	binary.LittleEndian.PutUint64(data, uint64(len(decPol.Poly)))

	var pointer, inc uint64

	pointer = 8
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

	size := binary.LittleEndian.Uint64(data[:8])
	decPol.Poly = make([]*ring.Poly, size)

	var pointer, inc uint64
	pointer = 8
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
