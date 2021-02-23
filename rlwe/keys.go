package rlwe

import (
	"encoding/binary"

	"github.com/ldsec/lattigo/v2/ring"
)

// SecretKey is a structure that stores the SecretKey.
type SecretKey struct {
	Value *ring.Poly
}

// PublicKey is a structure that stores the PublicKey.
type PublicKey struct {
	Value [2]*ring.Poly
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	Value [][2]*ring.Poly
}

// RelinearizationKey is a structure that stores the switching-keys required during the relinearization.
type RelinearizationKey struct {
	Keys []*SwitchingKey
}

// RotationKeySet is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeySet struct {
	Keys map[uint64]*SwitchingKey
}

// EvaluationKey is a structure representing the complete set of switching-keys required for evaluation
type EvaluationKey struct {
	Rlk  *RelinearizationKey
	Rtks *RotationKeySet
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(ringDegree, moduliCount uint64) (pk *PublicKey) {
	return &PublicKey{Value: [2]*ring.Poly{ring.NewPoly(ringDegree, moduliCount), ring.NewPoly(ringDegree, moduliCount)}}
}

// Get returns the polynomials of the PublicKey.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.Value
}

// Set sets the polynomial of the PublicKey as the input polynomials.
func (pk *PublicKey) Set(p [2]*ring.Poly) {
	pk.Value[0] = p[0].CopyNew()
	pk.Value[1] = p[1].CopyNew()
}

// NewRotationKeySet returns a new RotationKeySet with pre-allocated switching keys for each distinct galoisElement value.
func NewRotationKeySet(galoisElement []uint64, ringDegree, moduliCount, decompSize uint64) (rotKey *RotationKeySet) {
	rotKey = new(RotationKeySet)
	rotKey.Keys = make(map[uint64]*SwitchingKey, len(galoisElement))
	for _, galEl := range galoisElement {
		rotKey.Keys[galEl] = NewSwitchingKey(ringDegree, moduliCount, decompSize)
	}
	return
}

// Set stores a copy of the rotKey SwitchingKey inside the receiver RotationKeySet
func (rtks *RotationKeySet) Set(galEl uint64, rotKey *SwitchingKey) {
	s := new(SwitchingKey)
	s.Copy(rotKey)
	rtks.Keys[galEl] = s
}

// Delete empties the set of rotation keys
func (rtks *RotationKeySet) Delete() {
	for k := range rtks.Keys {
		delete(rtks.Keys, k)
	}
}

// GetRotationKey return the rotation key for the given galois element or nil if such key is not in the set. The
// second argument is true  iff the first one is non-nil.
func (rtks *RotationKeySet) GetRotationKey(galoisEl uint64) (*SwitchingKey, bool) {
	rotKey, inSet := rtks.Keys[galoisEl]
	return rotKey, inSet
}

// Get returns the switching key backing slice.
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.Value
}

// Copy copies the other SwitchingKey inside the receiver.
func (swk *SwitchingKey) Copy(other *SwitchingKey) {
	if other == nil {
		return
	}
	if len(swk.Value) == 0 {
		swk.Value = make([][2]*ring.Poly, len(other.Value), len(other.Value))
		for i, o := range other.Value {
			n, q := uint64(o[0].GetDegree()), uint64(o[0].GetLenModuli())
			swk.Value[i] = [2]*ring.Poly{ring.NewPoly(n, q), ring.NewPoly(n, q)}
		}
	}
	for i, o := range other.Value {
		swk.Value[i][0].Copy(o[0])
		swk.Value[i][1].Copy(o[1])
	}
}

func NewSwitchingKey(ringDegree, moduliCount, decompSize uint64) *SwitchingKey {

	swk := new(SwitchingKey)

	swk.Value = make([][2]*ring.Poly, int(decompSize))

	for i := uint64(0); i < decompSize; i++ {
		swk.Value[i][0] = ring.NewPoly(ringDegree, moduliCount)
		swk.Value[i][1] = ring.NewPoly(ringDegree, moduliCount)
	}

	return swk
}

// NewRelinKey creates a new EvaluationKey with zero values.
func NewRelinKey(maxRelinDegree, ringDegree, moduliCount, decompSize uint64) (evakey *RelinearizationKey) {

	evakey = new(RelinearizationKey)

	evakey.Keys = make([]*SwitchingKey, maxRelinDegree)

	for d := uint64(0); d < maxRelinDegree; d++ {
		evakey.Keys[d] = NewSwitchingKey(ringDegree, moduliCount, decompSize)
	}

	return
}

// Get returns the slice of SwitchingKeys of the target EvaluationKey.
func (evk *RelinearizationKey) Get() []*SwitchingKey {
	return evk.Keys
}

// Set sets the polynomial of the target EvaluationKey as the input polynomials.
func (evk *RelinearizationKey) Set(rlk [][][2]*ring.Poly) {

	evk.Keys = make([]*SwitchingKey, len(rlk))
	for i := range rlk {
		evk.Keys[i] = new(SwitchingKey)
		evk.Keys[i].Value = make([][2]*ring.Poly, len(rlk[i]))
		for j := range rlk[i] {
			evk.Keys[i].Value[j][0] = rlk[i][j][0].CopyNew()
			evk.Keys[i].Value[j][1] = rlk[i][j][1].CopyNew()
		}
	}
}

// GetDataLen returns the length in bytes of the target SecretKey.
func (sk *SecretKey) GetDataLen(WithMetadata bool) (dataLen uint64) {
	return sk.Value.GetDataLen(WithMetadata)
}

// MarshalBinary encodes a secret key in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, sk.GetDataLen(true))

	if _, err = sk.Value.WriteTo(data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	sk.Value = new(ring.Poly)

	if _, err = sk.Value.DecodePolyNew(data); err != nil {
		return err
	}

	return nil
}

// GetDataLen returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	for _, el := range pk.Value {
		dataLen += el.GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes a PublicKey in a byte slice.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {

	dataLen := pk.GetDataLen(true)

	data = make([]byte, dataLen)

	var pointer, inc uint64

	if inc, err = pk.Value[0].WriteTo(data[pointer:]); err != nil {
		return nil, err
	}

	if _, err = pk.Value[1].WriteTo(data[pointer+inc:]); err != nil {
		return nil, err
	}

	return data, err

}

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {

	var pointer, inc uint64

	pk.Value[0] = new(ring.Poly)
	pk.Value[1] = new(ring.Poly)

	if inc, err = pk.Value[0].DecodePolyNew(data[pointer:]); err != nil {
		return err
	}

	if _, err = pk.Value[1].DecodePolyNew(data[pointer+inc:]); err != nil {
		return err
	}

	return nil
}

// GetDataLen returns the length in bytes of the target EvaluationKey.
func (evaluationkey *RelinearizationKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	if WithMetadata {
		dataLen++
	}

	for _, evakey := range evaluationkey.Keys {
		dataLen += (*SwitchingKey)(evakey).GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an EvaluationKey key in a byte slice.
func (evaluationkey *RelinearizationKey) MarshalBinary() (data []byte, err error) {

	var pointer uint64

	dataLen := evaluationkey.GetDataLen(true)

	data = make([]byte, dataLen)

	data[0] = uint8(len(evaluationkey.Keys))

	pointer++

	for _, evakey := range evaluationkey.Keys {

		if pointer, err = (*SwitchingKey)(evakey).encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled EvaluationKey in the target EvaluationKey.
func (evaluationkey *RelinearizationKey) UnmarshalBinary(data []byte) (err error) {

	deg := uint64(data[0])

	evaluationkey.Keys = make([]*SwitchingKey, deg)

	pointer := uint64(1)
	var inc uint64
	for i := uint64(0); i < deg; i++ {
		evaluationkey.Keys[i] = new(SwitchingKey)
		if inc, err = evaluationkey.Keys[i].decode(data[pointer:]); err != nil {
			return err
		}
		pointer += inc
	}

	return nil
}

// GetDataLen returns the length in bytes of the target SwitchingKey.
func (switchkey *SwitchingKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	if WithMetadata {
		dataLen++
	}

	for j := uint64(0); j < uint64(len(switchkey.Value)); j++ {
		dataLen += switchkey.Value[j][0].GetDataLen(WithMetadata)
		dataLen += switchkey.Value[j][1].GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an SwitchingKey in a byte slice.
func (switchkey *SwitchingKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, switchkey.GetDataLen(true))

	if _, err = switchkey.encode(0, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decode a previously marshaled SwitchingKey in the target SwitchingKey.
func (switchkey *SwitchingKey) UnmarshalBinary(data []byte) (err error) {

	if _, err = switchkey.decode(data); err != nil {
		return err
	}

	return nil
}

func (switchkey *SwitchingKey) encode(pointer uint64, data []byte) (uint64, error) {

	var err error

	var inc uint64

	data[pointer] = uint8(len(switchkey.Value))

	pointer++

	for j := uint64(0); j < uint64(len(switchkey.Value)); j++ {

		if inc, err = switchkey.Value[j][0].WriteTo(data[pointer : pointer+switchkey.Value[j][0].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = switchkey.Value[j][1].WriteTo(data[pointer : pointer+switchkey.Value[j][1].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc
	}

	return pointer, nil
}

func (switchkey *SwitchingKey) decode(data []byte) (pointer uint64, err error) {

	decomposition := uint64(data[0])

	pointer = uint64(1)

	switchkey.Value = make([][2]*ring.Poly, decomposition)

	var inc uint64

	for j := uint64(0); j < decomposition; j++ {

		switchkey.Value[j][0] = new(ring.Poly)
		if inc, err = switchkey.Value[j][0].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

		switchkey.Value[j][1] = new(ring.Poly)
		if inc, err = switchkey.Value[j][1].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

	}

	return pointer, nil
}

// GetDataLen returns the length in bytes of the target RotationKeys.
func (rotationkey *RotationKeySet) GetDataLen(WithMetaData bool) (dataLen uint64) {
	for _, k := range rotationkey.Keys {
		if WithMetaData {
			dataLen += 4
		}
		dataLen += k.GetDataLen(WithMetaData)
	}
	return
}

// MarshalBinary encodes a RotationKeys struct in a byte slice.
func (rotationkey *RotationKeySet) MarshalBinary() (data []byte, err error) {

	data = make([]byte, rotationkey.GetDataLen(true))

	pointer := uint64(0)

	for galEL, key := range rotationkey.Keys {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(galEL))
		pointer += 4

		if pointer, err = key.encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled RotationKeys in the target RotationKeys.
func (rotationkey *RotationKeySet) UnmarshalBinary(data []byte) (err error) {

	if rotationkey.Keys == nil {
		rotationkey.Keys = make(map[uint64]*SwitchingKey)
	} else {
		rotationkey.Delete()
	}

	for len(data) > 0 {

		galEl := uint64(binary.BigEndian.Uint32(data))
		data = data[4:]

		swk := new(SwitchingKey)
		var inc uint64
		if inc, err = swk.decode(data); err != nil {
			return err
		}
		data = data[inc:]
		rotationkey.Keys[galEl] = swk

	}

	return nil
}
