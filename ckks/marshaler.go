package ckks

import (
	"encoding/binary"
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math/bits"
)

// MarshalBinary encodes a ciphertext on a byte slice. The total size
// in byte is 4 + 2 * 8 * N * (level + 1).
func (ciphertext *Ciphertext) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(ciphertext.Value()[0].Coeffs[0]))
	level := uint64(len(ciphertext.Value()[0].Coeffs))
	degree := ciphertext.Degree()
	scale := ciphertext.Scale()

	if level > 0xFF {
		return nil, errors.New("error : ciphertext level overflow uint8")
	}

	if degree > 0xFF {
		return nil, errors.New("error : ciphertext degree overflow uint8")
	}

	if scale > 0xFF {
		return nil, errors.New("error : ciphertext scale overflow uint8")
	}

	data := make([]byte, 4+((N*level*(degree+1))<<3))

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(degree)
	data[3] = uint8(scale)

	pointer := uint64(4)

	for i := uint64(0); i < degree+1; i++ {
		if pointer, err = ring.WriteCoeffsTo(pointer, N, level, ciphertext.Value()[i].Coeffs, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled ciphertext on the target ciphertext.
// The target ciphertext must be of the appropriate format and size, it can be created with the
// methode NewCiphertext(uint64, uin64t, uint64).
func (ciphertext *Ciphertext) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	degree := uint64(data[2])

	pointer := uint64(4)

	if ciphertext.Degree() != degree {
		return errors.New("error : invalid ciphertext encoding (unexpected degree)")
	}

	if ciphertext.Level() != level-1 {
		return errors.New("error : invalid ciphertext encoding (unexpected level)")
	}

	if ((uint64(len(data)) - pointer) >> 3) != N*level*(degree+1) {
		return errors.New("error : invalid ciphertext encoding (unexpected data length)")
	}

	ciphertext.SetScale(uint64(data[3]))

	for x := uint64(0); x < degree+1; x++ {
		pointer, _ = ring.DecodeCoeffs(pointer, N, level, ciphertext.Value()[x].Coeffs, data)
	}

	return nil
}

// MarshalBinary encodes a secret-key on a byte slice. The total size in byte is 1 + N/4.
func (sk *SecretKey) MarshalBinary(ckkscontext *CkksContext) ([]byte, error) {

	var Q, x uint64

	tmp := ckkscontext.contextLevel[0].NewPoly()

	ckkscontext.contextLevel[0].InvNTT(sk.sk, tmp)
	ckkscontext.contextLevel[0].InvMForm(tmp, tmp)

	N := uint64(len(sk.sk.Coeffs[0]))

	data := make([]byte, 1+(N>>2))

	Q = ckkscontext.contextLevel[0].Modulus[0]

	data[0] = uint8(bits.Len64(N) - 1)

	for i, coeff := range tmp.Coeffs[0] {

		x = ((coeff + 1) - Q)

		data[1+(i>>2)] <<= 2
		// encode 0 = 0b00, 1 = 0b01, -1 = 0b10
		data[1+(i>>2)] |= uint8(((x>>63)^1)<<1 | x&1)
	}

	return data, nil
}

// UnMarshalBinary decode a previously marshaled secret-key on the target secret-key.
// The target secret-key must be of the appropriate format, it can be created with the methode NewSecretKeyEmpty().
func (sk *SecretKey) UnMarshalBinary(data []byte, ckkscontext *CkksContext) error {

	var N, coeff uint64

	N = uint64(1 << data[0])

	if uint64(len(sk.sk.Coeffs[0])) != N {
		return errors.New("error : invalid secret key encoding (logN do not match)")
	}

	if uint64(len(data)) != 1+(N>>2) {
		return errors.New("error : invalid secret key encoding (unexpected data length)")
	}

	for i := uint64(0); i < N; i++ {

		coeff = uint64(data[1+(i>>2)]>>(6-((i<<1)&7))) & 3

		for j := range ckkscontext.moduli {
			sk.sk.Coeffs[j][i] = (ckkscontext.moduli[j]-1)*(coeff>>1) | (coeff & 1)
		}
	}

	ckkscontext.keyscontext.NTT(sk.sk, sk.sk)
	ckkscontext.keyscontext.MForm(sk.sk, sk.sk)

	return nil
}

// MarshalBinary encodes a public-key on a byte slice. The total size is 2 + 16 * N * (level + 1).
func (pk *PublicKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(pk.pk[0].Coeffs[0]))
	levels := uint64(len(pk.pk[0].Coeffs))

	if levels > 0xFF {
		return nil, errors.New("error : max degree uint8 overflow")
	}

	data := make([]byte, 2+((N*levels)<<4))

	data[0] = uint8((bits.Len64(uint64(N)) - 1))
	data[1] = uint8(levels)

	pointer := uint64(2)

	if pointer, err = ring.WriteCoeffsTo(pointer, N, levels, pk.pk[0].Coeffs, data); err != nil {
		return nil, err
	}

	if pointer, err = ring.WriteCoeffsTo(pointer, N, levels, pk.pk[1].Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled public-key on the target public-key.
// The target public-key must have the appropriate format and size, it can be created with
// the methode NewPublicKeyEmpty().
func (pk *PublicKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	levels := uint64(data[1])

	pointer := uint64(2)

	if uint64(len(pk.pk[0].Coeffs[0])) != N {
		return errors.New("error : invalid publickey[0] receiver (logN do not match)")
	}

	if uint64(len(pk.pk[0].Coeffs[1])) != N {
		return errors.New("error : invalid publickey[1] receiver (logN do not match)")
	}

	if uint64(len(pk.pk[0].Coeffs)) != levels {
		return errors.New("error : invalid publickey[0] receiver (level do not match data)")
	}

	if uint64(len(pk.pk[1].Coeffs)) != levels {
		return errors.New("error : invalid publickey[1] receiver (level do not match data)")
	}

	if ((uint64(len(data)) - pointer) >> 4) != (N * levels) {
		return errors.New("error : invalid PublicKey encoding")
	}

	pointer, _ = ring.DecodeCoeffs(pointer, N, levels, pk.pk[0].Coeffs, data)
	pointer, _ = ring.DecodeCoeffs(pointer, N, levels, pk.pk[1].Coeffs, data)

	return nil
}

// MarshalBinary encodes an evaluation key on a byte slice. The total size depends on each modulus size and the bit decomp.
// It will approximately be 5 + (level + 1) * ( 1 + 2 * 8 * N * (level + 1) * logQi/bitDecomp).
func (evaluationkey *EvaluationKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(evaluationkey.evakey.evakey[0][0][0].Coeffs[0]))
	levels := uint64(len(evaluationkey.evakey.evakey[0][0][0].Coeffs))
	decomposition := levels
	bitDecomp := evaluationkey.evakey.bitDecomp

	if levels > 0xFF {
		return nil, errors.New("error : max number modulis uint8 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint8 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint8 overflow")
	}

	var dataLen uint64
	dataLen = 4

	for j := uint64(0); j < decomposition; j++ {
		dataLen += 1                                                                                //Information about the size of the bitdecomposition
		dataLen += 2 * 8 * N * levels * decomposition * uint64(len(evaluationkey.evakey.evakey[j])) // nb coefficients * 8
	}

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(levels)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)

	pointer := uint64(4)

	var bitLog uint8

	for j := uint64(0); j < decomposition; j++ {
		bitLog = uint8(len(evaluationkey.evakey.evakey[j]))
		data[pointer] = bitLog
		pointer += 1
		for x := uint8(0); x < bitLog; x++ {
			if pointer, err = ring.WriteCoeffsTo(pointer, N, levels, evaluationkey.evakey.evakey[j][x][0].Coeffs, data); err != nil {
				return nil, err
			}

			if pointer, err = ring.WriteCoeffsTo(pointer, N, levels, evaluationkey.evakey.evakey[j][x][1].Coeffs, data); err != nil {
				return nil, err
			}
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled evaluation-key on the target evaluation-key. The target evaluation-key
// must have the appropriate format and size, it can be created with the methode NewRelinKeyEmpty(uint64, uint64).
func (evaluationkey *EvaluationKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	levels := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])

	pointer := uint64(4)
	var bitLog uint64

	evaluationkey.evakey.bitDecomp = bitDecomp

	for j := uint64(0); j < decomposition; j++ {

		bitLog = uint64(data[pointer])
		pointer += 1

		for x := uint64(0); x < bitLog; x++ {

			if uint64(len(evaluationkey.evakey.evakey[j][x][0].Coeffs)) != levels {
				return errors.New("error : evaluationkey receiver (level do not match data)")
			}

			if uint64(len(evaluationkey.evakey.evakey[j][x][1].Coeffs)) != levels {
				return errors.New("error : evaluationkey receiver (level do not match data)")
			}

			pointer, _ = ring.DecodeCoeffs(pointer, N, levels, evaluationkey.evakey.evakey[j][x][0].Coeffs, data)
			pointer, _ = ring.DecodeCoeffs(pointer, N, levels, evaluationkey.evakey.evakey[j][x][1].Coeffs, data)
		}
	}

	return nil
}

// MarshalBinary encodes an switching-key on a byte slice. The total size in byte will be approximately 5 + (level + 1) * ( 1 + 2 * 8 * N * (level + 1) * logQi/bitDecomp).
func (switchingkey *SwitchingKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(switchingkey.evakey[0][0][0].Coeffs[0]))
	level := uint64(len(switchingkey.evakey[0][0][0].Coeffs))
	decomposition := level
	bitDecomp := switchingkey.bitDecomp

	if level > 0xFF {
		return nil, errors.New("error : max number modulis uint8 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint8 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint8 overflow")
	}

	var dataLen uint64
	dataLen = 4

	for j := uint64(0); j < decomposition; j++ {
		dataLen += 1                                                                       //Information about the size of the bitdecomposition
		dataLen += 2 * 8 * N * level * decomposition * uint64(len(switchingkey.evakey[j])) // nb coefficients * 8
	}

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)

	pointer := uint64(4)

	var bitLog uint8

	for j := uint64(0); j < decomposition; j++ {
		bitLog = uint8(len(switchingkey.evakey[j]))
		data[pointer] = bitLog
		pointer += 1
		for x := uint8(0); x < bitLog; x++ {
			if pointer, err = ring.WriteCoeffsTo(pointer, N, level, switchingkey.evakey[j][x][0].Coeffs, data); err != nil {
				return nil, err
			}

			if pointer, err = ring.WriteCoeffsTo(pointer, N, level, switchingkey.evakey[j][x][1].Coeffs, data); err != nil {
				return nil, err
			}
		}
	}

	return data, nil
}

// UnMarshalBinary decode a previously marshaled switching-key on the target switching-key.
// The target switching-key must have the appropriate format and size, it can be created with the methode NewSwitchingKeyEmpty(uint64).
func (switchingkey *SwitchingKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])

	pointer := uint64(4)
	var bitLog uint64

	switchingkey.bitDecomp = bitDecomp

	for j := uint64(0); j < decomposition; j++ {

		bitLog = uint64(data[pointer])
		pointer += 1

		for x := uint64(0); x < bitLog; x++ {
			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchingkey.evakey[j][x][0].Coeffs, data)
			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchingkey.evakey[j][x][1].Coeffs, data)
		}
	}

	return nil
}

// MarshalBinary encodes a rotationkeys structure on a byte slice. The total size in byte is approximately
// 5 + 4*(nb left rot + num right rot) + (nb left rot + num right rot + 1 (if rotate row)) * (level + 1) * ( 1 + 2 * 8 * N * (level + 1) * logQi/bitDecomp).
func (rotationkey *RotationKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(rotationkey.ckkscontext.n)
	level := uint64(len(rotationkey.ckkscontext.keyscontext.Modulus))
	decomposition := level
	bitDecomp := rotationkey.bitDecomp
	mappingRow := 0
	mappingColL := []uint64{}
	mappingColR := []uint64{}

	if level > 0xFF {
		return nil, errors.New("error : max number modulis uint8 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint8 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint8 overflow")
	}

	var dataLen uint64
	dataLen = 13

	for i := uint64(1); i < N>>1; i++ {
		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColL = append(mappingColL, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                          //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}

		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColR = append(mappingColR, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                          //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}
	}

	if rotationkey.evakey_rot_row != nil {
		mappingRow = 1
		for j := uint64(0); j < decomposition; j++ {
			dataLen += 1                                                                                     //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakey_rot_row.evakey[j])) // nb coefficients * 8
		}
	}

	dataLen += uint64(len(mappingColL)+len(mappingColR)) << 2 // size needed to encode what rotation are present

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)
	data[4] = uint8(mappingRow)

	pointer := uint64(5)

	binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(len(mappingColL)))
	pointer += 4

	binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(len(mappingColR)))
	pointer += 4

	for _, i := range mappingColL {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(i))

		pointer += 4
	}

	for _, i := range mappingColR {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(i))

		pointer += 4
	}

	// Encodes the different rotation key indexes
	var bitLog uint8
	if mappingRow == 1 {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(rotationkey.evakey_rot_row.evakey[j]))
			data[pointer] = bitLog
			pointer += 1
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	for _, i := range mappingColL {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(rotationkey.evakey_rot_col_L[i].evakey[j]))
			data[pointer] = bitLog
			pointer += 1
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	for _, i := range mappingColR {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(rotationkey.evakey_rot_col_R[i].evakey[j]))
			data[pointer] = bitLog
			pointer += 1
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled rotation-keys on the target rotation-keys. In contrary to all
// the other structures, the unmarshaling for rotationkeys only need an empty receiver, as it is not possible to
// create receiver of the correct format and size without knowing all the content of the marshaled rotationkeys. The memory
// will be allocated on the fly.
func (rotationkey *RotationKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])
	mappingRow := uint64(data[4])
	mappingColL := make([]uint64, binary.BigEndian.Uint32(data[5:9]))
	mappingColR := make([]uint64, binary.BigEndian.Uint32(data[9:13]))

	rotationkey.bitDecomp = uint64(bitDecomp)

	rotationkey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
	//rotationkey.evakey_rot_col_R = make(map[uint64][][][2]*ring.Poly)

	pointer := uint64(13)

	for i := 0; i < len(mappingColL); i++ {
		mappingColL[i] = uint64(binary.BigEndian.Uint32(data[pointer : pointer+4]))
		pointer += 4
	}

	for i := 0; i < len(mappingColR); i++ {
		mappingColR[i] = uint64(binary.BigEndian.Uint32(data[pointer : pointer+4]))
		pointer += 4
	}

	var bitLog uint64
	if mappingRow == 1 {

		rotationkey.evakey_rot_row = new(SwitchingKey)
		rotationkey.evakey_rot_row.bitDecomp = bitDecomp
		rotationkey.evakey_rot_row.evakey = make([][][2]*ring.Poly, decomposition)

		for j := uint64(0); j < decomposition; j++ {

			bitLog = uint64(data[pointer])
			pointer += 1

			rotationkey.evakey_rot_row.evakey[j] = make([][2]*ring.Poly, bitLog)

			for x := uint64(0); x < bitLog; x++ {

				rotationkey.evakey_rot_row.evakey[j][x][0] = new(ring.Poly)
				rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs = make([][]uint64, level)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data)

				rotationkey.evakey_rot_row.evakey[j][x][1] = new(ring.Poly)
				rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs = make([][]uint64, level)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data)
			}
		}
	}

	if len(mappingColL) > 0 {

		rotationkey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)

		for _, i := range mappingColL {

			rotationkey.evakey_rot_col_L[i] = new(SwitchingKey)
			rotationkey.evakey_rot_col_L[i].bitDecomp = bitDecomp
			rotationkey.evakey_rot_col_L[i].evakey = make([][][2]*ring.Poly, decomposition)

			for j := uint64(0); j < decomposition; j++ {

				bitLog = uint64(data[pointer])
				pointer += 1

				rotationkey.evakey_rot_col_L[i].evakey[j] = make([][2]*ring.Poly, bitLog)

				for x := uint64(0); x < bitLog; x++ {

					rotationkey.evakey_rot_col_L[i].evakey[j][x][0] = new(ring.Poly)
					rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_L[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	if len(mappingColR) > 0 {

		rotationkey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)

		for _, i := range mappingColR {

			rotationkey.evakey_rot_col_R[i] = new(SwitchingKey)
			rotationkey.evakey_rot_col_R[i].bitDecomp = bitDecomp
			rotationkey.evakey_rot_col_R[i].evakey = make([][][2]*ring.Poly, decomposition)

			for j := uint64(0); j < decomposition; j++ {

				bitLog = uint64(data[pointer])
				pointer += 1

				rotationkey.evakey_rot_col_R[i].evakey[j] = make([][2]*ring.Poly, bitLog)

				for x := uint64(0); x < bitLog; x++ {

					rotationkey.evakey_rot_col_R[i].evakey[j][x][0] = new(ring.Poly)
					rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_R[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	return nil
}
