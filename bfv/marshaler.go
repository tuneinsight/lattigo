package bfv

import (
	"encoding/binary"
	"errors"
	"github.com/lca1/lattigo/ring"
	"math"
	"math/bits"
)

// MarshalBinary marshals a bfvcontext into bytes. It only encodes the minimum parameters
// to be able to reconstruct it. The total size in byte is 3 + 8 * (2 + numberModuliQ + numberModuliP).
func (bfvContext *BfvContext) MarshalBinary() ([]byte, error) {

	N := bfvContext.n
	numberModuliQ := len(bfvContext.contextQ.Modulus)
	numberModuliP := len(bfvContext.contextP.Modulus)

	if numberModuliQ > 0xFF {
		return nil, errors.New("cannot marshal bfvcontext -> numberModuliQ overflow uint8")
	}

	if numberModuliP > 0xFF {
		return nil, errors.New("cannot marshal bfvcontext -> numberModuliP overflow uint8")
	}

	data := make([]byte, 3+((2+numberModuliQ+numberModuliP)<<3))

	data[0] = uint8(bits.Len64(N) - 1)
	data[1] = uint8(numberModuliQ)
	data[2] = uint8(numberModuliP)

	pointer := 3

	binary.BigEndian.PutUint64(data[pointer:pointer+8], uint64(bfvContext.sigma*(1<<32)))
	pointer += 8

	binary.BigEndian.PutUint64(data[pointer:pointer+8], bfvContext.contextT.Modulus[0])
	pointer += 8

	for i := 0; i < numberModuliQ; i++ {
		binary.BigEndian.PutUint64(data[pointer:pointer+8], bfvContext.contextQ.Modulus[i])
		pointer += 8
	}

	for i := 0; i < numberModuliP; i++ {
		binary.BigEndian.PutUint64(data[pointer:pointer+8], bfvContext.contextP.Modulus[i])
		pointer += 8
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled bfvcontext on the target bfvcontext.
// Since only the minimum amount of information is encoded in the bytes, it will need to 
// re-compute all the internal variables and parameters.
func (bfvContext *BfvContext) UnMarshalBinary(data []byte) error {
	//var err error

	N := uint64(1 << data[0])
	numberModuliQ := int(data[1])
	numberModuliP := int(data[2])

	pointer := 3

	if ((len(data) - pointer) >> 3) != (2 + numberModuliQ + numberModuliP) {
		return errors.New("cannot unmarshal bfvcontext -> invlid encoding")
	}

	sigma := math.Round((float64(binary.BigEndian.Uint64(data[pointer:pointer+8]))/float64(1<<32))*100) / 100
	pointer += 8

	t := binary.BigEndian.Uint64(data[pointer : pointer+8])
	pointer += 8

	ModuliesQ := make([]uint64, numberModuliQ)
	for i := 0; i < numberModuliQ; i++ {
		ModuliesQ[i] = binary.BigEndian.Uint64(data[pointer : pointer+8])
		pointer += 8
	}

	ModuliesP := make([]uint64, numberModuliP)
	for i := 0; i < numberModuliP; i++ {
		ModuliesP[i] = binary.BigEndian.Uint64(data[pointer : pointer+8])
		pointer += 8
	}

	bfvContext.SetParameters(N, t, ModuliesQ, ModuliesP, sigma)

	return nil
}

// MarshalBinary encodes a ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() ([]byte, error) {

	if ciphertext.IsNTT() {
		return nil, errors.New("cannot marshal ciphertext -> ciphertext is in the NTT domain")
	}

	var err error

	N := uint64(len(ciphertext.Value()[0].Coeffs[0]))
	level := uint64(len(ciphertext.Value()[0].Coeffs))
	degree := ciphertext.Degree()

	if level > 0xFF {
		return nil, errors.New("cannot marshal ciphertext -> ciphertext numberModuli overflow uint8")
	}

	if degree > 0xFF {
		return nil, errors.New("cannot marshal ciphertext -> ciphertext degree overflow uint8")
	}

	data := make([]byte, 4+((N*level*(degree+1))<<3))

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(degree)

	pointer := uint64(3)

	for i := uint64(0); i < degree+1; i++ {
		if pointer, err = ring.WriteCoeffsTo(pointer, N, level, ciphertext.Value()[i].Coeffs, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled ciphertext on the target ciphertext.
// The target ciphertext must be of the appropriate format and size, it can be created with the
// methode NewCiphertext(uint64).
func (ciphertext *Ciphertext) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	degree := uint64(data[2])

	pointer := uint64(3)

	if ciphertext.Degree() != degree {
		return errors.New("cannot unmarshal ciphertext -> invalid ciphertext encoding (unexpected degree)")
	}

	if uint64(len(ciphertext.Value()[0].Coeffs)) != level {
		return errors.New("cannot unmarshal ciphertext -> invalid ciphertext encoding (unexpected number of moduli)")
	}

	if ((uint64(len(data)) - pointer) >> 3) != N*level*(degree+1) {
		return errors.New("cannot unmarshal ciphertext -> invalid ciphertext encoding (unexpected data length)")
	}

	for x := uint64(0); x < degree+1; x++ {
		pointer, _ = ring.DecodeCoeffs(pointer, N, level, ciphertext.Value()[x].Coeffs, data)
	}

	return nil
}

// MarshalBinary encodes a secret-key on a byte slice. The total size in byte is 1 + N/4.
func (sk *SecretKey) MarshalBinary(bfvcontext *BfvContext) ([]byte, error) {

	var Q, x uint64

	tmp := bfvcontext.contextQ.NewPoly()

	bfvcontext.contextQ.InvNTT(sk.sk, tmp)
	bfvcontext.contextQ.InvMForm(tmp, tmp)

	N := uint64(len(sk.sk.Coeffs[0]))

	data := make([]byte, 1+(N>>2))

	Q = bfvcontext.contextQ.Modulus[0]

	data[0] = uint8(bits.Len64(N) - 1)

	for i, coeff := range tmp.Coeffs[0] {

		x = ((coeff + 1) - Q)

		data[1+(i>>2)] <<= 2
		// encodes 0 = 0b00, 1 = 0b01, -1 = 0b10
		data[1+(i>>2)] |= uint8(((x>>63)^1)<<1 | x&1)
	}

	return data, nil
}

// UnMarshalBinary decode a previously marshaled secret-key on the target secret-key.
// The target secret-key must be of the appropriate format, it can be created with the methode NewSecretKeyEmpty().
func (sk *SecretKey) UnMarshalBinary(data []byte, bfvcontext *BfvContext) error {

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

		for j := range bfvcontext.contextQ.Modulus {
			sk.sk.Coeffs[j][i] = (bfvcontext.contextQ.Modulus[j]-1)*(coeff>>1) | (coeff & 1)
		}
	}

	bfvcontext.contextQ.NTT(sk.sk, sk.sk)
	bfvcontext.contextQ.MForm(sk.sk, sk.sk)

	return nil
}

// MarshalBinary encodes a public-key on a byte slice. The total size is 2 + 16 * N * numberModuliQ.
func (pk *PublicKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(pk.pk[0].Coeffs[0]))
	numberModuli := uint64(len(pk.pk[0].Coeffs))

	if numberModuli > 0xFF {
		return nil, errors.New("error : max degree uint16 overflow")
	}

	data := make([]byte, 2+((N*numberModuli)<<4))

	data[0] = uint8((bits.Len64(uint64(N)) - 1))
	data[1] = uint8(numberModuli)

	pointer := uint64(2)

	if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, pk.pk[0].Coeffs, data); err != nil {
		return nil, err
	}

	if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, pk.pk[1].Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled public-key on the target public-key.
// The target public-key must have the appropriate format and size, it can be created with
// the methode NewPublicKeyEmpty().
func (pk *PublicKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModuli := uint64(data[1])

	pointer := uint64(2)

	if ((uint64(len(data)) - pointer) >> 4) != (N * numberModuli) {
		return errors.New("error : invalid publickey encoding")
	}

	pointer, _ = ring.DecodeCoeffs(pointer, N, numberModuli, pk.pk[0].Coeffs, data)
	pointer, _ = ring.DecodeCoeffs(pointer, N, numberModuli, pk.pk[1].Coeffs, data)

	return nil
}

// MarshalBinary encodes an evaluation key on a byte slice. The total size depends on each modulus size and the bit decomp.
// It will approximately be 5 + maxDegree * numberModuli * ( 1 + 2 * 8 * N * numberModuli * logQi/bitDecomp).
func (evaluationkey *EvaluationKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(evaluationkey.evakey[0].evakey[0][0][0].Coeffs[0]))
	numberModuli := uint64(len(evaluationkey.evakey[0].evakey[0][0][0].Coeffs))
	decomposition := numberModuli
	bitDecomp := evaluationkey.evakey[0].bitDecomp

	maxDegree := uint64(len(evaluationkey.evakey))

	if numberModuli > 0xFF {
		return nil, errors.New("cannot marshal evaluationkey -> max number moduli uint16 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("cannot marshal evaluationkey -> max decomposition uint16 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("cannot marshal evaluationkey -> max bitDecomp uint16 overflow")
	}

	if maxDegree > 0xFF {
		return nil, errors.New("cannot marshal evaluationkey -> max degree uint16 overflow")
	}

	var dataLen uint64
	dataLen = 5
	for i := uint64(0); i < maxDegree; i++ {
		for j := uint64(0); j < decomposition; j++ {
			dataLen += 1                                                                                                       //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * numberModuli * uint64(len(evaluationkey.evakey[i].evakey[j])) // nb coefficients * 8
		}
	}

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModuli)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)
	data[4] = uint8(maxDegree)

	pointer := uint64(5)

	var bitLog uint8
	for i := uint64(0); i < maxDegree; i++ {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(evaluationkey.evakey[i].evakey[j]))
			data[pointer] = bitLog
			pointer += 1
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, evaluationkey.evakey[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, evaluationkey.evakey[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled evaluation-key on the target evaluation-key. The target evaluation-key 
// must have the appropriate format and size, it can be created with the methode NewRelinKeyEmpty(uint64, uint64).
func (evaluationkey *EvaluationKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModuli := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])
	maxDegree := uint64(data[4])

	pointer := uint64(5)
	var bitLog uint64
	for i := uint64(0); i < maxDegree; i++ {

		evaluationkey.evakey[i].bitDecomp = bitDecomp

		for j := uint64(0); j < decomposition; j++ {

			bitLog = uint64(data[pointer])
			pointer += 1

			for x := uint64(0); x < bitLog; x++ {

				if uint64(len(evaluationkey.evakey[i].evakey[j][x][0].Coeffs)) != numberModuli {
					return errors.New("cannot unmarshal evaluation-key -> receiver (numberModuli does not match data)")
				}

				if uint64(len(evaluationkey.evakey[i].evakey[j][x][1].Coeffs)) != numberModuli {
					return errors.New("cannot unmarshal evaluation-key -> receiver (numberModuli does not match data)")
				}

				pointer, _ = ring.DecodeCoeffs(pointer, N, numberModuli, evaluationkey.evakey[i].evakey[j][x][0].Coeffs, data)
				pointer, _ = ring.DecodeCoeffs(pointer, N, numberModuli, evaluationkey.evakey[i].evakey[j][x][1].Coeffs, data)
			}
		}
	}

	return nil
}

// MarshalBinary encodes an switching-key on a byte slice. The total size in byte will be approximately 5 + numberModuli * ( 1 + 2 * 8 * N * numberModuli * logQi/bitDecomp).
func (switchkey *SwitchingKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(switchkey.evakey[0][0][0].Coeffs[0]))
	level := uint64(len(switchkey.evakey[0][0][0].Coeffs))
	decomposition := level
	bitDecomp := switchkey.bitDecomp

	if level > 0xFF {
		return nil, errors.New("cannot marshal switching-key -> max number modulie uint8 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("cannot marshal switching-key -> max decomposition uint8 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("cannot marshal switching-key -> max bitDecomp uint8 overflow")
	}

	var dataLen uint64
	dataLen = 4

	for j := uint64(0); j < decomposition; j++ {
		dataLen += 1                                                                    //Information about the size of the bitdecomposition
		dataLen += 2 * 8 * N * level * decomposition * uint64(len(switchkey.evakey[j])) // nb coefficients * 8
	}

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)

	pointer := uint64(4)

	var bitLog uint8

	for j := uint64(0); j < decomposition; j++ {
		bitLog = uint8(len(switchkey.evakey[j]))
		data[pointer] = bitLog
		pointer += 1
		for x := uint8(0); x < bitLog; x++ {
			if pointer, err = ring.WriteCoeffsTo(pointer, N, level, switchkey.evakey[j][x][0].Coeffs, data); err != nil {
				return nil, err
			}

			if pointer, err = ring.WriteCoeffsTo(pointer, N, level, switchkey.evakey[j][x][1].Coeffs, data); err != nil {
				return nil, err
			}
		}
	}

	return data, nil
}

// UnMarshalBinary decode a previously marshaled switching-key on the target switching-key.
// The target switching-key must have the appropriate format and size, it can be created with the methode NewSwitchingKeyEmpty(uint64).
func (switchkey *SwitchingKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])

	pointer := uint64(4)
	var bitLog uint64

	switchkey.bitDecomp = bitDecomp

	for j := uint64(0); j < decomposition; j++ {

		bitLog = uint64(data[pointer])
		pointer += 1

		for x := uint64(0); x < bitLog; x++ {

			if uint64(len(switchkey.evakey[j][x][0].Coeffs)) != decomposition {
				return errors.New("cannot unmarshal switching-key -> receiver (numberModuli does not match data)")
			}

			if uint64(len(switchkey.evakey[j][x][1].Coeffs)) != decomposition {
				return errors.New("cannot unmarshal switching-key -> receiver (numberModuli does not match data)")
			}


			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchkey.evakey[j][x][0].Coeffs, data)
			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchkey.evakey[j][x][1].Coeffs, data)
		}
	}

	return nil
}

// MarshalBinary encodes a rotationkeys structure on a byte slice. The total size in byte is approximately
// 5 + 4*(nb left rot + num right rot) + (nb left rot + num right rot + 1 (if rotate row)) * maxDegree * numberModuli * ( 1 + 2 * 8 * N * numberModuli * logQi/bitDecomp).
func (rotationkey *RotationKeys) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(rotationkey.bfvcontext.n)
	numberModuli := uint64(len(rotationkey.bfvcontext.contextQ.Modulus))
	decomposition := numberModuli
	bitDecomp := rotationkey.bitDecomp
	mappingRow := 0
	mappingColL := []uint64{}
	mappingColR := []uint64{}

	if numberModuli > 0xFF {
		return nil, errors.New("error : max number modulies uint16 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint16 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint16 overflow")
	}

	var dataLen uint64
	dataLen = 13

	for i := uint64(1); i < rotationkey.bfvcontext.n>>1; i++ {
		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColL = append(mappingColL, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                                   //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * numberModuli * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}

		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColR = append(mappingColR, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                                   //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * numberModuli * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}
	}

	if rotationkey.evakey_rot_row != nil {
		mappingRow = 1
		for j := uint64(0); j < decomposition; j++ {
			dataLen += 1                                                                                              //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * numberModuli * decomposition * uint64(len(rotationkey.evakey_rot_row.evakey[j])) // nb coefficients * 8
		}
	}

	dataLen += uint64(len(mappingColL)+len(mappingColR)) << 2 // size needed to encode what rotation are present

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModuli)
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data); err != nil {
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data); err != nil {
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data); err != nil {
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
func (rotationkey *RotationKeys) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModuli := uint64(data[1])
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
				rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs = make([][]uint64, numberModuli)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data)

				rotationkey.evakey_rot_row.evakey[j][x][1] = new(ring.Poly)
				rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs = make([][]uint64, numberModuli)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data)
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
					rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs = make([][]uint64, numberModuli)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_L[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs = make([][]uint64, numberModuli)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data)
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
					rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs = make([][]uint64, numberModuli)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_R[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs = make([][]uint64, numberModuli)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	return nil
}
