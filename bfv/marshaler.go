package bfv

import (
	"encoding/binary"
	"errors"
	"github.com/lca1/lattigo/ring"
	"math"
	"math/bits"
)

func (bfvContext *BfvContext) MarshalBinary() ([]byte, error) {

	N := bfvContext.n
	numberModuliesQ := len(bfvContext.contextQ.Modulus)
	numberModuliesP := len(bfvContext.contextP.Modulus)

	if numberModuliesQ > 0xFF {
		return nil, errors.New("error : bfvcontext numberModuliesQ overflow uint8")
	}

	if numberModuliesP > 0xFF {
		return nil, errors.New("error : bfvcontext numberModuliesP overflow uint8")
	}

	data := make([]byte, 3+((2+numberModuliesQ+numberModuliesP)<<3))

	data[0] = uint8(bits.Len64(N) - 1)
	data[1] = uint8(numberModuliesQ)
	data[2] = uint8(numberModuliesP)

	pointer := 3

	binary.BigEndian.PutUint64(data[pointer:pointer+8], uint64(bfvContext.sigma*(1<<32)))
	pointer += 8

	binary.BigEndian.PutUint64(data[pointer:pointer+8], bfvContext.contextT.Modulus[0])
	pointer += 8

	for i := 0; i < numberModuliesQ; i++ {
		binary.BigEndian.PutUint64(data[pointer:pointer+8], bfvContext.contextQ.Modulus[i])
		pointer += 8
	}

	for i := 0; i < numberModuliesP; i++ {
		binary.BigEndian.PutUint64(data[pointer:pointer+8], bfvContext.contextP.Modulus[i])
		pointer += 8
	}

	return data, nil
}

func (bfvContext *BfvContext) UnMarshalBinary(data []byte) error {
	//var err error

	N := uint64(1 << data[0])
	numberModuliesQ := int(data[1])
	numberModuliesP := int(data[2])

	pointer := 3

	if ((len(data) - pointer) >> 3) != (2 + numberModuliesQ + numberModuliesP) {
		return errors.New("error : invalid PublicKey encoding")
	}

	sigma := math.Round((float64(binary.BigEndian.Uint64(data[pointer:pointer+8]))/float64(1<<32))*100) / 100
	pointer += 8

	t := binary.BigEndian.Uint64(data[pointer : pointer+8])
	pointer += 8

	ModuliesQ := make([]uint64, numberModuliesQ)
	for i := 0; i < numberModuliesQ; i++ {
		ModuliesQ[i] = binary.BigEndian.Uint64(data[pointer : pointer+8])
		pointer += 8
	}

	ModuliesP := make([]uint64, numberModuliesP)
	for i := 0; i < numberModuliesP; i++ {
		ModuliesP[i] = binary.BigEndian.Uint64(data[pointer : pointer+8])
		pointer += 8
	}

	bfvContext.SetParameters(N, t, ModuliesQ, ModuliesP, sigma)

	return nil
}

func (P *Plaintext) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(P.value[0].Coeffs[0]))

	data := make([]byte, 1+(N<<3))

	data[0] = uint8(bits.Len64(uint64(N)) - 1)

	pointer := uint64(1)

	if _, err = ring.WriteCoeffsTo(pointer, N, 1, P.value[0].Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

func (P *Plaintext) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])

	P.value = make([]*ring.Poly, 1)
	P.value[0] = new(ring.Poly)
	P.value[0].Coeffs = make([][]uint64, 1)
	P.value[0].Coeffs[0] = make([]uint64, N)

	pointer := uint64(1)

	if ((uint64(len(data)) - pointer) >> 3) != N {
		return errors.New("error : invalid Plaintext encoding")
	}

	pointer, _ = ring.DecodeCoeffs(pointer, N, 1, P.value[0].Coeffs, data)

	return nil
}

func (ciphertext *Ciphertext) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(ciphertext.Value()[0].Coeffs[0]))
	numberModulies := uint64(len(ciphertext.Value()[0].Coeffs))
	degree := ciphertext.Degree()

	if numberModulies > 0xFF {
		return nil, errors.New("error : ciphertext numberModulies overflow uint8")
	}

	if degree > 0xFF {
		return nil, errors.New("error : ciphertext degree overflow uint8")
	}

	data := make([]byte, 3+((N*numberModulies*(degree+1))<<3))

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModulies)
	data[2] = uint8(degree)

	pointer := uint64(3)

	for i := uint64(0); i < degree+1; i++ {
		if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, ciphertext.Value()[i].Coeffs, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModulies := uint64(data[1])
	degree := uint64(data[2])

	ciphertext.SetValue(make([]*ring.Poly, degree+1))

	for i := uint64(0); i < degree+1; i++ {
		ciphertext.Value()[i] = new(ring.Poly)
		ciphertext.Value()[i].Coeffs = make([][]uint64, numberModulies)
	}

	pointer := uint64(3)

	if ((uint64(len(data)) - pointer) >> 3) != N*numberModulies*(degree+1) {
		return errors.New("error : invalid ciphertext encoding")
	}

	for x := uint64(0); x < degree+1; x++ {
		pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, ciphertext.Value()[x].Coeffs, data)
	}

	return nil
}

func (sk *SecretKey) MarshalBinary() ([]byte, error) {
	var err error

	N := uint64(len(sk.sk.Coeffs[0]))
	numberModulies := uint64(len(sk.sk.Coeffs))

	if numberModulies > 0xFF {
		return nil, errors.New("error : max degree uint16 overflow")
	}

	data := make([]byte, 2+((N*numberModulies)<<3))

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModulies)

	pointer := uint64(2)

	if _, err = ring.WriteCoeffsTo(pointer, N, numberModulies, sk.sk.Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

func (sk *SecretKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModulies := uint64(data[1])

	sk.sk = new(ring.Poly)
	sk.sk.Coeffs = make([][]uint64, numberModulies)

	pointer := uint64(2)

	if ((uint64(len(data)) - pointer) >> 3) != (N * numberModulies) {
		return errors.New("error : invalid SecretKey encoding")
	}

	_, _ = ring.DecodeCoeffs(pointer, N, numberModulies, sk.sk.Coeffs, data)

	return nil
}

// PK
func (pk *PublicKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(pk.pk[0].Coeffs[0]))
	numberModulies := uint64(len(pk.pk[0].Coeffs))

	if numberModulies > 0xFF {
		return nil, errors.New("error : max degree uint16 overflow")
	}

	data := make([]byte, 2+((N*numberModulies)<<4))

	data[0] = uint8((bits.Len64(uint64(N)) - 1))
	data[1] = uint8(numberModulies)

	pointer := uint64(2)

	if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, pk.pk[0].Coeffs, data); err != nil {
		return nil, err
	}

	if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, pk.pk[1].Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

func (pk *PublicKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModulies := uint64(data[1])

	pk.pk[0] = new(ring.Poly)
	pk.pk[0].Coeffs = make([][]uint64, numberModulies)
	pk.pk[1] = new(ring.Poly)
	pk.pk[1].Coeffs = make([][]uint64, numberModulies)

	pointer := uint64(2)

	if ((uint64(len(data)) - pointer) >> 4) != (N * numberModulies) {
		return errors.New("error : invalid PublicKey encoding")
	}

	pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, pk.pk[0].Coeffs, data)
	pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, pk.pk[1].Coeffs, data)

	return nil
}

func (evaluationkey *EvaluationKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(evaluationkey.evakey[0].evakey[0][0][0].Coeffs[0]))
	numberModulies := uint64(len(evaluationkey.evakey[0].evakey[0][0][0].Coeffs))
	decomposition := numberModulies
	bitDecomp := evaluationkey.evakey[0].bitDecomp

	maxDegree := uint64(len(evaluationkey.evakey))

	if numberModulies > 0xFF {
		return nil, errors.New("error : max number modulies uint16 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint16 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint16 overflow")
	}

	if maxDegree > 0xFF {
		return nil, errors.New("error : max degree uint16 overflow")
	}

	var dataLen uint64
	dataLen = 5
	for i := uint64(0); i < maxDegree; i++ {
		for j := uint64(0); j < decomposition; j++ {
			dataLen += 1                                                                                                       //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * numberModulies * decomposition * maxDegree * uint64(len(evaluationkey.evakey[i].evakey[j])) // nb coefficients * 8
		}
	}

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModulies)
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, evaluationkey.evakey[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, evaluationkey.evakey[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	return data, nil
}

func (evaluationkey *EvaluationKey) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModulies := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])
	maxDegree := uint64(data[4])

	evaluationkey.evakey = make([]*SwitchingKey, maxDegree)

	pointer := uint64(5)
	var bitLog uint64
	for i := uint64(0); i < maxDegree; i++ {

		evaluationkey.evakey[i] = new(SwitchingKey)
		evaluationkey.evakey[i].bitDecomp = bitDecomp
		evaluationkey.evakey[i].evakey = make([][][2]*ring.Poly, decomposition)

		for j := uint64(0); j < decomposition; j++ {

			bitLog = uint64(data[pointer])
			pointer += 1

			evaluationkey.evakey[i].evakey[j] = make([][2]*ring.Poly, bitLog)

			for x := uint64(0); x < bitLog; x++ {

				evaluationkey.evakey[i].evakey[j][x][0] = new(ring.Poly)
				evaluationkey.evakey[i].evakey[j][x][0].Coeffs = make([][]uint64, numberModulies)
				pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, evaluationkey.evakey[i].evakey[j][x][0].Coeffs, data)

				evaluationkey.evakey[i].evakey[j][x][1] = new(ring.Poly)
				evaluationkey.evakey[i].evakey[j][x][1].Coeffs = make([][]uint64, numberModulies)
				pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, evaluationkey.evakey[i].evakey[j][x][1].Coeffs, data)
			}
		}
	}

	return nil
}

func (rotationkey *RotationKeys) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(rotationkey.bfvcontext.n)
	numberModulies := uint64(len(rotationkey.bfvcontext.contextQ.Modulus))
	decomposition := numberModulies
	bitDecomp := rotationkey.bitDecomp
	mappingRow := 0
	mappingColL := []uint64{}
	mappingColR := []uint64{}

	if numberModulies > 0xFF {
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
				dataLen += 2 * 8 * N * numberModulies * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}

		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColR = append(mappingColR, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                                   //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * numberModulies * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}
	}

	if rotationkey.evakey_rot_row != nil {
		mappingRow = 1
		for j := uint64(0); j < decomposition; j++ {
			dataLen += 1                                                                                              //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * numberModulies * decomposition * uint64(len(rotationkey.evakey_rot_row.evakey[j])) // nb coefficients * 8
		}
	}

	dataLen += uint64(len(mappingColL)+len(mappingColR)) << 2 // size needed to encode what rotation are present

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModulies)
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data); err != nil {
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data); err != nil {
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModulies, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	return data, nil
}

func (rotationkey *RotationKeys) UnMarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModulies := uint64(data[1])
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
				rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs = make([][]uint64, numberModulies)
				pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data)

				rotationkey.evakey_rot_row.evakey[j][x][1] = new(ring.Poly)
				rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs = make([][]uint64, numberModulies)
				pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data)
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
					rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs = make([][]uint64, numberModulies)
					pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_L[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs = make([][]uint64, numberModulies)
					pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data)
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
					rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs = make([][]uint64, numberModulies)
					pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_R[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs = make([][]uint64, numberModulies)
					pointer, _ = ring.DecodeCoeffs(pointer, N, numberModulies, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	return nil
}
