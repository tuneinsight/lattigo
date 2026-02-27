package rlwe

import (
	"bytes"
	"encoding/binary"
	"math/bits"

	"github.com/cipherflow-fhe/lattigo/ring"
)

func component_to_bytes(src []uint64, bit_length uint8, writer *bytes.Buffer) {
	byte_mask := [9]byte{0b0, 0b1, 0b11, 0b111, 0b1111, 0b11111, 0b111111, 0b1111111, 0b11111111}
	N := len(src)
	var b byte             // the current byte, converted from an integer, but not yet written to writer
	bit_offset := uint8(0) // [0, 7], number of available bits in b
	for i := 0; i < N; i++ {
		n_low_bit := 8 - bit_offset                // [1, 8]
		n_high_bit := (bit_length - n_low_bit) % 8 // [0, 7]
		n_full_byte := (bit_length - n_low_bit - n_high_bit) / 8

		x := src[i]

		x_low := byte(x) & byte_mask[n_low_bit]
		b |= (x_low << byte(bit_offset))
		writer.WriteByte(b)
		x >>= uint64(n_low_bit)

		for j := uint8(0); j < n_full_byte; j++ {
			writer.WriteByte(byte(x))
			x >>= 8
		}

		b = byte(x)
		bit_offset = n_high_bit
	}
	if bit_offset > 0 {
		writer.WriteByte(b)
	}
}

func bytes_to_component(reader *bytes.Reader, bit_length uint8, dest []uint64) {
	byte_mask := [9]byte{0b0, 0b1, 0b11, 0b111, 0b1111, 0b11111, 0b111111, 0b1111111, 0b11111111}
	N := len(dest)
	var b byte             // the current byte, read from reader, but not fully converted to integer
	bit_offset := uint8(0) // [0, 7], number of processed bits in b
	b, _ = reader.ReadByte()
	for i := 0; i < N; i++ {
		n_low_bit := 8 - bit_offset                // [1, 8]
		n_high_bit := (bit_length - n_low_bit) % 8 // [0, 7]
		n_full_byte := (bit_length - n_low_bit - n_high_bit) / 8

		var x uint64

		x = uint64(b >> uint64(bit_offset))
		b, _ = reader.ReadByte()

		n_shift := n_low_bit
		for j := uint8(0); j < n_full_byte; j++ {
			x |= (uint64(b) << n_shift)
			b, _ = reader.ReadByte()
			n_shift += 8
		}

		x |= uint64(b&byte_mask[n_high_bit]) << n_shift
		bit_offset = n_high_bit

		dest[i] = x
	}
	if bit_offset == 0 {
		reader.UnreadByte()
	}
}

func poly_to_bytes(poly *ring.Poly, q_bit_lengths []uint8, n_drop_bit uint8, writer *bytes.Buffer) {
	N := poly.N()
	level := poly.Level()

	binary.Write(writer, binary.LittleEndian, poly.IsNTT)
	binary.Write(writer, binary.LittleEndian, poly.IsMForm)
	writer.WriteByte(n_drop_bit)

	for j, component := range poly.Coeffs {
		bit_length := q_bit_lengths[j]
		if level == 0 && n_drop_bit != 0 {
			bit_length -= n_drop_bit
			component_shift := make([]uint64, N)
			for k := 0; k < N; k++ {
				component_shift[k] = (component[k] >> n_drop_bit) + ((component[k] >> (n_drop_bit - 1)) & 1)
			}
			component_to_bytes(component_shift, bit_length, writer)
		} else {
			component_to_bytes(component, bit_length, writer)
		}
	}
}

func bytes_to_poly(reader *bytes.Reader, N int, level int, q_bit_lengths []uint8, poly *ring.Poly) {
	binary.Read(reader, binary.LittleEndian, &poly.IsNTT)
	binary.Read(reader, binary.LittleEndian, &poly.IsMForm)
	n_drop_bit, _ := reader.ReadByte()

	poly.Buff = make([]uint64, N*(level+1))
	poly.Coeffs = make([][]uint64, level+1)
	for j := 0; j < level+1; j++ {
		poly.Coeffs[j] = poly.Buff[int(j)*N : (int(j)+1)*N]
		bit_length := q_bit_lengths[j]
		if level == 0 && n_drop_bit != 0 {
			bit_length -= n_drop_bit
		}
		bytes_to_component(reader, bit_length, poly.Coeffs[j])
		if level == 0 && n_drop_bit != 0 {
			for k := 0; k < N; k++ {
				poly.Coeffs[j][k] <<= n_drop_bit
			}
		}
	}
}

func CiphertextToBytes(src *Ciphertext, param *Parameters, n_drop_bit_0 int, n_drop_bit_1 int, writer *bytes.Buffer) {
	// N: 4B, n_poly 1B, level 1B, q_bit_lengths (level+1)B
	// per poly: { IsNTT: 1B, IsMForm: 1B, n_drop_bit: 1B }
	// data_length := 2 + level + 1 + 8*n_poly + ct_data_length

	param_q := param.Q()
	N := param.N()
	n_poly := src.Degree() + 1
	level := src.Level()
	q_bit_lengths := make([]uint8, level+1)
	for j := 0; j < level+1; j++ {
		q_bit_lengths[j] = uint8(bits.Len64(param_q[j]))
	}
	n_drop_bits := []uint8{uint8(n_drop_bit_0), uint8(n_drop_bit_1)}

	binary.Write(writer, binary.LittleEndian, uint32(N))
	writer.WriteByte(byte(n_poly))
	writer.WriteByte(byte(level))
	for j := 0; j < level+1; j++ {
		writer.WriteByte(byte(q_bit_lengths[j]))
	}

	for i, poly := range src.Value {
		poly_to_bytes(poly, q_bit_lengths, n_drop_bits[i], writer)
	}
}

func BytesToCiphertext(reader *bytes.Reader) Ciphertext {
	var N32 uint32
	binary.Read(reader, binary.LittleEndian, &N32)
	N := int(N32)
	n_poly, _ := reader.ReadByte()
	level, _ := reader.ReadByte()
	q_bit_lengths := make([]uint8, level+1)
	for j := uint8(0); j < level+1; j++ {
		q_bit_lengths[j], _ = reader.ReadByte()
	}

	var ct Ciphertext
	ct.Value = make([]*ring.Poly, n_poly)
	for i := uint8(0); i < n_poly; i++ {
		poly := new(ring.Poly)
		ct.Value[i] = poly
		bytes_to_poly(reader, N, int(level), q_bit_lengths, poly)
	}

	return ct
}

func CompressedCiphertextToBytes(src *CompressedCiphertext, param *Parameters, writer *bytes.Buffer) {
	param_q := param.Q()
	N := param.N()
	level := src.Level()
	q_bit_lengths := make([]uint8, level+1)
	for j := 0; j < level+1; j++ {
		q_bit_lengths[j] = uint8(bits.Len64(param_q[j]))
	}

	binary.Write(writer, binary.LittleEndian, uint32(N))
	writer.WriteByte(byte(level))
	for j := 0; j < level+1; j++ {
		writer.WriteByte(byte(q_bit_lengths[j]))
	}

	for i := 0; i < 64; i++ {
		writer.WriteByte(src.Seed[i])
	}
	poly_to_bytes(src.Value, q_bit_lengths, 0, writer)
}

func BytesToCompressedCiphertext(reader *bytes.Reader) CompressedCiphertext {
	var N32 uint32
	binary.Read(reader, binary.LittleEndian, &N32)
	N := int(N32)
	level, _ := reader.ReadByte()
	q_bit_lengths := make([]uint8, level+1)
	for j := uint8(0); j < level+1; j++ {
		q_bit_lengths[j], _ = reader.ReadByte()
	}

	var ct CompressedCiphertext

	ct.Seed = make([]byte, 64)
	for i := 0; i < 64; i++ {
		ct.Seed[i], _ = reader.ReadByte()
	}
	ct.Value = new(ring.Poly)
	bytes_to_poly(reader, N, int(level), q_bit_lengths, ct.Value)

	return ct
}

func SecretKeyToBytes(src *SecretKey, param *Parameters, writer *bytes.Buffer) {
	param_q := param.Q()
	param_p := param.P()
	N := param.N()
	level_q := src.LevelQ()
	level_p := src.LevelP()
	q_bit_lengths := make([]uint8, level_q+1)
	p_bit_lengths := make([]uint8, level_p+1)
	for i := 0; i < level_q+1; i++ {
		q_bit_lengths[i] = uint8(bits.Len64(param_q[i]))
	}
	for i := 0; i < level_p+1; i++ {
		p_bit_lengths[i] = uint8(bits.Len64(param_p[i]))
	}

	binary.Write(writer, binary.LittleEndian, uint32(N))
	writer.WriteByte(byte(level_q))
	writer.WriteByte(byte(level_p))
	for _, q_bit_length := range q_bit_lengths {
		writer.WriteByte(byte(q_bit_length))
	}
	for _, p_bit_length := range p_bit_lengths {
		writer.WriteByte(byte(p_bit_length))
	}
	poly_to_bytes(src.Value.Q, q_bit_lengths, 0, writer)
	poly_to_bytes(src.Value.P, p_bit_lengths, 0, writer)
}

func BytesToSecretKey(reader *bytes.Reader) SecretKey {
	var N32 uint32
	binary.Read(reader, binary.LittleEndian, &N32)
	N := int(N32)
	level_q, _ := reader.ReadByte()
	level_p, _ := reader.ReadByte()
	q_bit_lengths := make([]uint8, level_q+1)
	for i := uint8(0); i < level_q+1; i++ {
		q_bit_lengths[i], _ = reader.ReadByte()
	}
	p_bit_lengths := make([]uint8, level_p+1)
	for i := uint8(0); i < level_p+1; i++ {
		p_bit_lengths[i], _ = reader.ReadByte()
	}

	var sk SecretKey
	poly_qp := &sk.Value
	poly_qp.Q = new(ring.Poly)
	bytes_to_poly(reader, N, int(level_q), q_bit_lengths, poly_qp.Q)
	poly_qp.P = new(ring.Poly)
	bytes_to_poly(reader, N, int(level_p), p_bit_lengths, poly_qp.P)

	return sk
}

func CiphertextQPToBytes(src *CiphertextQP, param *Parameters, writer *bytes.Buffer) {
	N := param.N()
	param_q := param.Q()
	param_p := param.P()
	level_q := src.LevelQ()
	level_p := src.LevelP()
	q_bit_lengths := make([]uint8, level_q+1)
	p_bit_lengths := make([]uint8, level_p+1)
	for i := 0; i < level_q+1; i++ {
		q_bit_lengths[i] = uint8(bits.Len64(param_q[i]))
	}
	for i := 0; i < level_p+1; i++ {
		p_bit_lengths[i] = uint8(bits.Len64(param_p[i]))
	}

	binary.Write(writer, binary.LittleEndian, uint32(N))
	writer.WriteByte(byte(level_q))
	writer.WriteByte(byte(level_p))
	for i := 0; i < level_q+1; i++ {
		writer.WriteByte(byte(q_bit_lengths[i]))
	}
	for i := 0; i < level_p+1; i++ {
		writer.WriteByte(byte(p_bit_lengths[i]))
	}

	for i := 0; i < 64; i++ {
		writer.WriteByte(src.Seed[i])
	}
	poly_to_bytes(src.Value[0].Q, q_bit_lengths, 0, writer)
	poly_to_bytes(src.Value[0].P, p_bit_lengths, 0, writer)
}

func BytesToCiphertextQP(reader *bytes.Reader) CiphertextQP {
	var N32 uint32
	binary.Read(reader, binary.LittleEndian, &N32)
	N := int(N32)
	level_q, _ := reader.ReadByte()
	level_p, _ := reader.ReadByte()
	q_bit_lengths := make([]uint8, level_q+1)
	for i := uint8(0); i < level_q+1; i++ {
		q_bit_lengths[i], _ = reader.ReadByte()
	}
	p_bit_lengths := make([]uint8, level_p+1)
	for i := uint8(0); i < level_p+1; i++ {
		p_bit_lengths[i], _ = reader.ReadByte()
	}

	var ct CiphertextQP

	ct.Compressed = true
	ct.Seed = make([]byte, 64)
	for i := 0; i < 64; i++ {
		ct.Seed[i], _ = reader.ReadByte()
	}
	ct.Value[0].Q = new(ring.Poly)
	bytes_to_poly(reader, N, int(level_q), q_bit_lengths, ct.Value[0].Q)
	ct.Value[0].P = new(ring.Poly)
	bytes_to_poly(reader, N, int(level_p), p_bit_lengths, ct.Value[0].P)

	return ct
}

func GadgetCiphertextToBytes(src *GadgetCiphertext, param *Parameters, writer *bytes.Buffer) {
	writer.WriteByte(byte(len(src.Value)))
	writer.WriteByte(byte(len(src.Value[0])))

	for _, x := range src.Value {
		for _, y := range x {
			CiphertextQPToBytes(&y, param, writer)
		}
	}
}

func BytesToGadgetCiphertext(reader *bytes.Reader) GadgetCiphertext {
	n0, _ := reader.ReadByte()
	n1, _ := reader.ReadByte()

	var gc GadgetCiphertext
	gc.Value = make([][]CiphertextQP, n0)
	for i := uint8(0); i < n0; i++ {
		gc.Value[i] = make([]CiphertextQP, n1)
		for j := uint8(0); j < n1; j++ {
			gc.Value[i][j] = BytesToCiphertextQP(reader)
		}
	}

	return gc
}

func SwitchingKeyToBytes(src *SwitchingKey, param *Parameters, writer *bytes.Buffer) {
	GadgetCiphertextToBytes(&src.GadgetCiphertext, param, writer)
}

func BytesToSwitchingKey(reader *bytes.Reader) *SwitchingKey {
	swk := new(SwitchingKey)
	swk.GadgetCiphertext = BytesToGadgetCiphertext(reader)
	swk.NMFormBits = 64
	return swk
}

func RelinearizationKeyToByte(src *RelinearizationKey, param *Parameters, writer *bytes.Buffer) {
	writer.WriteByte(byte(len(src.Keys)))

	for _, x := range src.Keys {
		SwitchingKeyToBytes(x, param, writer)
	}
}

func BytesToRelinearizationKey(reader *bytes.Reader) RelinearizationKey {
	n, _ := reader.ReadByte()

	var rlk RelinearizationKey
	rlk.Keys = make([]*SwitchingKey, n)
	for i := uint8(0); i < n; i++ {
		rlk.Keys[i] = BytesToSwitchingKey(reader)
	}

	return rlk
}

func RotationKeySetToBytes(src *RotationKeySet, param *Parameters, writer *bytes.Buffer) {
	n := uint16(len(src.Keys))
	binary.Write(writer, binary.LittleEndian, n)

	for step, key := range src.Keys {
		binary.Write(writer, binary.LittleEndian, step)
		GadgetCiphertextToBytes(&key.GadgetCiphertext, param, writer)
	}
}

func BytesToRotationKeySet(reader *bytes.Reader) RotationKeySet {
	var n uint16
	binary.Read(reader, binary.LittleEndian, &n)

	var glk RotationKeySet
	glk.Keys = make(map[uint64]*SwitchingKey)
	for i := uint16(0); i < n; i++ {
		var step uint64
		binary.Read(reader, binary.LittleEndian, &step)

		swk := new(SwitchingKey)
		swk.GadgetCiphertext = BytesToGadgetCiphertext(reader)
		glk.Keys[step] = swk
	}

	return glk
}
