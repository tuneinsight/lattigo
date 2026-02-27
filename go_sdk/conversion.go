package main

/*
#include "../../fhe_types_v2.h"
*/
import "C"
import (
	"unsafe"

	"github.com/cipherflow-fhe/lattigo/bfv"
	"github.com/cipherflow-fhe/lattigo/ckks"
	"github.com/cipherflow-fhe/lattigo/ring"
	"github.com/cipherflow-fhe/lattigo/rlwe"
)

//export BfvComponentNttInplace
func BfvComponentNttInplace(parameter_handle uint64, coeff *C.ulong, lvl_idx int) {
	param := get_object[bfv.Parameters](parameter_handle)
	ringq := param.RingQ()
	data_slice := unsafe.Slice((*uint64)(coeff), ringq.N)

	if lvl_idx < param.QCount() {
		ring.NTT(data_slice, data_slice, ringq.N, ringq.NttPsi[lvl_idx], ringq.Modulus[lvl_idx], ringq.MredParams[lvl_idx], ringq.BredParams[lvl_idx])
	} else {
		ringp := param.RingP()
		sp_lvl_idx := lvl_idx - param.QCount()
		ring.NTT(data_slice, data_slice, ringp.N, ringp.NttPsi[sp_lvl_idx], ringp.Modulus[sp_lvl_idx], ringp.MredParams[sp_lvl_idx], ringp.BredParams[sp_lvl_idx])
	}

}

//export BfvComponentInvNttInplace
func BfvComponentInvNttInplace(parameter_handle uint64, coeff *C.ulong, lvl_idx int) {
	param := get_object[bfv.Parameters](parameter_handle)
	ringq := param.RingQ()
	data_slice := unsafe.Slice((*uint64)(coeff), ringq.N)

	if lvl_idx < param.QCount() {
		ring.InvNTT(data_slice, data_slice, ringq.N, ringq.NttPsiInv[lvl_idx], ringq.NttNInv[lvl_idx], ringq.Modulus[lvl_idx], ringq.MredParams[lvl_idx])
	} else {
		ringp := param.RingP()
		sp_lvl_idx := lvl_idx - param.QCount()
		ring.InvNTT(data_slice, data_slice, ringp.N, ringp.NttPsiInv[sp_lvl_idx], ringp.NttNInv[sp_lvl_idx], ringp.Modulus[sp_lvl_idx], ringp.MredParams[sp_lvl_idx])
	}
}

//export CkksComponentNttInplace
func CkksComponentNttInplace(parameter_handle uint64, coeff *C.ulong, lvl_idx int) {
	param := get_object[ckks.Parameters](parameter_handle)
	ringq := param.RingQ()
	data_slice := unsafe.Slice((*uint64)(coeff), ringq.N)

	if lvl_idx < param.QCount() {
		ring.NTT(data_slice, data_slice, ringq.N, ringq.NttPsi[lvl_idx], ringq.Modulus[lvl_idx], ringq.MredParams[lvl_idx], ringq.BredParams[lvl_idx])
	} else {
		ringp := param.RingP()
		sp_lvl_idx := lvl_idx - param.QCount()
		ring.NTT(data_slice, data_slice, ringp.N, ringp.NttPsi[sp_lvl_idx], ringp.Modulus[sp_lvl_idx], ringp.MredParams[sp_lvl_idx], ringp.BredParams[sp_lvl_idx])
	}

}

//export CkksComponentInvNttInplace
func CkksComponentInvNttInplace(parameter_handle uint64, coeff *C.ulong, lvl_idx int) {
	param := get_object[ckks.Parameters](parameter_handle)
	ringq := param.RingQ()
	data_slice := unsafe.Slice((*uint64)(coeff), ringq.N)

	if lvl_idx < param.QCount() {
		ring.InvNTT(data_slice, data_slice, ringq.N, ringq.NttPsiInv[lvl_idx], ringq.NttNInv[lvl_idx], ringq.Modulus[lvl_idx], ringq.MredParams[lvl_idx])
	} else {
		ringp := param.RingP()
		sp_lvl_idx := lvl_idx - param.QCount()
		ring.InvNTT(data_slice, data_slice, ringp.N, ringp.NttPsiInv[sp_lvl_idx], ringp.NttNInv[sp_lvl_idx], ringp.Modulus[sp_lvl_idx], ringp.MredParams[sp_lvl_idx])
	}
}

//export BfvComponentMulByPow2Inplace
func BfvComponentMulByPow2Inplace(parameter_handle uint64, coeff *C.ulong, lvl_idx int, pow2 int) {
	param := get_object[bfv.Parameters](parameter_handle)
	ringq := param.RingQ()
	data_slice := unsafe.Slice((*uint64)(coeff), ringq.N)

	if lvl_idx < param.QCount() {
		ring.MFormVec(data_slice, data_slice, ringq.Modulus[lvl_idx], ringq.BredParams[lvl_idx])
		ring.MulByPow2Vec(data_slice, data_slice, pow2, ringq.Modulus[lvl_idx], ringq.MredParams[lvl_idx])
	} else {
		ringp := param.RingP()
		sp_lvl_idx := lvl_idx - param.QCount()
		ring.MFormVec(data_slice, data_slice, ringp.Modulus[sp_lvl_idx], ringp.BredParams[sp_lvl_idx])
		ring.MulByPow2Vec(data_slice, data_slice, pow2, ringp.Modulus[sp_lvl_idx], ringp.MredParams[sp_lvl_idx])
	}

}

//export CkksComponentMulByPow2Inplace
func CkksComponentMulByPow2Inplace(parameter_handle uint64, coeff *C.ulong, lvl_idx int, pow2 int) {
	param := get_object[ckks.Parameters](parameter_handle)
	ringq := param.RingQ()
	data_slice := unsafe.Slice((*uint64)(coeff), ringq.N)

	if lvl_idx < param.QCount() {
		ring.MFormVec(data_slice, data_slice, ringq.Modulus[lvl_idx], ringq.BredParams[lvl_idx])
		ring.MulByPow2Vec(data_slice, data_slice, pow2, ringq.Modulus[lvl_idx], ringq.MredParams[lvl_idx])
	} else {
		ringp := param.RingP()
		sp_lvl_idx := lvl_idx - param.QCount()
		ring.MFormVec(data_slice, data_slice, ringp.Modulus[sp_lvl_idx], ringp.BredParams[sp_lvl_idx])
		ring.MulByPow2Vec(data_slice, data_slice, pow2, ringp.Modulus[sp_lvl_idx], ringp.MredParams[sp_lvl_idx])
	}

}

//export BfvPlaintextMulInvMFormAndMulByPow2
func BfvPlaintextMulInvMFormAndMulByPow2(parameter_handle uint64, plaintext_mul_handle uint64, pow2 int) {
	param := get_object[bfv.Parameters](parameter_handle)
	plaintext_mul := get_object[bfv.PlaintextMul](plaintext_mul_handle)
	param.RingQ().InvMFormAndMulByPow2(plaintext_mul.Value, pow2, plaintext_mul.Value)
}

//export CkksPlaintextMulInvMFormAndMulByPow2
func CkksPlaintextMulInvMFormAndMulByPow2(parameter_handle uint64, plaintext_mul_handle uint64, pow2 int) {
	param := get_object[ckks.Parameters](parameter_handle)
	plaintext_mul := get_object[ckks.PlaintextMul](plaintext_mul_handle)
	param.RingQ().InvMFormAndMulByPow2(plaintext_mul.Value, pow2, plaintext_mul.Value)
}

//export BfvRlkInvMForm
func BfvRlkInvMForm(parameter_handle uint64, relin_key_handle uint64) {
	param := get_object[bfv.Parameters](parameter_handle)
	ringq := param.RingQ()
	ringp := param.RingP()
	relin_key := get_object[rlwe.RelinearizationKey](relin_key_handle)

	for _, swk := range relin_key.Keys {
		for j := range swk.Value {
			for k := range swk.Value[j][0].Value {
				ringq.InvMForm(swk.Value[j][0].Value[k].Q, swk.Value[j][0].Value[k].Q)
				ringp.InvMForm(swk.Value[j][0].Value[k].P, swk.Value[j][0].Value[k].P)
			}
		}
	}
}

//export BfvRlkInvMFormAndMulByPow2
func BfvRlkInvMFormAndMulByPow2(parameter_handle uint64, relin_key_handle uint64, pow2 int) {
	param := get_object[bfv.Parameters](parameter_handle)
	ringq := param.RingQ()
	ringp := param.RingP()
	relin_key := get_object[rlwe.RelinearizationKey](relin_key_handle)

	for _, swk := range relin_key.Keys {
		for j := range swk.Value {
			for k := range swk.Value[j][0].Value {
				ringq.InvMFormAndMulByPow2(swk.Value[j][0].Value[k].Q, pow2, swk.Value[j][0].Value[k].Q)
				ringp.InvMFormAndMulByPow2(swk.Value[j][0].Value[k].P, pow2, swk.Value[j][0].Value[k].P)
			}
		}
	}
}

//export BfvGlkInvMForm
func BfvGlkInvMForm(parameter_handle uint64, galois_key_handle uint64) {
	param := get_object[bfv.Parameters](parameter_handle)
	ringq := param.RingQ()
	ringp := param.RingP()
	galois_key := get_object[rlwe.RotationKeySet](galois_key_handle)

	for _, swk := range galois_key.Keys {
		for i := 0; i < len(swk.Value); i++ {
			for j := 0; j < len(swk.Value[i][0].Value); j++ {
				ringq.InvMForm(swk.Value[i][0].Value[j].Q, swk.Value[i][0].Value[j].Q)
				ringp.InvMForm(swk.Value[i][0].Value[j].P, swk.Value[i][0].Value[j].P)
			}
		}
	}
}

//export BfvGlkInvMFormAndMulByPow2
func BfvGlkInvMFormAndMulByPow2(parameter_handle uint64, galois_key_handle uint64, pow2 int) {
	param := get_object[bfv.Parameters](parameter_handle)
	ringq := param.RingQ()
	ringp := param.RingP()
	galois_key := get_object[rlwe.RotationKeySet](galois_key_handle)

	for _, swk := range galois_key.Keys {
		for i := 0; i < len(swk.Value); i++ {
			for j := 0; j < len(swk.Value[i][0].Value); j++ {
				ringq.InvMFormAndMulByPow2(swk.Value[i][0].Value[j].Q, pow2, swk.Value[i][0].Value[j].Q)
				ringp.InvMFormAndMulByPow2(swk.Value[i][0].Value[j].P, pow2, swk.Value[i][0].Value[j].P)
			}
		}
	}
}

//export CkksRlkInvMForm
func CkksRlkInvMForm(parameter_handle uint64, relin_key_handle uint64) {
	param := get_object[ckks.Parameters](parameter_handle)
	ringq := param.RingQ()
	ringp := param.RingP()
	relin_key := get_object[rlwe.RelinearizationKey](relin_key_handle)

	for _, swk := range relin_key.Keys {
		for j := range swk.Value {
			for k := range swk.Value[j][0].Value {
				ringq.InvMForm(swk.Value[j][0].Value[k].Q, swk.Value[j][0].Value[k].Q)
				ringp.InvMForm(swk.Value[j][0].Value[k].P, swk.Value[j][0].Value[k].P)
			}
		}
	}
}

//export CkksRlkInvMFormAndMulByPow2
func CkksRlkInvMFormAndMulByPow2(parameter_handle uint64, relin_key_handle uint64, pow2 int) {
	param := get_object[ckks.Parameters](parameter_handle)
	ringq := param.RingQ()
	ringp := param.RingP()
	relin_key := get_object[rlwe.RelinearizationKey](relin_key_handle)

	for _, swk := range relin_key.Keys {
		for j := range swk.Value {
			for k := range swk.Value[j][0].Value {
				ringq.InvMFormAndMulByPow2(swk.Value[j][0].Value[k].Q, pow2, swk.Value[j][0].Value[k].Q)
				ringp.InvMFormAndMulByPow2(swk.Value[j][0].Value[k].P, pow2, swk.Value[j][0].Value[k].P)
			}
		}
	}
}

//export CkksGlkInvMForm
func CkksGlkInvMForm(parameter_handle uint64, galois_key_handle uint64) {
	param := get_object[ckks.Parameters](parameter_handle)
	ringq := param.RingQ()
	ringp := param.RingP()
	galois_key := get_object[rlwe.RotationKeySet](galois_key_handle)

	for _, swk := range galois_key.Keys {
		for i := 0; i < len(swk.Value); i++ {
			for j := 0; j < len(swk.Value[i][0].Value); j++ {
				ringq.InvMForm(swk.Value[i][0].Value[j].Q, swk.Value[i][0].Value[j].Q)
				ringp.InvMForm(swk.Value[i][0].Value[j].P, swk.Value[i][0].Value[j].P)
			}
		}
	}
}

//export CkksGlkInvMFormAndMulByPow2
func CkksGlkInvMFormAndMulByPow2(parameter_handle uint64, galois_key_handle uint64, pow2 int) {
	param := get_object[ckks.Parameters](parameter_handle)
	ringq := param.RingQ()
	ringp := param.RingP()
	galois_key := get_object[rlwe.RotationKeySet](galois_key_handle)

	for _, swk := range galois_key.Keys {
		for i := 0; i < len(swk.Value); i++ {
			for j := 0; j < len(swk.Value[i][0].Value); j++ {
				ringq.InvMFormAndMulByPow2(swk.Value[i][0].Value[j].Q, pow2, swk.Value[i][0].Value[j].Q)
				ringp.InvMFormAndMulByPow2(swk.Value[i][0].Value[j].P, pow2, swk.Value[i][0].Value[j].P)
			}
		}
	}
}

func set_switching_key_n_mform_bits(param *rlwe.Parameters, swk *rlwe.SwitchingKey, n_mform_bits int) {
	ringq := param.RingQ()
	ringp := param.RingP()
	diff := n_mform_bits - swk.NMFormBits

	if diff == 0 {
		return
	} else if diff == -64 {
		for j := range swk.Value {
			for k := range swk.Value[j][0].Value {
				ringq.InvMForm(swk.Value[j][0].Value[k].Q, swk.Value[j][0].Value[k].Q)
				ringp.InvMForm(swk.Value[j][0].Value[k].P, swk.Value[j][0].Value[k].P)
			}
		}
	} else if diff > 0 {
		for j := range swk.Value {
			for k := range swk.Value[j][0].Value {
				ringq.MulByPow2(swk.Value[j][0].Value[k].Q, diff, swk.Value[j][0].Value[k].Q)
				ringp.MulByPow2(swk.Value[j][0].Value[k].P, diff, swk.Value[j][0].Value[k].P)
			}
		}
	} else {
		for j := range swk.Value {
			for k := range swk.Value[j][0].Value {
				ringq.InvMFormAndMulByPow2(swk.Value[j][0].Value[k].Q, 64+diff, swk.Value[j][0].Value[k].Q)
				ringp.InvMFormAndMulByPow2(swk.Value[j][0].Value[k].P, 64+diff, swk.Value[j][0].Value[k].P)
			}
		}
	}
	swk.NMFormBits = n_mform_bits
}

//export SetBfvRlkNMFormBits
func SetBfvRlkNMFormBits(parameter_handle uint64, relin_key_handle uint64, n_mform_bits int) {
	param := get_object[bfv.Parameters](parameter_handle)
	relin_key := get_object[rlwe.RelinearizationKey](relin_key_handle)

	for _, swk := range relin_key.Keys {
		set_switching_key_n_mform_bits(&param.Parameters, swk, n_mform_bits)
	}
}

//export SetCkksRlkNMFormBits
func SetCkksRlkNMFormBits(parameter_handle uint64, relin_key_handle uint64, n_mform_bits int) {
	param := get_object[ckks.Parameters](parameter_handle)
	relin_key := get_object[rlwe.RelinearizationKey](relin_key_handle)

	for _, swk := range relin_key.Keys {
		set_switching_key_n_mform_bits(&param.Parameters, swk, n_mform_bits)
	}
}

//export SetBfvGlkNMFormBits
func SetBfvGlkNMFormBits(parameter_handle uint64, galois_key_handle uint64, n_mform_bits int) {
	param := get_object[bfv.Parameters](parameter_handle)
	galois_key_set := get_object[rlwe.RotationKeySet](galois_key_handle)

	for _, swk := range galois_key_set.Keys {
		set_switching_key_n_mform_bits(&param.Parameters, swk, n_mform_bits)
	}
}

//export SetBfvGlkNMFormBitsForGaloisElement
func SetBfvGlkNMFormBitsForGaloisElement(parameter_handle uint64, galois_key_handle uint64, galois_element uint64, n_mform_bits int) {
	param := get_object[bfv.Parameters](parameter_handle)
	galois_key_set := get_object[rlwe.RotationKeySet](galois_key_handle)

	set_switching_key_n_mform_bits(&param.Parameters, galois_key_set.Keys[galois_element], n_mform_bits)
}

//export SetCkksSwkNMFormBits
func SetCkksSwkNMFormBits(parameter_handle uint64, switching_key_handle uint64, n_mform_bits int) {
	param := get_object[ckks.Parameters](parameter_handle)
	switching_key := get_object[rlwe.SwitchingKey](switching_key_handle)
	set_switching_key_n_mform_bits(&param.Parameters, switching_key, n_mform_bits)
}

//export SetCkksGlkNMFormBits
func SetCkksGlkNMFormBits(parameter_handle uint64, galois_key_handle uint64, n_mform_bits int) {
	param := get_object[ckks.Parameters](parameter_handle)
	galois_key_set := get_object[rlwe.RotationKeySet](galois_key_handle)

	for _, swk := range galois_key_set.Keys {
		set_switching_key_n_mform_bits(&param.Parameters, swk, n_mform_bits)
	}
}

//export SetCkksGlkNMFormBitsForGaloisElement
func SetCkksGlkNMFormBitsForGaloisElement(parameter_handle uint64, galois_key_handle uint64, galois_element uint64, n_mform_bits int) {
	param := get_object[ckks.Parameters](parameter_handle)
	galois_key_set := get_object[rlwe.RotationKeySet](galois_key_handle)

	set_switching_key_n_mform_bits(&param.Parameters, galois_key_set.Keys[galois_element], n_mform_bits)

}
