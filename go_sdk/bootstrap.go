package main

/*
#include "../../fhe_types_v2.h"
*/
import "C"
import (
	"bytes"
	"encoding/binary"
	"unsafe"

	"github.com/cipherflow-fhe/lattigo/ckks"
	"github.com/cipherflow-fhe/lattigo/ckks/bootstrapping"
	"github.com/cipherflow-fhe/lattigo/rlwe"
)

type BtpParameterSet struct {
	SchemeParam ckks.Parameters
	BtpParam    bootstrapping.Parameters
}

type CkksBtpContext struct {
	CkksContext
	btp_parameter *bootstrapping.Parameters
	evk           *bootstrapping.EvaluationKeys
	bootstrapper  *bootstrapping.Bootstrapper
}

//export CreateCkksBtpParameter
func CreateCkksBtpParameter() uint64 {
	literal := bootstrapping.N16QP1546H192H32
	ckks_params := literal.SchemeParams
	btpParams := literal.BootstrappingParams

	params, err := ckks.NewParametersFromLiteral(ckks_params)
	if err != nil {
		panic(err)
	}

	btp_param := BtpParameterSet{params, btpParams}
	id := insert_object(&btp_param)
	return id
}

//export CreateCkksToyBtpParameter
func CreateCkksToyBtpParameter() uint64 {
	literal := bootstrapping.N16QP1546H192H32
	ckks_params := literal.SchemeParams
	btpParams := literal.BootstrappingParams

	// Important: the following N value is not secure, and should be used only for development.
	ckks_params.LogN = 13
	ckks_params.LogSlots = 12

	params, err := ckks.NewParametersFromLiteral(ckks_params)
	if err != nil {
		panic(err)
	}

	btp_param := BtpParameterSet{params, btpParams}
	id := insert_object(&btp_param)
	return id
}

//export GetCkksParameterFromBtpParameter
func GetCkksParameterFromBtpParameter(parameter_handle uint64) uint64 {
	param := get_object[BtpParameterSet](parameter_handle)
	id := insert_object(&param.SchemeParam)
	return id
}

//export CreateRandomCkksBtpContext
func CreateRandomCkksBtpContext(parameter_handle uint64) uint64 {
	param := get_object[BtpParameterSet](parameter_handle)
	var context CkksBtpContext

	context.parameter = &param.SchemeParam
	context.btp_parameter = &param.BtpParam
	context.encoder = ckks.NewEncoder(*context.parameter)
	context.kgen = ckks.NewKeyGenerator(*context.parameter)
	context.sk, context.pk = context.kgen.GenKeyPair()

	// Generate bootstrap evaluation keys first
	context.evk = new(bootstrapping.EvaluationKeys)
	*context.evk = bootstrapping.GenEvaluationKeys(*context.btp_parameter, *context.parameter, context.sk)

	// Point rlk and gk to evk's keys (shared references, not copies)
	context.rlk = context.evk.Rlk
	context.gk = context.evk.Rtks

	context.encryptor_pk = ckks.NewEncryptor(*context.parameter, context.pk)
	context.encryptor_sk = ckks.NewEncryptor(*context.parameter, context.sk)
	context.decryptor = ckks.NewDecryptor(*context.parameter, context.sk)
	context.evaluator = ckks.NewEvaluator(*context.parameter, rlwe.EvaluationKey{
		Rlk: context.rlk,
	})
	var err error
	context.bootstrapper, err = bootstrapping.NewBootstrapper(*context.parameter, *context.btp_parameter, *context.evk)
	if err != nil {
		panic(err)
	}

	id := insert_object(&context)
	return id
}

//export GenCkksBtpContextRotationKeys
func GenCkksBtpContextRotationKeys(context_handle uint64) {
	context := get_object[CkksBtpContext](context_handle)
	rots := make([]int, 2*context.parameter.LogN()-3)
	for i := 0; i < context.parameter.LogN()-1; i++ {
		rots[i] = (1 << i)
	}
	for i := 0; i < context.parameter.LogN()-2; i++ {
		rots[i+context.parameter.LogN()-1] = -1 * (1 << i)
	}

	for i := range rots {
		galEl := context.parameter.GaloisElementForColumnRotationBy(rots[i])
		if _, ok := context.gk.Keys[galEl]; !ok {
			context.gk.Keys[galEl] = context.kgen.GenSwitchingKeyForRotationBy(rots[i], context.sk)
		}
	}

	galEl := context.parameter.GaloisElementForRowRotation()
	if _, ok := context.gk.Keys[galEl]; !ok {
		context.gk.Keys[galEl] = context.kgen.GenSwitchingKeyForRowRotation(context.sk)
	}

	context.evaluator = context.evaluator.WithKey(rlwe.EvaluationKey{
		Rlk:  context.rlk,
		Rtks: context.gk,
	})
}

//export GenCkksBtpContextRotationKeysForRotations
func GenCkksBtpContextRotationKeysForRotations(context_handle uint64, rots *int32, rots_length int, include_swap_rows bool) {
	context := get_object[CkksBtpContext](context_handle)
	rots_slice := convert_slice(unsafe.Slice((*int32)(unsafe.Pointer(rots)), rots_length))

	for i := range rots_slice {
		galEl := context.parameter.GaloisElementForColumnRotationBy(rots_slice[i])
		if _, ok := context.gk.Keys[galEl]; !ok {
			context.gk.Keys[galEl] = context.kgen.GenSwitchingKeyForRotationBy(rots_slice[i], context.sk)
		}
	}

	if include_swap_rows {
		galEl := context.parameter.GaloisElementForRowRotation()
		if _, ok := context.gk.Keys[galEl]; !ok {
			context.gk.Keys[galEl] = context.kgen.GenSwitchingKeyForRowRotation(context.sk)
		}
	}

	context.evaluator = context.evaluator.WithKey(rlwe.EvaluationKey{
		Rlk:  context.rlk,
		Rtks: context.gk,
	})
}

//export ShallowCopyCkksBtpContext
func ShallowCopyCkksBtpContext(context_handle uint64) uint64 {
	var context_dest CkksBtpContext
	context_src := get_object[CkksBtpContext](context_handle)
	context_dest.parameter = context_src.parameter
	context_dest.btp_parameter = context_src.btp_parameter
	context_dest.sk = context_src.sk
	context_dest.pk = context_src.pk
	context_dest.rlk = context_src.rlk
	context_dest.gk = context_src.gk
	context_dest.evk = context_src.evk

	context_dest.encoder = context_src.encoder.ShallowCopy()

	if context_src.encryptor_pk != nil {
		context_dest.encryptor_pk = context_src.encryptor_pk.ShallowCopy()
	} else {
		context_dest.encryptor_pk = nil
	}

	if context_src.evaluator != nil {
		context_dest.evaluator = context_src.evaluator.ShallowCopy()
	}

	if context_src.sk != nil {
		context_dest.encryptor_sk = context_src.encryptor_sk.ShallowCopy()
		context_dest.decryptor = context_src.decryptor.ShallowCopy()
	} else {
		context_dest.encryptor_sk = nil
		context_dest.decryptor = nil
	}
	if context_src.evk != nil {
		context_dest.bootstrapper = context_src.bootstrapper.ShallowCopy()
	} else {
		context_dest.bootstrapper = nil
	}

	id := insert_object(&context_dest)
	return id
}

//export MakePublicCkksBtpContext
func MakePublicCkksBtpContext(context_handle uint64) uint64 {
	var context_dest CkksBtpContext
	context_src := get_object[CkksBtpContext](context_handle)
	context_dest.parameter = context_src.parameter
	context_dest.btp_parameter = context_src.btp_parameter
	context_dest.sk = nil
	context_dest.pk = context_src.pk
	context_dest.rlk = context_src.rlk
	context_dest.gk = context_src.gk
	context_dest.evk = context_src.evk

	context_dest.encoder = context_src.encoder.ShallowCopy()
	context_dest.encryptor_sk = nil
	if context_src.encryptor_pk != nil {
		context_dest.encryptor_pk = context_src.encryptor_pk.ShallowCopy()
	} else {
		context_dest.encryptor_pk = nil
	}
	context_dest.decryptor = nil
	if context_src.evaluator != nil {
		context_dest.evaluator = context_src.evaluator.ShallowCopy()
	} else {
		context_dest.evaluator = nil
	}
	if context_src.bootstrapper != nil {
		context_dest.bootstrapper = context_src.bootstrapper.ShallowCopy()
	} else {
		context_dest.bootstrapper = nil
	}

	id := insert_object(&context_dest)
	return id
}

//export GetCkksBtpParameter
func GetCkksBtpParameter(context_handle uint64) uint64 {
	context := get_object[CkksBtpContext](context_handle)
	var param BtpParameterSet
	param.SchemeParam = *context.parameter
	param.BtpParam = *context.btp_parameter
	id := insert_object(&param)
	return id
}

//export GetCkksSchemeParameter
func GetCkksSchemeParameter(context_handle uint64) uint64 {
	context := get_object[CkksBtpContext](context_handle)
	param := context.parameter
	id := insert_object(param)
	return id
}

//export CkksBootstrap
func CkksBootstrap(context_handle uint64, x_ciphertext_handle uint64) uint64 {
	context := get_object[CkksBtpContext](context_handle)
	x_ciphertext := get_object[ckks.Ciphertext](x_ciphertext_handle)
	y_ciphertext := context.bootstrapper.Bootstrapp(x_ciphertext)
	id := insert_object(y_ciphertext)
	return id
}

//export ExtractCkksBtpSwkDtS
func ExtractCkksBtpSwkDtS(context_handle uint64) uint64 {
	context := get_object[CkksBtpContext](context_handle)
	id := insert_object(context.evk.SwkDtS)
	return id
}

//export ExtractCkksBtpSwkStD
func ExtractCkksBtpSwkStD(context_handle uint64) uint64 {
	context := get_object[CkksBtpContext](context_handle)
	id := insert_object(context.evk.SwkStD)
	return id
}

//export CreateEmptyCkksBtpContext
func CreateEmptyCkksBtpContext(parameter_handle uint64) uint64 {
	param := get_object[BtpParameterSet](parameter_handle)
	var context CkksBtpContext

	context.parameter = &param.SchemeParam
	context.btp_parameter = &param.BtpParam

	context.encoder = ckks.NewEncoder(*context.parameter)

	context.evk = nil

	context.sk = nil
	context.pk = nil

	context.rlk = nil
	context.gk = nil

	context.bootstrapper = nil
	context.evaluator = nil

	id := insert_object(&context)
	return id
}

//export SetCkksBtpContextRelinKey
func SetCkksBtpContextRelinKey(context_handle uint64, relin_key_handle uint64) {
	context := get_object[CkksBtpContext](context_handle)
	if context.evk == nil {
		context.evk = new(bootstrapping.EvaluationKeys)
	}
	context.evk.Rlk = get_object[rlwe.RelinearizationKey](relin_key_handle)
	context.rlk = context.evk.Rlk
}

//export SetCkksBtpContextGaloisKey
func SetCkksBtpContextGaloisKey(context_handle uint64, galois_key_handle uint64) {
	context := get_object[CkksBtpContext](context_handle)
	if context.evk == nil {
		context.evk = new(bootstrapping.EvaluationKeys)
	}
	context.evk.Rtks = get_object[rlwe.RotationKeySet](galois_key_handle)
	context.gk = context.evk.Rtks
}

//export SetCkksBtpContextSwitchkeyDts
func SetCkksBtpContextSwitchkeyDts(context_handle uint64, switch_key_handle uint64) {
	context := get_object[CkksBtpContext](context_handle)
	if context.evk == nil {
		context.evk = new(bootstrapping.EvaluationKeys)
	}
	context.evk.SwkDtS = get_object[rlwe.SwitchingKey](switch_key_handle)
}

//export SetCkksBtpContextSwitchkeyStd
func SetCkksBtpContextSwitchkeyStd(context_handle uint64, switch_key_handle uint64) {
	context := get_object[CkksBtpContext](context_handle)
	if context.evk == nil {
		context.evk = new(bootstrapping.EvaluationKeys)
	}
	context.evk.SwkStD = get_object[rlwe.SwitchingKey](switch_key_handle)
}

//export CreateCkksBtpContextBootstrapper
func CreateCkksBtpContextBootstrapper(context_handle uint64) {
	context := get_object[CkksBtpContext](context_handle)

	context.evaluator = ckks.NewEvaluator(*context.parameter, rlwe.EvaluationKey{
		Rlk:  context.rlk,
		Rtks: context.gk,
	})

	var err error
	context.bootstrapper, err = bootstrapping.NewBootstrapper(*context.parameter, *context.btp_parameter, *context.evk)
	if err != nil {
		panic(err)
	}
}

//export SerializeCkksBtpContextAdvanced
func SerializeCkksBtpContextAdvanced(context_handle uint64, raw_data **byte, length *C.uint64_t) uint64 {
	context := get_object[CkksBtpContext](context_handle)
	var data_slice []byte
	writer := new(bytes.Buffer)

	param_data, _ := context.parameter.MarshalBinary()
	binary.Write(writer, binary.LittleEndian, uint32(len(param_data)))
	writer.Write(param_data)

	binary.Write(writer, binary.LittleEndian, context.sk != nil)
	if context.sk != nil {
		rlwe.SecretKeyToBytes(context.sk, &context.parameter.Parameters, writer)
	}

	binary.Write(writer, binary.LittleEndian, context.pk != nil)
	if context.pk != nil {
		rlwe.CiphertextQPToBytes(context.pk, &context.parameter.Parameters, writer)
	}

	btp_param_data, _ := context.btp_parameter.MarshalBinary()
	writer.WriteByte(byte(len(btp_param_data)))
	writer.Write(btp_param_data)

	binary.Write(writer, binary.LittleEndian, context.evk.Rlk != nil)
	if context.evk.Rlk != nil {
		rlwe.RelinearizationKeyToByte(context.evk.Rlk, &context.parameter.Parameters, writer)
	}

	binary.Write(writer, binary.LittleEndian, context.evk.Rtks != nil)
	if context.evk.Rtks != nil {
		rlwe.RotationKeySetToBytes(context.evk.Rtks, &context.parameter.Parameters, writer)
	}

	paramsSparse, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN: context.parameter.Parameters.LogN(),
		Q:    context.parameter.Parameters.Q()[:1],
		P:    context.parameter.Parameters.P(),
	})

	binary.Write(writer, binary.LittleEndian, context.evk.SwkDtS != nil)
	if context.evk.SwkDtS != nil {
		rlwe.GadgetCiphertextToBytes(&context.evk.SwkDtS.GadgetCiphertext, &paramsSparse, writer)
	}

	binary.Write(writer, binary.LittleEndian, context.evk.SwkStD != nil)
	if context.evk.SwkStD != nil {
		rlwe.GadgetCiphertextToBytes(&context.evk.SwkStD.GadgetCiphertext, &context.parameter.Parameters, writer)
	}

	data_slice = writer.Bytes()
	*raw_data = (*byte)(unsafe.Pointer(&data_slice[0]))
	*length = (C.uint64_t)(len(data_slice))
	id := insert_object(&data_slice)
	return id
}

//export DeserializeCkksBtpContextAdvanced
func DeserializeCkksBtpContextAdvanced(raw_data *byte, length uint64) uint64 {
	data_slice := unsafe.Slice(raw_data, length)
	var context CkksBtpContext
	var exist bool
	reader := bytes.NewReader(data_slice)

	param := new(ckks.Parameters)
	var param_size uint32
	binary.Read(reader, binary.LittleEndian, &param_size)
	param_data := make([]byte, param_size)
	reader.Read(param_data)
	param.UnmarshalBinary(param_data)
	context.parameter = param

	binary.Read(reader, binary.LittleEndian, &exist)
	if exist {
		sk := rlwe.BytesToSecretKey(reader)
		context.sk = &sk
	} else {
		context.sk = nil
	}

	binary.Read(reader, binary.LittleEndian, &exist)
	if exist {
		pk := rlwe.BytesToCiphertextQP(reader)
		context.pk = &pk
		context.pk.Decompress(&param.Parameters)
	} else {
		context.pk = nil
	}

	btp_param := new(bootstrapping.Parameters)
	btp_param_size, _ := reader.ReadByte()
	btp_param_data := make([]byte, btp_param_size)
	reader.Read(btp_param_data)
	btp_param.UnmarshalBinary(btp_param_data)
	context.btp_parameter = btp_param

	context.evk = new(bootstrapping.EvaluationKeys)

	binary.Read(reader, binary.LittleEndian, &exist)
	if exist {
		btp_rlk := rlwe.BytesToRelinearizationKey(reader)
		context.evk.Rlk = &btp_rlk

		for _, swk := range context.evk.Rlk.Keys {
			swk.Decompress(&context.parameter.Parameters)
		}
	} else {
		context.evk.Rlk = nil
	}

	binary.Read(reader, binary.LittleEndian, &exist)
	if exist {
		btp_gk := rlwe.BytesToRotationKeySet(reader)
		context.evk.Rtks = &btp_gk

		for _, swk := range context.evk.Rtks.Keys {
			swk.Decompress(&context.parameter.Parameters)
		}
	} else {
		context.evk.Rtks = nil
	}

	binary.Read(reader, binary.LittleEndian, &exist)
	if exist {
		context.evk.SwkDtS = new(rlwe.SwitchingKey)
		context.evk.SwkDtS.GadgetCiphertext = rlwe.BytesToGadgetCiphertext(reader)

		context.evk.SwkDtS.Decompress(&context.parameter.Parameters)
	} else {
		context.evk.SwkDtS = nil
	}

	binary.Read(reader, binary.LittleEndian, &exist)
	if exist {
		context.evk.SwkStD = new(rlwe.SwitchingKey)
		context.evk.SwkStD.GadgetCiphertext = rlwe.BytesToGadgetCiphertext(reader)

		context.evk.SwkStD.Decompress(&context.parameter.Parameters)
	} else {
		context.evk.SwkStD = nil
	}

	context.rlk = context.evk.Rlk
	context.gk = context.evk.Rtks

	init_ckks_context(&context.CkksContext)
	context.bootstrapper, _ = bootstrapping.NewBootstrapper(*context.parameter, *context.btp_parameter, *context.evk)

	id := insert_object(&context)
	return id
}
