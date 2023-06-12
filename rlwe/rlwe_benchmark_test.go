package rlwe

import (
	"bufio"
	"bytes"
	"encoding/json"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
)

func BenchmarkRLWE(b *testing.B) {

	var err error

	defaultParamsLiteral := TestParamsLiteral[:1]

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParamsLiteral = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, paramsLit := range defaultParamsLiteral {

		var params Parameters
		if params, err = NewParametersFromLiteral(paramsLit); err != nil {
			b.Fatal(err)
		}

		tc := NewTestContext(params)

		for _, testSet := range []func(tc *TestContext, b *testing.B){
			benchKeyGenerator,
			benchEncryptor,
			benchDecryptor,
			benchEvaluator,
			benchMarshalling,
		} {
			testSet(tc, b)
			runtime.GC()
		}
	}
}

func benchKeyGenerator(tc *TestContext, b *testing.B) {

	params := tc.params
	kgen := tc.kgen

	b.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenSecretKey"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenSecretKey(tc.sk)
		}
	})

	b.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenPublicKey"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenPublicKey(tc.sk, tc.pk)
		}

	})

	b.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenEvaluationKey"), func(b *testing.B) {
		sk0, sk1 := tc.sk, kgen.GenSecretKeyNew()
		evk := NewEvaluationKey(params, params.MaxLevelQ(), params.MaxLevelP())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			kgen.GenEvaluationKey(sk0, sk1, evk)
		}
	})
}

func benchEncryptor(tc *TestContext, b *testing.B) {

	params := tc.params

	b.Run(testString(params, params.MaxLevel(), "Encryptor/EncryptZero/SecretKey"), func(b *testing.B) {
		ct := NewCiphertext(params, 1, params.MaxLevel())
		enc := tc.enc.WithKey(tc.sk)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			enc.EncryptZero(ct)
		}

	})

	b.Run(testString(params, params.MaxLevel(), "Encryptor/EncryptZero/PublicKey"), func(b *testing.B) {
		ct := NewCiphertext(params, 1, params.MaxLevel())
		enc := tc.enc.WithKey(tc.pk)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			enc.EncryptZero(ct)
		}
	})
}

func benchDecryptor(tc *TestContext, b *testing.B) {

	params := tc.params

	b.Run(testString(params, params.MaxLevel(), "Decryptor/Decrypt"), func(b *testing.B) {
		dec := tc.dec
		ct := tc.enc.EncryptZeroNew(params.MaxLevel())
		pt := NewPlaintext(params, ct.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			dec.Decrypt(ct, pt)
		}
	})
}

func benchEvaluator(tc *TestContext, b *testing.B) {

	params := tc.params
	kgen := tc.kgen
	sk := tc.sk
	eval := tc.eval

	b.Run(testString(params, params.MaxLevel(), "Evaluator/GadgetProduct"), func(b *testing.B) {

		ct := NewEncryptor(params, sk).EncryptZeroNew(params.MaxLevel())
		evk := kgen.GenEvaluationKeyNew(sk, kgen.GenSecretKeyNew())

		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			eval.GadgetProduct(ct.Level(), &ct.Value[1], &evk.GadgetCiphertext, ct)
		}
	})
}

func benchMarshalling(tc *TestContext, b *testing.B) {
	params := tc.params
	sk := tc.sk

	ctf := NewEncryptor(params, sk).EncryptZeroNew(params.MaxLevel())
	ct := ctf.Value

	badbuf := bytes.NewBuffer(make([]byte, ct.BinarySize()))
	b.Run(testString(params, params.MaxLevel(), "Marshalling/WriteToBadBuf"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			_, err := ct.WriteTo(badbuf)

			b.StopTimer()
			if err != nil {
				b.Fatal(err)
			}
			badbuf.Reset()
			b.StartTimer()
		}
	})

	runtime.GC()

	bytebuff := bytes.NewBuffer(make([]byte, ct.BinarySize()))
	bufiobuf := bufio.NewWriter(bytebuff)
	b.Run(testString(params, params.MaxLevel(), "Marshalling/WriteToIOBuf"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			_, err := ct.WriteTo(bufiobuf)

			b.StopTimer()
			if err != nil {
				b.Fatal(err)
			}
			bytebuff.Reset()
			bufiobuf.Reset(bytebuff)
			b.StartTimer()
		}
	})

	runtime.GC()

	bsliceour := make([]byte, ct.BinarySize())
	ourbuf := buffer.NewBuffer(bsliceour)
	b.Run(testString(params, params.MaxLevel(), "Marshalling/WriteToOurBuf"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			_, err := ct.WriteTo(ourbuf)

			b.StopTimer()
			if err != nil {
				b.Fatal(err)
			}
			ourbuf.Reset()
			b.StartTimer()
		}
	})

	runtime.GC()
	require.Equal(b, ct.BinarySize(), len(ourbuf.Bytes()))

	rdr := bytes.NewReader(ourbuf.Bytes())
	//bufiordr := bufio.NewReaderSize(rdr, len(ourbuf.Bytes()))
	bufiordr := bufio.NewReader(rdr)
	ct2f := NewCiphertext(tc.params, 1, tc.params.MaxLevel())
	ct2 := ct2f.Value
	b.Run(testString(params, params.MaxLevel(), "Marshalling/ReadFromIO"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {

			_, err := ct2.ReadFrom(bufiordr)

			b.StopTimer()
			if err != nil {
				b.Fatal(err)
			}
			rdr.Seek(0, 0)
			bufiordr.Reset(rdr)
			b.StartTimer()
		}
	})

	// require.True(b, ct.Equal(ct2))

	ct3f := NewCiphertext(tc.params, 1, tc.params.MaxLevel())
	ct3 := ct3f.Value
	b.Run(testString(params, params.MaxLevel(), "Marshalling/ReadFromOur"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			_, err := ct3.ReadFrom(ourbuf)

			b.StopTimer()
			if err != nil {
				b.Fatal(err)
			}
			ourbuf.Reset()
			b.StartTimer()
		}
	})
	require.True(b, ct.Equal(ct3))
}
