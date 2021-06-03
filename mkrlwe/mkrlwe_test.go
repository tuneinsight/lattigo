package mkrlwe

import (
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/rlwe"
)

func Test_MKRLWE(t *testing.T) {

	testMarshaler(t)

}

func testMarshaler(t *testing.T) {

	t.Run("Test marshaler MKDecomposedPoly", func(t *testing.T) {

		params, err1 := bfv.NewParametersFromLiteral(bfv.PN12QP101pq)
		if err1 != nil {
			t.Fail()
		}

		ringQ := GetRingQ(&params.Parameters)

		decPol := NewDecomposedPoly(ringQ, 3)
		decPol.Poly[0].Coeffs[0][0] = 23
		bytes, err2 := decPol.MarshalBinary()

		if err2 != nil {
			t.Fail()
		}
		resPol := new(MKDecomposedPoly)

		resPol.UnmarshalBinary(bytes)

		if resPol.Poly[0].Coeffs[0][0] != 23 {
			t.Error("Marshaling and Unmarshaling decomposed polynomial fails")
		}
	})

	t.Run("Test marshaler MKEvaluationKey", func(t *testing.T) {

		params, err1 := bfv.NewParametersFromLiteral(bfv.PN12QP101pq)
		if err1 != nil {
			t.Fail()
		}

		ringQ := GetRingQ(&params.Parameters)

		evalKey := NewMKEvaluationKey(ringQ, 3, &params.Parameters)
		evalKey.Key01.Value[0][0].Coeffs[0][0] = 23
		bytes, err2 := evalKey.MarshalBinary()

		if err2 != nil {
			t.Fail()
		}

		resKey := new(MKEvaluationKey)
		resKey.UnmarshalBinary(bytes)
		if resKey.PeerID != 3 || resKey.Key01.Value[0][0].Coeffs[0][0] != 23 {
			t.Error("Marshaling and Unmarshaling evaluation key fails")
		}

	})

	t.Run("Test marshaler MKEvalGalKey", func(t *testing.T) {

		params, err1 := bfv.NewParametersFromLiteral(bfv.PN12QP101pq)
		if err1 != nil {
			t.Fail()
		}

		evalKey := &MKEvalGalKey{Key: rlwe.NewSwitchingKey(params.Parameters), PeerID: 3}

		evalKey.Key.Value[0][0].Coeffs[0][0] = 23

		bytes, err2 := evalKey.MarshalBinary()

		if err2 != nil {
			t.Fail()
		}

		resKey := new(MKEvalGalKey)

		resKey.UnmarshalBinary(bytes)

		if resKey.PeerID != 3 || resKey.Key.Value[0][0].Coeffs[0][0] != 23 {
			t.Error("Marshaling and Unmarshaling galois evaluation key fails")
		}

	})

	t.Run("Test marshaler MKPublicKey", func(t *testing.T) {

		params, err1 := bfv.NewParametersFromLiteral(bfv.PN12QP101pq)
		if err1 != nil {
			t.Fail()
		}

		ringQ := GetRingQ(&params.Parameters)

		pubKey := &MKPublicKey{[2]*MKDecomposedPoly{NewDecomposedPoly(ringQ, 2), NewDecomposedPoly(ringQ, 2)}, 3}

		pubKey.Key[0].Poly[0].Coeffs[0][0] = 23

		bytes, err2 := pubKey.MarshalBinary()

		if err2 != nil {
			t.Fail()
		}

		resKey := new(MKPublicKey)

		resKey.UnmarshalBinary(bytes)

		if resKey.PeerID != 3 || resKey.Key[0].Poly[0].Coeffs[0][0] != 23 {
			t.Error("Marshaling and Unmarshaling public key fails")
		}

	})
}
