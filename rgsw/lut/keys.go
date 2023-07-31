package lut

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/rgsw"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

const (
	// Parameter w of Algorithm 3 in https://eprint.iacr.org/2022/198
	windowSize = 10
)

// BlindRotatationEvaluationKeySet is a interface implementing methods
// to load the blind rotation keys (RGSW) and automorphism keys
// (via the rlwe.EvaluationKeySet interface).
// Implementation of this interface must be safe for concurrent use.
type BlindRotatationEvaluationKeySet interface {

	// GetBlindRotationKey should return RGSW(X^{s[i]})
	GetBlindRotationKey(i int) (brk *rgsw.Ciphertext, err error)

	// GetEvaluationKeySet should return an rlwe.EvaluationKeySet
	// providing access to all the required automorphism keys.
	GetEvaluationKeySet() (evk rlwe.EvaluationKeySet, err error)
}

// MemBlindRotatationEvaluationKeySet is a basic in-memory implementation of the BlindRotatationEvaluationKeySet interface.
type MemBlindRotatationEvaluationKeySet struct {
	BlindRotationKeys []*rgsw.Ciphertext
	AutomorphismKeys  []*rlwe.GaloisKey
}

func (evk MemBlindRotatationEvaluationKeySet) GetBlindRotationKey(i int) (*rgsw.Ciphertext, error) {
	return evk.BlindRotationKeys[i], nil
}

func (evk MemBlindRotatationEvaluationKeySet) GetEvaluationKeySet() (rlwe.EvaluationKeySet, error) {
	return rlwe.NewMemEvaluationKeySet(nil, evk.AutomorphismKeys...), nil
}

// GenEvaluationKeyNew generates a new LUT evaluation key
func GenEvaluationKeyNew(paramsRLWE rlwe.Parameters, skRLWE *rlwe.SecretKey, paramsLWE rlwe.Parameters, skLWE *rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (key MemBlindRotatationEvaluationKeySet, err error) {

	skLWECopy := skLWE.CopyNew()
	paramsLWE.RingQ().AtLevel(0).INTT(skLWECopy.Value.Q, skLWECopy.Value.Q)
	paramsLWE.RingQ().AtLevel(0).IMForm(skLWECopy.Value.Q, skLWECopy.Value.Q)
	sk := make([]*big.Int, paramsLWE.N())
	for i := range sk {
		sk[i] = new(big.Int)
	}
	paramsLWE.RingQ().AtLevel(0).PolyToBigintCentered(skLWECopy.Value.Q, 1, sk)

	encryptor, err := rgsw.NewEncryptor(paramsRLWE, skRLWE)
	if err != nil {
		return key, err
	}

	levelQ, levelP, BaseTwoDecomposition := rlwe.ResolveEvaluationKeyParameters(paramsRLWE, evkParams)

	skiRGSW := make([]*rgsw.Ciphertext, paramsLWE.N())

	ptXi := make(map[int]*rlwe.Plaintext)

	for i, si := range sk {

		siInt := int(si.Int64())

		if _, ok := ptXi[siInt]; !ok {

			pt := &rlwe.Plaintext{}
			pt.MetaData = &rlwe.MetaData{}
			pt.IsNTT = true
			pt.Value = paramsRLWE.RingQ().NewMonomialXi(siInt)
			paramsRLWE.RingQ().NTT(pt.Value, pt.Value)

			ptXi[siInt] = pt
		}

		skiRGSW[i] = rgsw.NewCiphertext(paramsRLWE, levelQ, levelP, BaseTwoDecomposition)

		if err = encryptor.Encrypt(ptXi[siInt], skiRGSW[i]); err != nil {
			return
		}
	}

	kgen := rlwe.NewKeyGenerator(paramsRLWE)

	galEls := make([]uint64, windowSize)
	for i := 0; i < windowSize; i++ {
		galEls[i] = paramsRLWE.GaloisElement(i + 1)
	}

	galEls = append(galEls, paramsRLWE.RingQ().NthRoot()-ring.GaloisGen)

	gks, err := kgen.GenGaloisKeysNew(galEls, skRLWE, rlwe.EvaluationKeyParameters{
		LevelQ:               utils.Pointy(levelQ),
		LevelP:               utils.Pointy(levelP),
		BaseTwoDecomposition: utils.Pointy(BaseTwoDecomposition),
	})

	if err != nil {
		return MemBlindRotatationEvaluationKeySet{}, err
	}

	return MemBlindRotatationEvaluationKeySet{BlindRotationKeys: skiRGSW, AutomorphismKeys: gks}, nil
}
