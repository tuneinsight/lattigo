package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// RingPackingEvaluationKey is a struct storing the
// ring packing evaluation keys.
// All fields of this struct are public, enabling
// custom instantiations.
type RingPackingEvaluationKey struct {
	// Parameters are the different Parameters among
	// which a ciphertext will be switched during the
	// procedure. These parameters share the same primes
	// but support different ring degrees.
	Parameters map[int]ParameterProvider

	// RingSwitchingKeys are the ring degree switching keys
	// indexed as map[inputLogN][outputLogN]
	RingSwitchingKeys map[int]map[int]*EvaluationKey

	// RepackKeys are the [EvaluationKey] used for the
	// RLWE repacking.
	RepackKeys map[int]EvaluationKeySet

	// ExtractKeys are the [EvaluationKey] used for the
	// RLWE extraction.
	ExtractKeys map[int]EvaluationKeySet
}

// MinLogN returns the minimum Log(N) among the supported ring degrees.
// This method requires that the field Parameters of [RingPackingEvaluationKey]
// has been populated.
func (rpk RingPackingEvaluationKey) MinLogN() (minLogN int) {
	return utils.GetSortedKeys(rpk.Parameters)[0]
}

// MaxLogN returns the maximum Log(N) among the supported ring degrees.
// This method requires that the field Parameters of [RingPackingEvaluationKey]
// has been populated.
func (rpk RingPackingEvaluationKey) MaxLogN() (maxLogN int) {
	return utils.GetSortedKeys(rpk.Parameters)[len(rpk.Parameters)-1]
}

// GenRingSwitchingKeys generates the [Parameter]s and [EvaluationKey]s
// to be able to split an [Ciphertext] into two [Ciphertext]s of half
// the ring degree and merge two [Ciphertext]s into one [Ciphertext]
// of twice the ring degree.
//
// The method returns the [Parameter]s, [EvaluationKey]s and ephemeral
// [SecretKey]s used to generate the ring-switching [EvaluationKey]s.
//
// See the methods [RingPackingEvaluator.Split] and [RingPackingEvaluator.Repack].
//
// This function will return an error if minLogN >= params.LogN().
func (rpk *RingPackingEvaluationKey) GenRingSwitchingKeys(params ParameterProvider, sk *SecretKey, minLogN int, evkParams EvaluationKeyParameters) (ski map[int]*SecretKey, err error) {

	p := *params.GetRLWEParameters()

	if minLogN >= p.LogN() {
		return nil, fmt.Errorf("invalid minLogN: cannot be equal or larger than params.LogN()")
	}

	LevelQ, LevelP, _, _ := ResolveEvaluationKeyParameters(p, []EvaluationKeyParameters{evkParams})

	Q := p.Q()
	P := p.P()

	parameters := map[int]ParameterProvider{}
	parameters[p.LogN()] = &p

	ski = map[int]*SecretKey{}
	ski[p.LogN()] = sk

	kgen := map[int]*KeyGenerator{}
	kgen[p.LogN()] = NewKeyGenerator(p)

	for i := minLogN; i < p.LogN(); i++ {

		var pi Parameters
		if pi, err = NewParametersFromLiteral(ParametersLiteral{
			LogN:         i,
			Q:            Q[:LevelQ+1],
			P:            P[:LevelP+1],
			NTTFlag:      p.NTTFlag(),
			DefaultScale: p.DefaultScale(),
		}); err != nil {
			return nil, fmt.Errorf("NewParametersFromLiteral: %w", err)
		}

		kgen[i] = NewKeyGenerator(pi)
		ski[i] = kgen[i].GenSecretKeyNew()
		parameters[i] = &pi
	}

	// Ring switching evaluation keys
	RingSwitchingKeys := map[int]map[int]*EvaluationKey{}

	for i := minLogN; i < p.LogN()+1; i++ {
		RingSwitchingKeys[i] = map[int]*EvaluationKey{}
	}

	for i := minLogN; i < p.LogN(); i++ {
		RingSwitchingKeys[i][i+1] = kgen[i+1].GenEvaluationKeyNew(ski[i], ski[i+1], evkParams)
		RingSwitchingKeys[i+1][i] = kgen[i+1].GenEvaluationKeyNew(ski[i+1], ski[i], evkParams)
	}

	rpk.Parameters = parameters
	rpk.RingSwitchingKeys = RingSwitchingKeys

	return ski, nil
}

// GenRepackEvaluationKeys generates the set of params.LogN() [EvaluationKey]s necessary to perform the repacking operation.
// See [RingPackingEvaluator.Repack] for additional information.
func (rpk *RingPackingEvaluationKey) GenRepackEvaluationKeys(params ParameterProvider, sk *SecretKey, evkParams EvaluationKeyParameters) {
	p := *params.GetRLWEParameters()

	if rpk.RepackKeys == nil {
		rpk.RepackKeys = map[int]EvaluationKeySet{}
	}

	rpk.RepackKeys[p.LogN()] = NewMemEvaluationKeySet(nil, NewKeyGenerator(p).GenGaloisKeysNew(GaloisElementsForPack(p, p.LogN()), sk, evkParams)...)
}

// GenExtractEvaluationKeys generates the set of params.LogN() [EvaluationKey]s necessary to perform the extraction operation.
// See [RingPackingEvaluator.Extract] for additional information.
func (rpk *RingPackingEvaluationKey) GenExtractEvaluationKeys(params ParameterProvider, sk *SecretKey, evkParams EvaluationKeyParameters) {
	p := *params.GetRLWEParameters()

	if rpk.ExtractKeys == nil {
		rpk.ExtractKeys = map[int]EvaluationKeySet{}
	}

	rpk.ExtractKeys[p.LogN()] = NewMemEvaluationKeySet(nil, NewKeyGenerator(p).GenGaloisKeysNew(GaloisElementsForExpand(p, p.LogN()), sk, evkParams)...)
}

// GaloisElementsForExpand returns the list of Galois elements required
// to perform the `Expand` operation with parameter `logN`.
func GaloisElementsForExpand(params ParameterProvider, logN int) (galEls []uint64) {
	galEls = make([]uint64, logN)

	NthRoot := params.GetRLWEParameters().RingQ().NthRoot()

	for i := 0; i < logN; i++ {
		galEls[i] = uint64(NthRoot/(2<<i) + 1)
	}

	return
}

// GaloisElementsForPack returns the list of Galois elements required to perform the `Pack` operation.
func GaloisElementsForPack(params ParameterProvider, logGap int) (galEls []uint64) {

	p := params.GetRLWEParameters()

	// Sanity check
	if logGap > p.LogN() || logGap < 0 {
		panic(fmt.Errorf("cannot GaloisElementsForPack: logGap > logN || logGap < 0"))
	}

	galEls = make([]uint64, 0, logGap)
	for i := 0; i < logGap; i++ {
		galEls = append(galEls, p.GaloisElement(1<<i))
	}

	switch p.RingType() {
	case ring.Standard:
		if logGap == p.LogN() {
			galEls = append(galEls, p.GaloisElementOrderTwoOrthogonalSubgroup())
		}
	default:
		panic("cannot GaloisElementsForPack: invalid ring type")
	}

	return
}
