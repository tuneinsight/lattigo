package rlwe

import "fmt"

// EvaluationKeySetInterface is an interface implementing methods
// to load the RelinearizationKey and GaloisKeys in the Evaluator.
// This interface must support concurrent calls on the methods
// GetGaloisKey and GetRelinearizationKey.
type EvaluationKeySetInterface interface {

	// GetGaloisKey retrieves the Galois key for the automorphism X^{i} -> X^{i*galEl}.
	GetGaloisKey(galEl uint64) (evk *GaloisKey, err error)

	// GetGaloisKeysList returns the list of all the Galois elements
	// for which a Galois key exists in the object.
	GetGaloisKeysList() (galEls []uint64)

	// GetRelinearizationKey retrieves the RelinearizationKey.
	GetRelinearizationKey() (evk *RelinearizationKey, err error)
}

// EvaluationKeySet is a generic struct that complies to the EvaluationKeySetInterface interface.
// This interface can be re-implemented by users to suit application specific requirement.
type EvaluationKeySet struct {
	*RelinearizationKey
	GaloisKeys map[uint64]*GaloisKey
}

// NewEvaluationKeySet returns a new EvaluationKeySet with nil RelinearizationKey and empty GaloisKeys map.
func NewEvaluationKeySet() (evk *EvaluationKeySet) {
	return &EvaluationKeySet{
		RelinearizationKey: nil,
		GaloisKeys:         make(map[uint64]*GaloisKey),
	}
}

// GetGaloisKey retrieves the Galois key for the automorphism X^{i} -> X^{i*galEl}.
func (evk *EvaluationKeySet) GetGaloisKey(galEl uint64) (gk *GaloisKey, err error) {
	var ok bool
	if gk, ok = evk.GaloisKeys[galEl]; !ok {
		return nil, fmt.Errorf("GaloiKey[%d] is nil", galEl)
	}

	return
}

// GetGaloisKeysList returns the list of all the Galois elements
// for which a Galois key exists in the object.
func (evk *EvaluationKeySet) GetGaloisKeysList() (galEls []uint64) {

	if evk.GaloisKeys == nil {
		return []uint64{}
	}

	galEls = make([]uint64, len(evk.GaloisKeys))

	var i int
	for galEl := range evk.GaloisKeys {
		galEls[i] = galEl
		i++
	}

	return
}

// GetRelinearizationKey retrieves the RelinearizationKey.
func (evk *EvaluationKeySet) GetRelinearizationKey() (rk *RelinearizationKey, err error) {
	if evk.RelinearizationKey != nil {
		return evk.RelinearizationKey, nil
	}

	return nil, fmt.Errorf("RelinearizationKey is nil")
}
