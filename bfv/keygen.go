package bfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeyGenerator is an interface implementing the methods of the keyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretkeyWithDistrib(p float64) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenRelinKey(sk *SecretKey, maxDegree uint64) (evk *EvaluationKey)
	GenSwitchingKey(skIn, skOut *SecretKey) (evk *SwitchingKey)
	GenRot(rotType RotationType, sk *SecretKey, k uint64, rotKey *RotationKeys)
	GenRotationKeysPow2(sk *SecretKey) (rotKey *RotationKeys)
}

// keyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keyGenerator struct {
	params          *Parameters
	ringQP          *ring.Ring
	pBigInt         *big.Int
	polypool        [2]*ring.Poly
	gaussianSampler *ring.GaussianSampler
	uniformSampler  *ring.UniformSampler
}

// SecretKey is a structure that stores the SecretKey.
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the PublicKey.
type PublicKey struct {
	pk [2]*ring.Poly
}

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	keys map[uint64]*SwitchingKey
}

func (rtk RotationKeys) Delete() {
	for k := range rtk.keys {
		delete(rtk.keys, k)
	}
}

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey []*SwitchingKey
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	evakey [][2]*ring.Poly
}

// Get returns the switching key backing slice.
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.evakey
}

func (swk *SwitchingKey) Copy(other *SwitchingKey) {
	if other == nil {
		return
	}
	if len(swk.evakey) == 0 {
		swk.evakey = make([][2]*ring.Poly, len(other.evakey), len(other.evakey))
		for i, o := range other.evakey {
			n, q := uint64(o[0].GetDegree()), uint64(o[0].GetLenModuli())
			swk.evakey[i] = [2]*ring.Poly{ring.NewPoly(n, q), ring.NewPoly(n, q)}
		}
	}
	for i, o := range other.evakey {
		swk.evakey[i][0].Copy(o[0])
		swk.evakey[i][1].Copy(o[1])
	}
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params *Parameters) KeyGenerator {

	var ringQP *ring.Ring
	var err error
	if ringQP, err = ring.NewRing(params.N(), append(params.qi, params.pi...)); err != nil {
		panic(err)
	}

	var pBigInt *big.Int
	if len(params.pi) != 0 {
		pBigInt = ring.NewUint(1)
		for _, pi := range params.pi {
			pBigInt.Mul(pBigInt, ring.NewUint(pi))
		}
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &keyGenerator{
		params:          params.Copy(),
		ringQP:          ringQP,
		pBigInt:         pBigInt,
		polypool:        [2]*ring.Poly{ringQP.NewPoly(), ringQP.NewPoly()},
		gaussianSampler: ring.NewGaussianSampler(prng, ringQP, params.Sigma(), uint64(6*params.Sigma())),
		uniformSampler:  ring.NewUniformSampler(prng, ringQP),
	}
}

// GenSecretKey creates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretkeyWithDistrib(1.0 / 3)
}

// GenSecretkeyWithDistrib creates a new SecretKey with the distribution [(1-p)/2, p, (1-p)/2].
func (keygen *keyGenerator) GenSecretkeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.ringQP, p, true)

	sk = new(SecretKey)
	sk.sk = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.sk, sk.sk)
	return sk
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params *Parameters) *SecretKey {

	sk := new(SecretKey)
	sk.sk = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	return sk
}

// Get returns the polynomial of the target SecretKey.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the polynomial of the target secret key as the input polynomial.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// GenPublicKey generates a new PublicKey from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringQP := keygen.ringQP

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]

	pk.pk[0] = keygen.gaussianSampler.ReadNew()
	ringQP.NTT(pk.pk[0], pk.pk[0])
	pk.pk[1] = keygen.uniformSampler.ReadNew()

	ringQP.MulCoeffsMontgomeryAndSub(sk.sk, pk.pk[1], pk.pk[0])

	return pk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {

	pk = new(PublicKey)

	pk.pk[0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	pk.pk[1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))

	return
}

// Get returns the polynomials of the PublicKey.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the polynomial of the PublicKey as the input polynomials.
func (pk *PublicKey) Set(p [2]*ring.Poly) {
	pk.pk[0] = p[0].CopyNew()
	pk.pk[1] = p[1].CopyNew()
}

// NewKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding PublicKey.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// NewRelinKey generates a new evaluation key from the provided SecretKey. It will be used to relinearize a ciphertext (encrypted under a PublicKey generated from the provided SecretKey)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize.
func (keygen *keyGenerator) GenRelinKey(sk *SecretKey, maxDegree uint64) (evk *EvaluationKey) {

	if keygen.ringQP == nil {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	evk = new(EvaluationKey)
	evk.evakey = make([]*SwitchingKey, maxDegree)
	for i := range evk.evakey {
		evk.evakey[i] = NewSwitchingKey(keygen.params)
	}

	keygen.polypool[0].Copy(sk.Get()) // TODO Remove ?

	ringQP := keygen.ringQP

	keygen.polypool[1].Copy(sk.Get())
	for i := uint64(0); i < maxDegree; i++ {
		ringQP.MulCoeffsMontgomery(keygen.polypool[1], sk.Get(), keygen.polypool[1])
		keygen.newSwitchingKey(keygen.polypool[1], sk.Get(), evk.evakey[i])
	}

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

// NewRelinKey creates a new EvaluationKey with zero values.
func NewRelinKey(params *Parameters, maxDegree uint64) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)

	beta := params.Beta()

	evakey.evakey = make([]*SwitchingKey, maxDegree)

	for w := uint64(0); w < maxDegree; w++ {

		evakey.evakey[w] = new(SwitchingKey)

		evakey.evakey[w].evakey = make([][2]*ring.Poly, beta)

		for i := uint64(0); i < beta; i++ {

			evakey.evakey[w].evakey[i][0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
			evakey.evakey[w].evakey[i][1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
		}
	}

	return
}

// Get returns the slice of SwitchingKeys of the target EvaluationKey.
func (evk *EvaluationKey) Get() []*SwitchingKey {
	return evk.evakey
}

// Set sets the polynomial of the target EvaluationKey as the input polynomials.
func (evk *EvaluationKey) Set(rlk [][][2]*ring.Poly) {

	evk.evakey = make([]*SwitchingKey, len(rlk))
	for i := range rlk {
		evk.evakey[i] = new(SwitchingKey)
		evk.evakey[i].evakey = make([][2]*ring.Poly, len(rlk[i]))
		for j := range rlk[i] {
			evk.evakey[i].evakey[j][0] = rlk[i][j][0].CopyNew()
			evk.evakey[i].evakey[j][1] = rlk[i][j][1].CopyNew()
		}
	}
}

// GenSwitchingKey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (swkOut *SwitchingKey) {

	if keygen.ringQP == nil {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	swkOut = NewSwitchingKey(keygen.params)

	keygen.ringQP.Copy(skInput.Get(), keygen.polypool[0]) // TODO: remove and pass skInput directly ?
	keygen.newSwitchingKey(keygen.polypool[0], skOutput.Get(), swkOut)
	keygen.polypool[0].Zero()
	return
}

// NewSwitchingKey returns a new SwitchingKey with zero values.
func NewSwitchingKey(params *Parameters) (evakey *SwitchingKey) {

	evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey = make([][2]*ring.Poly, params.Beta())

	for i := uint64(0); i < params.Beta(); i++ {
		evakey.evakey[i][0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
		evakey.evakey[i][1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	}

	return
}

// NewRotationKeys returns a new empty RotationKeys struct.
func NewRotationKeys(params *Parameters) (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	rotKey.keys = make(map[uint64]*SwitchingKey, 0)
	return
}

// GenRot populates the target RotationKeys with a SwitchingKey for the desired rotation type and amount.
func (keygen *keyGenerator) GenRot(rotType RotationType, sk *SecretKey, k uint64, rotKeys *RotationKeys) {

	if keygen.ringQP == nil {
		panic("cannot generate rotation keys when modulus P is empty")
	}

	galEl := getGaloisElementForRotation(rotType, k, keygen.ringQP.N)
	galElRev := getGaloisElementForRotationRev(rotType, k, keygen.ringQP.N)

	rotKey, inSet := rotKeys.keys[galEl]
	if !inSet {
		rotKey = NewSwitchingKey(keygen.params)
		rotKeys.keys[galEl] = rotKey
	}

	keygen.genrotKey(sk.sk, galElRev, rotKey)
}

// GenRotationKeysPow2 generates a new rotation key with all the power-of-two rotations to the left and right, as well as the conjugation.
func (keygen *keyGenerator) GenRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys) {

	if keygen.ringQP == nil {
		panic("Cannot GenRotationKeysPow2: modulus P is empty")
	}

	rotKey = NewRotationKeys(keygen.params)

	for n := uint64(1); n < 1<<(keygen.params.LogN()-1); n <<= 1 {
		keygen.GenRot(RotationLeft, skOutput, n, rotKey)
		keygen.GenRot(RotationRight, skOutput, n, rotKey)
	}

	keygen.GenRot(RotationRow, skOutput, 0, rotKey)

	return
}

func (rotKeys *RotationKeys) SetRotKeyGalEl(galoisEl uint64, swk *SwitchingKey) {
	rotKey, inSet := rotKeys.keys[galoisEl]
	if !inSet {
		rotKey = new(SwitchingKey)
		rotKeys.keys[galoisEl] = rotKey
	}
	if rotKey != swk {
		rotKey.Copy(swk)
	}
}

// SetRotKey sets the target RotationKeys' SwitchingKey for the specified rotation type and amount with the input polynomials.
func (rotKeys *RotationKeys) SetRotKey(rotType RotationType, k uint64, swk *SwitchingKey) {

	galEl := getGaloisElementForRotation(rotType, k, uint64(swk.evakey[0][0].GetDegree()))

	rotKeys.SetRotKeyGalEl(galEl, swk)
}

func (rotKeys *RotationKeys) GetRotKey(galoisEl uint64) (*SwitchingKey, bool) {
	rotKey, inSet := rotKeys.keys[galoisEl]
	return rotKey, inSet
}

func (keygen *keyGenerator) genrotKey(sk *ring.Poly, gen uint64, swkOut *SwitchingKey) {

	skIn := sk
	skOut := keygen.polypool[1]

	ring.PermuteNTT(skIn, gen, skOut)

	keygen.newSwitchingKey(skIn, skOut, swkOut)

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly, swkOut *SwitchingKey) {

	ringQP := keygen.ringQP

	alpha := keygen.params.Alpha()
	beta := keygen.params.Beta()

	var index uint64

	// delta_sk = skIn - skOut = GaloisEnd(skOut, rotation) - skOut

	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	for i := uint64(0); i < beta; i++ {

		// e
		keygen.gaussianSampler.Read(swkOut.evakey[i][0])
		ringQP.NTTLazy(swkOut.evakey[i][0], swkOut.evakey[i][0])
		ringQP.MForm(swkOut.evakey[i][0], swkOut.evakey[i][0])
		// a
		keygen.uniformSampler.Read(swkOut.evakey[i][1])

		// e + skIn * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0

		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := keygen.polypool[0].Coeffs[index]
			p1tmp := swkOut.evakey[i][0].Coeffs[index]

			for w := uint64(0); w < ringQP.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// Handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QiCount() {
				break
			}

		}

		// skIn * (qiBarre*qiStar) * 2^w - a*sk + e
		ringQP.MulCoeffsMontgomeryAndSub(swkOut.evakey[i][1], skOut, swkOut.evakey[i][0])
	}

	return
}
