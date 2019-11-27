package ckks

import (
	"math"

	"github.com/ldsec/lattigo/ring"
)

type keyGeneratorContext struct {
	// Context parameters
	logN uint64
	n    uint64

	// Number of available levels
	levels uint64

	// Contexts
	specialPrimes []uint64
	alpha         uint64
	beta          uint64
	contextKeys   *ring.Context

	// Samplers
	gaussianSampler *ring.KYSampler

	// Rotation params
	galElRotRow      uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64
}

func newKeyGeneratorContext(params *Parameters) *keyGeneratorContext {
	logN := uint64(params.LogN)
	n := uint64(1 << uint64(params.LogN))
	maxSlots := uint64(1 << (uint64(params.LogN) - 1))
	levels := uint64(len(params.Modulichain))

	alpha := uint64(len(params.P))
	beta := uint64(math.Ceil(float64(levels) / float64(alpha)))

	scalechain := make([]float64, len(params.Modulichain))

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[uint64]uint64)
	for i, qi := range params.Modulichain {

		primesbitlen[uint64(qi)]++

		if uint64(params.Modulichain[i]) > 60 {
			panic("provided moduli must be smaller than 61")
		}
	}

	for _, pj := range params.P {
		primesbitlen[uint64(pj)]++

		if uint64(pj) > 60 {
			panic("provided P must be smaller than 61")
		}
	}

	// For each bitsize, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = GenerateCKKSPrimes(key, uint64(params.LogN), value)
	}

	// Assigns the primes to the ckks moduli chain
	moduli := make([]uint64, len(params.Modulichain))
	for i, qi := range params.Modulichain {
		moduli[i] = primes[uint64(params.Modulichain[i])][0]
		primes[uint64(qi)] = primes[uint64(qi)][1:]

		scalechain[i] = float64(moduli[i])
	}

	// Assigns the primes to the special primes list for the the keyscontext
	specialPrimes := make([]uint64, len(params.P))
	for i, pj := range params.P {
		specialPrimes[i] = primes[uint64(pj)][0]
		primes[uint64(pj)] = primes[uint64(pj)][1:]
	}

	contextKeys := ring.NewContext()
	contextKeys.SetParameters(1<<params.LogN, append(moduli, specialPrimes...))

	err := contextKeys.GenNTTParams()
	if err != nil {
		panic(err)
	}

	gaussianSampler := contextKeys.NewKYSampler(params.Sigma, int(6*params.Sigma))

	var m, mask uint64

	m = n << 1

	mask = m - 1

	gen := uint64(5) // Any integer equal to 1 mod 4 and comprime to 2N will do fine
	genInv := ring.ModExp(gen, mask, m)

	galElRotColLeft := make([]uint64, maxSlots)
	galElRotColRight := make([]uint64, maxSlots)

	galElRotColRight[0] = 1
	galElRotColLeft[0] = 1

	for i := uint64(1); i < maxSlots; i++ {
		galElRotColLeft[i] = (galElRotColLeft[i-1] * gen) & mask
		galElRotColRight[i] = (galElRotColRight[i-1] * genInv) & mask
	}

	galElRotRow := mask

	return &keyGeneratorContext{
		logN:             logN,
		n:                n,
		levels:           levels,
		specialPrimes:    specialPrimes,
		alpha:            alpha,
		beta:             beta,
		contextKeys:      contextKeys,
		gaussianSampler:  gaussianSampler,
		galElRotRow:      galElRotRow,
		galElRotColLeft:  galElRotColLeft,
		galElRotColRight: galElRotColRight,
	}
}

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	context     *keyGeneratorContext
	ringContext *ring.Context
	polypool    *ring.Poly
}

// SecretKey is a structure that stores the SecretKey
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the PublicKey
type PublicKey struct {
	pk [2]*ring.Poly
}

// Rotation is a type used to represent the rotations types.
type Rotation int

const (
	RotationRight = iota + 1
	RotationLeft
	Conjugate
)

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	permuteNTTRightIndex     map[uint64][]uint64
	permuteNTTLeftIndex      map[uint64][]uint64
	permuteNTTConjugateIndex []uint64

	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
	evakeyConjugate   *SwitchingKey
}

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey *SwitchingKey
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	evakey [][2]*ring.Poly
}

// Get returns the switching key backing slice
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.evakey
}

// NewKeyGenerator creates a new keygenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params *Parameters) (keygen *KeyGenerator) {
	context := newKeyGeneratorContext(params)

	keygen = new(KeyGenerator)
	keygen.context = context
	keygen.ringContext = context.contextKeys
	keygen.polypool = keygen.ringContext.NewPoly()
	return
}

// NewSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) NewSecretKey() (sk *SecretKey) {
	return keygen.NewSecretKeyWithDistrib(1.0 / 3)
}

// NewSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *KeyGenerator) NewSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.context.contextKeys.SampleTernaryMontgomeryNTTNew(p)
	return sk
}

// NewSecretKeySparse generates a new SecretKey with exactly hw non zero coefficients.
func (keygen *KeyGenerator) NewSecretKeySparse(hw uint64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.context.contextKeys.SampleTernarySparseMontgomeryNTTNew(hw)
	return sk
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params *Parameters) *SecretKey {
	context := newKeyGeneratorContext(params)

	sk := new(SecretKey)
	sk.sk = context.contextKeys.NewPoly()
	return sk
}

// Get returns the SecretKey value of the SecretKey.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the value of the SecretKey to the provided value.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// NewPublicKey generates a new public key from the provided SecretKey.
func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.context.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.ringContext.NewUniformPoly()

	keygen.ringContext.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.ringContext.Neg(pk.pk[0], pk.pk[0])

	return pk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {
	context := newKeyGeneratorContext(params)

	pk = new(PublicKey)

	pk.pk[0] = context.contextKeys.NewPoly()
	pk.pk[1] = context.contextKeys.NewPoly()

	return
}

// Get returns the value of the the public key.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the value of the public key to the provided value.
func (pk *PublicKey) Set(poly [2]*ring.Poly) {
	pk.pk[0] = poly[0].CopyNew()
	pk.pk[1] = poly[1].CopyNew()
}

// NewKeyPair generates a new secretkey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *KeyGenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.NewSecretKey()
	return sk, keygen.NewPublicKey(sk)
}

// NewKeyPairSparse generates a new secretkey with exactly hw non zero coefficients [1/2, 0, 1/2].
func (keygen *KeyGenerator) NewKeyPairSparse(hw uint64) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.NewSecretKeySparse(hw)
	return sk, keygen.NewPublicKey(sk)
}

// NewRelinKey generates a new evaluation key that will be used to relinearize the ciphertexts during multiplication.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)
	keygen.polypool.Copy(sk.Get())
	keygen.ringContext.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
	evakey.evakey = keygen.newSwitchingKey(keygen.polypool, sk.Get())
	keygen.polypool.Zero()

	return
}

// NewRelinKey returns  new EvaluationKey with zero values.
func NewRelinKey(params *Parameters) (evakey *EvaluationKey) {
	context := newKeyGeneratorContext(params)

	evakey = new(EvaluationKey)
	evakey.evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey.evakey = make([][2]*ring.Poly, context.beta)
	for i := uint64(0); i < context.beta; i++ {

		evakey.evakey.evakey[i][0] = context.contextKeys.NewPoly()
		evakey.evakey.evakey[i][1] = context.contextKeys.NewPoly()
	}

	return
}

// Get returns the slice of switchintkeys of the evaluation-key.
func (evk *EvaluationKey) Get() *SwitchingKey {
	return evk.evakey
}

// Set sets the target Evaluation key with the input polynomials.
func (evk *EvaluationKey) Set(rlk [][2]*ring.Poly) {

	evk.evakey = new(SwitchingKey)
	evk.evakey.evakey = make([][2]*ring.Poly, len(rlk))
	for j := range rlk {
		evk.evakey.evakey[j][0] = rlk[j][0].CopyNew()
		evk.evakey.evakey[j][1] = rlk[j][1].CopyNew()
	}
}

// NewSwitchingKey generated a new keyswitching key, that will re-encrypt a ciphertext encrypted under the input key to the output key.
func (keygen *KeyGenerator) NewSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey) {
	keygen.ringContext.Sub(skInput.Get(), skOutput.Get(), keygen.polypool)
	newevakey = keygen.newSwitchingKey(keygen.polypool, skOutput.Get())
	keygen.polypool.Zero()
	return
}

// NewSwitchingKey returns a new SwitchingKey with zero values.
func NewSwitchingKey(params *Parameters) (evakey *SwitchingKey) {
	context := newKeyGeneratorContext(params)

	evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey = make([][2]*ring.Poly, context.beta)

	for i := uint64(0); i < context.beta; i++ {
		evakey.evakey[i][0] = context.contextKeys.NewPoly()
		evakey.evakey[i][1] = context.contextKeys.NewPoly()
	}

	return
}

func (keygen *KeyGenerator) newSwitchingKey(skIn, skOut *ring.Poly) (switchingkey *SwitchingKey) {

	switchingkey = new(SwitchingKey)

	context := keygen.context.contextKeys

	// Computes P * skIn
	for _, pj := range keygen.context.specialPrimes {
		context.MulScalar(skIn, pj, skIn)
	}

	alpha := keygen.context.alpha
	beta := keygen.context.beta

	var index uint64

	switchingkey.evakey = make([][2]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {

		// e
		switchingkey.evakey[i][0] = keygen.context.gaussianSampler.SampleNTTNew()
		context.MForm(switchingkey.evakey[i][0], switchingkey.evakey[i][0])

		// a (since a is uniform, we consider we already sample it in the NTT and montgomery domain)
		switchingkey.evakey[i][1] = keygen.ringContext.NewUniformPoly()

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := context.Modulus[index]
			p0tmp := skIn.Coeffs[index]
			p1tmp := switchingkey.evakey[i][0].Coeffs[index]

			for w := uint64(0); w < context.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= keygen.context.levels-1 {
				break
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		context.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][1], skOut, switchingkey.evakey[i][0])
	}

	return
}

// NewRotationKeys generates a new instance of rotationkeys, with the provided rotation to the left, right and conjugation if asked.
func NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	return
}

// GenRot populates input RotationKeys with a SwitchingKey for the given rotation type and amount.
func (keygen *KeyGenerator) GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {
	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(keygen.context.galElRotColLeft[k], 1<<keygen.context.logN)
			rotKey.evakeyRotColLeft[k] = keygen.genrotKey(sk.Get(), keygen.context.galElRotColLeft[k])
		}

	case RotationRight:

		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTRightIndex == nil {
			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(keygen.context.galElRotColRight[k], 1<<keygen.context.logN)
			rotKey.evakeyRotColRight[k] = keygen.genrotKey(sk.Get(), keygen.context.galElRotColRight[k])
		}

	case Conjugate:
		rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(keygen.context.galElRotRow, 1<<keygen.context.logN)
		rotKey.evakeyConjugate = keygen.genrotKey(sk.Get(), keygen.context.galElRotRow)
	}
}

// NewRotationKeysPow2 generates a new rotation key with all the power of two rotation to the left and right, as well as the conjugation.
func (keygen *KeyGenerator) NewRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)

	rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
	rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)

	rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
	rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)

	for n := uint64(1); n < keygen.context.n>>1; n <<= 1 {

		rotKey.permuteNTTLeftIndex[n] = ring.PermuteNTTIndex(keygen.context.galElRotColLeft[n], 1<<keygen.context.logN)
		rotKey.permuteNTTRightIndex[n] = ring.PermuteNTTIndex(keygen.context.galElRotColRight[n], 1<<keygen.context.logN)

		rotKey.evakeyRotColLeft[n] = keygen.genrotKey(skOutput.Get(), keygen.context.galElRotColLeft[n])
		rotKey.evakeyRotColRight[n] = keygen.genrotKey(skOutput.Get(), keygen.context.galElRotColRight[n])
	}

	rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(keygen.context.galElRotRow, 1<<keygen.context.logN)
	rotKey.evakeyConjugate = keygen.genrotKey(skOutput.Get(), keygen.context.galElRotRow)
	return
}

// SetRotKey sets the target RotationKeys' SwitchingKey for the specified rotation type and amount with the input polynomials.
func (ckkscontext *Context) SetRotKey(evakey [][2]*ring.Poly, rotType Rotation, k uint64, rotKey *RotationKeys) {
	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {

			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(ckkscontext.galElRotColLeft[k], 1<<ckkscontext.logN)

			rotKey.evakeyRotColLeft[k] = new(SwitchingKey)
			rotKey.evakeyRotColLeft[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColLeft[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColLeft[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}

	case RotationRight:

		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {

			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(ckkscontext.galElRotColRight[k], 1<<ckkscontext.logN)

			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}

	case Conjugate:

		if rotKey.evakeyConjugate == nil {

			rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(ckkscontext.galElRotRow, 1<<ckkscontext.logN)

			rotKey.evakeyConjugate = new(SwitchingKey)
			rotKey.evakeyConjugate.evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyConjugate.evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyConjugate.evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	}
}

func (keygen *KeyGenerator) genrotKey(skOutput *ring.Poly, gen uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(skOutput, gen, keygen.polypool)
	keygen.ringContext.Sub(keygen.polypool, skOutput, keygen.polypool)
	switchingkey = keygen.newSwitchingKey(keygen.polypool, skOutput)
	keygen.polypool.Zero()

	return
}
