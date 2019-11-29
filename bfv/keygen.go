package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	params     *Parameters
	bfvContext *Context
	polypool   *ring.Poly
}

// SecretKey is a structure that stores the SecretKey
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the PublicKey
type PublicKey struct {
	pk [2]*ring.Poly
}

// Rotation is a constant used to represent the rotations types.
type Rotation int

const (
	RotationRight = iota + 1
	RotationLeft
	RotationRow
)

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
	evakeyRotRow      *SwitchingKey
}

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey []*SwitchingKey
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	evakey [][2]*ring.Poly
}

// Get returns the switching key backing slice
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.evakey
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params *Parameters) (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	params = params.Copy()
	keygen.bfvContext = NewContext(params)
	keygen.polypool = keygen.bfvContext.contextQ1P.NewPoly()
	return
}

// NewSecretKey creates a new SecretKey with the distribution [1/3, 1/3, 1/3]
func (keygen *KeyGenerator) NewSecretKey() (sk *SecretKey) {
	return keygen.NewSecretkeyWithDistrib(1.0 / 3)
}

// NewSecretkeyWithDistrib creates a new SecretKey with the distribution [(1-p)/2, p, (1-p)/2]
func (keygen *KeyGenerator) NewSecretkeyWithDistrib(p float64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.bfvContext.contextQ1P.SampleTernaryMontgomeryNTTNew(p)
	return sk
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params *Parameters) *SecretKey {
	sk := new(SecretKey)
	sk.sk = ring.NewPoly(uint64(1<<params.LogN), uint64(len(params.Q1)+len(params.P)))
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

// NewPublicKey generates a new publickkey from the provided SecretKey
func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringContext := keygen.bfvContext.contextQ1P

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.bfvContext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = ringContext.NewUniformPoly()

	ringContext.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	ringContext.Neg(pk.pk[0], pk.pk[0])

	return pk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {
	pk = new(PublicKey)

	pk.pk[0] = ring.NewPoly(uint64(1<<params.LogN), uint64(len(params.Q1)+len(params.P)))
	pk.pk[1] = ring.NewPoly(uint64(1<<params.LogN), uint64(len(params.Q1)+len(params.P)))

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
func (keygen *KeyGenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.NewSecretKey()
	return sk, keygen.NewPublicKey(sk)
}

// NewRelinKey generates a new evaluation key from the provided SecretKey. It will be used to relinearize a ciphertext (encrypted under a PublicKey generated from the provided SecretKey)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey, maxDegree uint64) (evk *EvaluationKey) {

	evk = new(EvaluationKey)

	evk.evakey = make([]*SwitchingKey, maxDegree)

	keygen.polypool.Copy(sk.Get())

	ringContext := keygen.bfvContext.contextQ1P

	ringContext.MulScalarBigint(keygen.polypool, keygen.bfvContext.contextP.ModulusBigint, keygen.polypool)

	for i := uint64(0); i < maxDegree; i++ {
		ringContext.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
		evk.evakey[i] = newswitchintkey(keygen.bfvContext, keygen.polypool, sk.Get())
	}

	keygen.polypool.Zero()

	return
}

// NewRelinKey creates a new EvaluationKey with zero values.
func NewRelinKey(params *Parameters, maxDegree uint64) (evakey *EvaluationKey) {
	evakey = new(EvaluationKey)

	beta := uint64(math.Ceil(float64(len(params.Q1)) / float64(len(params.P))))

	evakey.evakey = make([]*SwitchingKey, maxDegree)

	for w := uint64(0); w < maxDegree; w++ {

		evakey.evakey[w] = new(SwitchingKey)

		evakey.evakey[w].evakey = make([][2]*ring.Poly, beta)

		for i := uint64(0); i < beta; i++ {

			evakey.evakey[w].evakey[i][0] = ring.NewPoly(uint64(1<<params.LogN), uint64(len(params.Q1)+len(params.P)))
			evakey.evakey[w].evakey[i][1] = ring.NewPoly(uint64(1<<params.LogN), uint64(len(params.Q1)+len(params.P)))
		}
	}

	return
}

// Get returns the slice of SwitchingKeys of the target EvaluationKey.
func (evk *EvaluationKey) Get() []*SwitchingKey {
	return evk.evakey
}

// SetRelinKeys sets the polynomial of the target EvaluationKey as the input polynomials.
func (evk *EvaluationKey) SetRelinKeys(rlk [][][2]*ring.Poly) {

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

// NewSwitchingKey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key.
func (keygen *KeyGenerator) NewSwitchingKey(skIn, skOut *SecretKey) (evk *SwitchingKey) {

	ringContext := keygen.bfvContext.contextQ1P

	ringContext.Sub(skIn.Get(), skOut.Get(), keygen.polypool)

	ringContext.MulScalarBigint(keygen.polypool, keygen.bfvContext.contextP.ModulusBigint, keygen.polypool)

	evk = newswitchintkey(keygen.bfvContext, keygen.polypool, skOut.Get())
	keygen.polypool.Zero()

	return
}

// NewSwitchingKey returns a new SwitchingKey with zero values.
func NewSwitchingKey(params *Parameters) (evakey *SwitchingKey) {
	evakey = new(SwitchingKey)

	beta := uint64(math.Ceil(float64(len(params.Q1)) / float64(len(params.P))))

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey = make([][2]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {
		evakey.evakey[i][0] = ring.NewPoly(uint64(1<<params.LogN), uint64(len(params.Q1)+len(params.P)))
		evakey.evakey[i][1] = ring.NewPoly(uint64(1<<params.LogN), uint64(len(params.Q1)+len(params.P)))
	}

	return
}

func newswitchintkey(bfvContext *Context, skIn, skOut *ring.Poly) (switchkey *SwitchingKey) {

	switchkey = new(SwitchingKey)

	ringContext := bfvContext.contextQ1P

	var index uint64

	// delta_sk = skIn - skOut = GaloisEnd(skOut, rotation) - skOut

	switchkey.evakey = make([][2]*ring.Poly, bfvContext.beta)

	for i := uint64(0); i < bfvContext.beta; i++ {

		// e
		switchkey.evakey[i][0] = bfvContext.gaussianSampler.SampleNTTNew()
		ringContext.MForm(switchkey.evakey[i][0], switchkey.evakey[i][0])
		// a
		switchkey.evakey[i][1] = ringContext.NewUniformPoly()

		// e + skIn * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0

		for j := uint64(0); j < bfvContext.alpha; j++ {

			index = i*bfvContext.alpha + j

			qi := ringContext.Modulus[index]
			p0tmp := skIn.Coeffs[index]
			p1tmp := switchkey.evakey[i][0].Coeffs[index]

			for w := uint64(0); w < ringContext.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= uint64(len(ringContext.Modulus)-1) {
				break
			}

		}

		// skIn * (qiBarre*qiStar) * 2^w - a*sk + e
		ringContext.MulCoeffsMontgomeryAndSub(switchkey.evakey[i][1], skOut, switchkey.evakey[i][0])
	}

	return
}

// NewRotationKeys returns a new empty RotationKeys struct.
func NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	return
}

// GenRot populates the target RotationKeys with a SwitchingKey for the desired rotation type and amount.
func (keygen *KeyGenerator) GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {

	k &= ((keygen.bfvContext.n >> 1) - 1)

	switch rotType {
	case RotationLeft:
		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
			rotKey.evakeyRotColLeft[k] = genrotkey(keygen, sk.Get(), keygen.bfvContext.galElRotColLeft[k])
		}
	case RotationRight:
		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.evakeyRotColRight[k] = genrotkey(keygen, sk.Get(), keygen.bfvContext.galElRotColRight[k])
		}
	case RotationRow:
		rotKey.evakeyRotRow = genrotkey(keygen, sk.Get(), keygen.bfvContext.galElRotRow)
	}
}

// NewRotationKeysPow2 generates a new struct of rotationkeys storing the keys of all the left and right powers of two rotations. The provided SecretKey must be the SecretKey used to generate the PublicKey under
// which the ciphertexts to rotate are encrypted under. rows is a boolean value indicatig if the keys for the row rotation have to be generated.
func (keygen *KeyGenerator) NewRotationKeysPow2(sk *SecretKey) (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)

	rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
	rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)

	for n := uint64(1); n < keygen.bfvContext.n>>1; n <<= 1 {

		rotKey.evakeyRotColLeft[n] = genrotkey(keygen, sk.Get(), keygen.bfvContext.galElRotColLeft[n])
		rotKey.evakeyRotColRight[n] = genrotkey(keygen, sk.Get(), keygen.bfvContext.galElRotColRight[n])
	}

	rotKey.evakeyRotRow = genrotkey(keygen, sk.Get(), keygen.bfvContext.galElRotRow)

	return
}

// SetRotKey populates the target RotationKeys with a new SwitchingKey using the input polynomials.
func (rotKey *RotationKeys) SetRotKey(rotType Rotation, k uint64, evakey [][2]*ring.Poly) {
	switch rotType {
	case RotationLeft:
		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
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
		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	case RotationRow:
		if rotKey.evakeyRotRow == nil {
			rotKey.evakeyRotRow = new(SwitchingKey)
			rotKey.evakeyRotRow.evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotRow.evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotRow.evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	}
}

func genrotkey(keygen *KeyGenerator, sk *ring.Poly, gen uint64) (switchkey *SwitchingKey) {

	ringContext := keygen.bfvContext.contextQ1P

	ring.PermuteNTT(sk, gen, keygen.polypool)
	ringContext.Sub(keygen.polypool, sk, keygen.polypool)

	ringContext.MulScalarBigint(keygen.polypool, keygen.bfvContext.contextP.ModulusBigint, keygen.polypool)

	switchkey = newswitchintkey(keygen.bfvContext, keygen.polypool, sk)
	keygen.polypool.Zero()

	return
}
