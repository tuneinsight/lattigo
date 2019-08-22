package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	"math"
	"math/bits"
)

type SecretKey struct {
	sk *ring.Poly
}

type PublicKey struct {
	pk [2]*ring.Poly
}

type KeyGenerator struct {
	bfvcontext *BfvContext
	context    *ring.Context
	polypool   *ring.Poly
}

type RotationKeys struct {
	bfvcontext       *BfvContext
	bitDecomp        uint64
	evakey_rot_col_L map[uint64]*SwitchingKey
	evakey_rot_col_R map[uint64]*SwitchingKey
	evakey_rot_row   *SwitchingKey
}

type EvaluationKey struct {
	evakey []*SwitchingKey
}

type SwitchingKey struct {
	bitDecomp uint64
	evakey    [][][2]*ring.Poly
}

func (bfvcontext *BfvContext) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.bfvcontext = bfvcontext
	keygen.context = bfvcontext.contextQ
	keygen.polypool = keygen.context.NewPoly()
	return
}

func (keygen *KeyGenerator) NewSecretKey() *SecretKey {

	sk := new(SecretKey)
	sk.sk = keygen.bfvcontext.ternarySampler.SampleMontgomeryNTTNew()

	return sk
}

func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

func (keygen *KeyGenerator) check_sk(sk_output *SecretKey) error {

	if sk_output.Get().GetDegree() != int(keygen.context.N) {
		return errors.New("error : pol degree sk != bfvcontext.n")
	}

	if len(sk_output.Get().Coeffs) != len(keygen.context.Modulus) {
		return errors.New("error : nb modulus sk != nb modulus bfvcontext")
	}

	return nil
}

func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (*PublicKey, error) {

	if err := keygen.check_sk(sk); err != nil {
		return nil, err
	}

	pk := new(PublicKey)

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.bfvcontext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.context.NewUniformPoly()

	keygen.context.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.context.Neg(pk.pk[0], pk.pk[0])

	return pk, nil
}

func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

func (pk *PublicKey) Set(p [2]*ring.Poly) {
	pk.pk[0] = p[0].CopyNew()
	pk.pk[1] = p[1].CopyNew()
}

func (keygen *KeyGenerator) SetPublicKey(p [2]*ring.Poly) (*PublicKey, error) {

	pk := new(PublicKey)

	pk.pk[0] = p[0].CopyNew()
	pk.pk[1] = p[1].CopyNew()

	return pk, nil
}

func (keygen *KeyGenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey, err error) {
	sk = keygen.NewSecretKey()
	pk, err = keygen.NewPublicKey(sk)
	return
}

func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey, maxDegree, bitDecomp uint64) (newEvakey *EvaluationKey, err error) {

	newEvakey = new(EvaluationKey)
	newEvakey.evakey = make([]*SwitchingKey, maxDegree)
	sk.Get().Copy(keygen.polypool)
	for i := uint64(0); i < maxDegree; i++ {
		keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
		newEvakey.evakey[i] = newswitchingkey(keygen.bfvcontext, keygen.polypool, sk.Get(), bitDecomp)
	}
	keygen.polypool.Zero()

	return newEvakey, nil
}

func (evk *EvaluationKey) Get() []*SwitchingKey {
	return evk.evakey
}

func (keygen *KeyGenerator) SetRelinKeys(rlk [][][][2]*ring.Poly, bitDecomp uint64) (*EvaluationKey, error) {

	newevakey := new(EvaluationKey)

	newevakey.evakey = make([]*SwitchingKey, len(rlk))
	for i := range rlk {
		newevakey.evakey[i] = new(SwitchingKey)
		newevakey.evakey[i].bitDecomp = bitDecomp
		newevakey.evakey[i].evakey = make([][][2]*ring.Poly, len(rlk[i]))
		for j := range rlk[i] {
			newevakey.evakey[i].evakey[j] = make([][2]*ring.Poly, len(rlk[i][j]))
			for u := range rlk[i][j] {
				newevakey.evakey[i].evakey[j][u][0] = rlk[i][j][u][0].CopyNew()
				newevakey.evakey[i].evakey[j][u][1] = rlk[i][j][u][1].CopyNew()
			}
		}
	}

	return newevakey, nil
}

func (keygen *KeyGenerator) NewSwitchingKey(sk_input, sk_output *SecretKey, bitDecomp uint64) (newevakey *SwitchingKey, err error) {

	if err = keygen.check_sk(sk_input); err != nil {
		return nil, err
	}

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	keygen.context.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)
	newevakey = newswitchingkey(keygen.bfvcontext, keygen.polypool, sk_output.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) NewRotationKeys(sk_output *SecretKey, bitDecomp uint64, rotLeft []uint64, rotRight []uint64, conjugate bool) (rotKey *RotationKeys, err error) {

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	rotKey = new(RotationKeys)
	rotKey.bfvcontext = keygen.bfvcontext
	rotKey.bitDecomp = bitDecomp

	if rotLeft != nil {
		rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
		for _, n := range rotLeft {
			if rotKey.evakey_rot_col_L[n] == nil && n != 0 {
				rotKey.evakey_rot_col_L[n] = genrotkey(keygen, sk_output.Get(), keygen.bfvcontext.galElRotColLeft[n], bitDecomp)
			}
		}
	}

	if rotRight != nil {
		rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)
		for _, n := range rotRight {
			if rotKey.evakey_rot_col_R[n] == nil && n != 0 {
				rotKey.evakey_rot_col_R[n] = genrotkey(keygen, sk_output.Get(), keygen.bfvcontext.galElRotColRight[n], bitDecomp)
			}
		}
	}

	if conjugate {
		rotKey.evakey_rot_row = genrotkey(keygen, sk_output.Get(), keygen.bfvcontext.galElRotRow, bitDecomp)
	}

	return rotKey, nil

}

func (keygen *KeyGenerator) NewRotationKeysPow2(sk_output *SecretKey, bitDecomp uint64, conjugate bool) (rotKey *RotationKeys, err error) {

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	rotKey = new(RotationKeys)
	rotKey.bfvcontext = keygen.bfvcontext
	rotKey.bitDecomp = bitDecomp

	rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
	rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)

	for n := uint64(1); n < rotKey.bfvcontext.n>>1; n <<= 1 {

		rotKey.evakey_rot_col_L[n] = genrotkey(keygen, sk_output.Get(), keygen.bfvcontext.galElRotColLeft[n], bitDecomp)
		rotKey.evakey_rot_col_R[n] = genrotkey(keygen, sk_output.Get(), keygen.bfvcontext.galElRotColRight[n], bitDecomp)
	}

	if conjugate {
		rotKey.evakey_rot_row = genrotkey(keygen, sk_output.Get(), keygen.bfvcontext.galElRotRow, bitDecomp)
	}

	return
}

func genrotkey(keygen *KeyGenerator, sk_output *ring.Poly, gen, bitDecomp uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(sk_output, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, sk_output, keygen.polypool)
	switchingkey = newswitchingkey(keygen.bfvcontext, keygen.polypool, sk_output, bitDecomp)
	keygen.polypool.Zero()

	return
}

func newswitchingkey(bfvcontext *BfvContext, sk_in, sk_out *ring.Poly, bitDecomp uint64) (switchingkey *SwitchingKey) {

	if bitDecomp > bfvcontext.maxBit || bitDecomp == 0 {
		bitDecomp = bfvcontext.maxBit
	}

	switchingkey = new(SwitchingKey)

	context := bfvcontext.contextQ

	switchingkey.bitDecomp = uint64(bitDecomp)

	mredParams := context.GetMredParams()

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output

	var bitLog uint64

	switchingkey.evakey = make([][][2]*ring.Poly, len(context.Modulus))

	for i, qi := range context.Modulus {

		bitLog = uint64(math.Ceil(float64(bits.Len64(qi)) / float64(bitDecomp)))

		switchingkey.evakey[i] = make([][2]*ring.Poly, bitLog)

		for j := uint64(0); j < bitLog; j++ {

			// e
			switchingkey.evakey[i][j][0] = bfvcontext.gaussianSampler.SampleNTTNew()
			// a
			switchingkey.evakey[i][j][1] = context.NewUniformPoly()

			// e + sk_in * (qiBarre*qiStar) * 2^w
			// (qiBarre*qiStar)%qi = 1, else 0
			for w := uint64(0); w < context.N; w++ {
				switchingkey.evakey[i][j][0].Coeffs[i][w] += PowerOf2(sk_in.Coeffs[i][w], bitDecomp*j, qi, mredParams[i])
			}

			// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
			context.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][j][1], sk_out, switchingkey.evakey[i][j][0])

			context.MForm(switchingkey.evakey[i][j][0], switchingkey.evakey[i][j][0])
			context.MForm(switchingkey.evakey[i][j][1], switchingkey.evakey[i][j][1])
		}
	}

	return
}
