package bfv

import (
	"github.com/lca1/lattigo-private/ring"
	"math/bits"
)

type BfvContext struct {

	// Polynomial degree
	n uint64

	// Plaintext Modulus
	t uint64

	// floor(Q/T) mod each Qi in montgomeryform
	DeltaMont []uint64
	Delta     []uint64

	// Ternary and Gaussian samplers
	sigma           float64
	gaussianSampler *ring.KYSampler
	ternarySampler  *ring.TernarySampler

	// Maximum bitLength of the modulies
	maxBit uint64

	// Polynomial contexts
	contextT  *ring.Context
	contextQ  *ring.Context
	contextP  *ring.Context
	contextQP *ring.Context

	// Galois elements used to permute the batched plaintext in the encrypted domain
	gen    uint64
	genInv uint64

	galElRotRow      uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64

	// Checksum of [N, [modulies]]
	checksum []byte
}

func NewBfvContext() *BfvContext {
	return new(BfvContext)
}

func NewBfvContextWithParam(N, t uint64, ModulieQ, ModulieP []uint64, sigma float64) (newbfvcontext *BfvContext, err error) {
	newbfvcontext = new(BfvContext)
	if err := newbfvcontext.SetParameters(N, t, ModulieQ, ModulieP, sigma); err != nil {
		return nil, err
	}
	return
}

func (bfvContext *BfvContext) SetParameters(N, t uint64, ModulieQ, ModulieP []uint64, sigma float64) (err error) {

	contextT := ring.NewContext()
	contextQ := ring.NewContext()
	contextP := ring.NewContext()
	contextQP := ring.NewContext()

	// Plaintext NTT Parameters
	// We do not check for an error since the plaintext NTT is optional
	// it will still compute the other relevant parameters
	contextT.SetParameters(N, []uint64{t})
	contextT.ValidateParameters()
	// ========================

	if err := contextQ.SetParameters(N, ModulieQ); err != nil {
		return err
	}

	if err := contextQ.ValidateParameters(); err != nil {
		return err
	}

	if err := contextP.SetParameters(N, ModulieP); err != nil {
		return err
	}

	if err := contextP.ValidateParameters(); err != nil {
		return err
	}

	if err := contextQP.Merge(contextQ, contextP); err != nil {
		return err
	}

	bfvContext.n = N

	bfvContext.t = t

	delta := ring.NewUint(1).Div(contextQ.ModulusBigint, ring.NewUint(t))
	tmpBig := ring.NewUint(1)
	bfvContext.DeltaMont = make([]uint64, len(ModulieQ))
	bfvContext.Delta = make([]uint64, len(ModulieQ))
	for i, Qi := range ModulieQ {
		bfvContext.Delta[i] = tmpBig.Mod(delta, ring.NewUint(Qi)).Uint64()
		bfvContext.DeltaMont[i] = ring.MForm(bfvContext.Delta[i], Qi, contextQ.GetBredParams()[i])
	}

	bfvContext.sigma = sigma

	bfvContext.gaussianSampler = contextQ.NewKYSampler(sigma, int(6*sigma))
	bfvContext.ternarySampler = contextQ.NewTernarySampler()

	bfvContext.maxBit = 0

	for i := range ModulieQ {
		tmp := uint64(bits.Len64(ModulieQ[i]))
		if tmp > bfvContext.maxBit {
			bfvContext.maxBit = tmp
		}
	}

	bfvContext.contextT = contextT
	bfvContext.contextQ = contextQ
	bfvContext.contextP = contextP
	bfvContext.contextQP = contextQP

	bfvContext.gen = 5
	bfvContext.genInv = modexp(bfvContext.gen, (N<<1)-1, N<<1)

	mask := (N << 1) - 1

	bfvContext.galElRotColLeft = make([]uint64, N>>1)
	bfvContext.galElRotColRight = make([]uint64, N>>1)

	bfvContext.galElRotColRight[0] = 1
	bfvContext.galElRotColLeft[0] = 1

	for i := uint64(1); i < N>>1; i++ {
		bfvContext.galElRotColLeft[i] = (bfvContext.galElRotColLeft[i-1] * bfvContext.gen) & mask
		bfvContext.galElRotColRight[i] = (bfvContext.galElRotColRight[i-1] * bfvContext.genInv) & mask

	}

	bfvContext.galElRotRow = (N << 1) - 1

	toHash := make([]uint64, len(ModulieQ)+1)
	toHash[0] = N
	for i := 1; i < len(toHash); i++ {
		toHash[i] = ModulieQ[i-1]
	}

	if bfvContext.checksum, err = Hash(toHash); err != nil {
		return err
	}

	return nil
}

func (bfvContext *BfvContext) GetN() uint64 {
	return bfvContext.n
}

func (bfvContext *BfvContext) GetPlaintextModulus() uint64 {
	return bfvContext.t
}

func (bfvContext *BfvContext) GetDelta() []uint64 {
	return bfvContext.DeltaMont
}

func (bfvContext *BfvContext) GetSigma() float64 {
	return bfvContext.sigma
}

func (bfvContext *BfvContext) GetContextT() *ring.Context {
	return bfvContext.contextT
}

func (bfvContext *BfvContext) GetContextQ() *ring.Context {
	return bfvContext.contextQ
}

func (bfvContext *BfvContext) GetContextP() *ring.Context {
	return bfvContext.contextP
}

func (bfvContext *BfvContext) GetContextQP() *ring.Context {
	return bfvContext.contextQP
}
