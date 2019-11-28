package dbfv

import (
	"math/big"

	"encoding/binary"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	//"fmt"
)

type refreshProtocolContext struct {
	// Polynomial degree
	n uint64

	// Plaintext Modulus
	t uint64

	// floor(Q/T) mod each Qi in Montgomery form
	deltaMont []uint64

	// Polynomial contexts
	contextT *ring.Context
	contextQ *ring.Context

	contextKeys       *ring.Context
	contextPKeys      *ring.Context
	specialPrimes     []uint64
	rescaleParamsKeys []uint64 // (P^-1) mod each qi
}

func newRefreshProtocolContext(params *bfv.Parameters) *refreshProtocolContext {
	n := params.N
	t := params.T

	contextT := ring.NewContext()
	contextT.SetParameters(n, []uint64{t})
	err := contextT.GenNTTParams()
	if err != nil {
		panic(err)
	}

	contextQ := ring.NewContext()
	contextQ.SetParameters(n, params.Qi)
	err = contextQ.GenNTTParams()
	if err != nil {
		panic(err)
	}

	contextKeys := ring.NewContext()
	contextKeys.SetParameters(n, append(params.Qi, params.KeySwitchPrimes...))
	err = contextKeys.GenNTTParams()
	if err != nil {
		panic(err)
	}

	contextPKeys := ring.NewContext()
	contextPKeys.SetParameters(n, params.KeySwitchPrimes)
	err = contextPKeys.GenNTTParams()
	if err != nil {
		panic(err)
	}

	specialPrimes := make([]uint64, len(params.KeySwitchPrimes))
	for i := range params.KeySwitchPrimes {
		specialPrimes[i] = params.KeySwitchPrimes[i]
	}

	rescaleParamsKeys := make([]uint64, len(params.Qi))

	PBig := ring.NewUint(1)
	for _, pj := range specialPrimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := new(big.Int)
	bredParams := contextQ.GetBredParams()
	for i, Qi := range params.Qi {
		tmp.Mod(PBig, ring.NewUint(Qi))
		rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	delta0 := new(big.Int).Quo(contextQ.ModulusBigint, ring.NewUint(t))
	tmpBig := new(big.Int)
	deltaMont := make([]uint64, len(params.Qi))
	delta := make([]uint64, len(params.Qi))
	for i, Qi := range params.Qi {
		delta[i] = tmpBig.Mod(delta0, ring.NewUint(Qi)).Uint64()
		deltaMont[i] = ring.MForm(delta[i], Qi, contextQ.GetBredParams()[i])
	}

	return &refreshProtocolContext{
		n:                 n,
		t:                 t,
		deltaMont:         deltaMont,
		contextT:          contextT,
		contextQ:          contextQ,
		contextKeys:       contextKeys,
		contextPKeys:      contextPKeys,
		specialPrimes:     specialPrimes,
		rescaleParamsKeys: rescaleParamsKeys,
	}
}

// RefreshProtocol is a struct storing the relevant parameters for the Refresh protocol.
type RefreshProtocol struct {
	context       *refreshProtocolContext
	tmp1          *ring.Poly
	tmp2          *ring.Poly
	hP            *ring.Poly
	baseconverter *ring.FastBasisExtender
}

// RefreshShareDecrypt is a struct storing the decrpytion share.
type RefreshShareDecrypt *ring.Poly

// RefreshShareRecrypt is a struct storing the recrpytion share.
type RefreshShareRecrypt *ring.Poly

// RefreshShare is a struct storing the decryption and recryption shares.
type RefreshShare struct {
	RefreshShareDecrypt RefreshShareDecrypt
	RefreshShareRecrypt RefreshShareRecrypt
}

// MarshalBinary encodes a RefreshShare on a slice of bytes.
func (share *RefreshShare) MarshalBinary() ([]byte, error) {
	lenDecrypt := (*share.RefreshShareDecrypt).GetDataLen(true)
	lenRecrypt := (*share.RefreshShareRecrypt).GetDataLen(true)

	data := make([]byte, lenDecrypt+lenRecrypt+2*8) // 2 * 3 to write the len of lenDecrypt and lenRecrypt.
	binary.BigEndian.PutUint64(data[0:8], lenDecrypt)
	binary.BigEndian.PutUint64(data[8:16], lenRecrypt)

	ptr := uint64(16)
	tmp, err := (*share.RefreshShareDecrypt).WriteTo(data[ptr : ptr+lenDecrypt])
	if err != nil {
		return []byte{}, err
	}

	ptr += tmp
	tmp, err = (*share.RefreshShareRecrypt).WriteTo(data[ptr : ptr+lenRecrypt])
	if err != nil {
		return []byte{}, err
	}

	return data, nil
}

// UnmarshalBinary decodes a marshaled RefreshShare on the target RefreshShare.
func (share *RefreshShare) UnmarshalBinary(data []byte) error {
	lenDecrypt := binary.BigEndian.Uint64(data[0:8])
	lenRecrypt := binary.BigEndian.Uint64(data[8:16])
	ptr := uint64(16)
	if share.RefreshShareRecrypt == nil || share.RefreshShareDecrypt == nil {
		share.RefreshShareRecrypt = new(ring.Poly)
		share.RefreshShareDecrypt = new(ring.Poly)

	}

	err := (*share.RefreshShareDecrypt).UnmarshalBinary(data[ptr : ptr+lenDecrypt])
	if err != nil {
		return err
	}
	ptr += lenDecrypt
	err = (*share.RefreshShareRecrypt).UnmarshalBinary(data[ptr : ptr+lenRecrypt])
	if err != nil {
		return err
	}
	return nil
}

// NewRefreshProtocol creates a new Refresh protocol instance.
func NewRefreshProtocol(params *bfv.Parameters) (refreshProtocol *RefreshProtocol) {
	context := newRefreshProtocolContext(params)

	refreshProtocol = new(RefreshProtocol)
	refreshProtocol.context = context
	refreshProtocol.tmp1 = context.contextKeys.NewPoly()
	refreshProtocol.tmp2 = context.contextKeys.NewPoly()
	refreshProtocol.hP = context.contextPKeys.NewPoly()

	refreshProtocol.baseconverter = ring.NewFastBasisExtender(context.contextQ.Modulus, context.specialPrimes)

	return
}

// AllocateShares allocates the shares of the Refresh protocol.
func (rfp *RefreshProtocol) AllocateShares() RefreshShare {
	return RefreshShare{rfp.context.contextQ.NewPoly(),
		rfp.context.contextQ.NewPoly()}
}

// GenShares generates a share for the Refresh protocol.
func (rfp *RefreshProtocol) GenShares(sk *ring.Poly, ciphertext *bfv.Ciphertext, crs *ring.Poly, share RefreshShare) {

	level := uint64(len(ciphertext.Value()[1].Coeffs) - 1)

	contextQ := rfp.context.contextQ
	contextT := rfp.context.contextT
	contextKeys := rfp.context.contextKeys
	contextP := rfp.context.contextPKeys
	sampler := rfp.context.contextKeys.NewKYSampler(3.19, 19) // TODO : add smudging noise

	// h0 = s*ct[1]
	contextQ.NTT(ciphertext.Value()[1], rfp.tmp1)
	contextQ.MulCoeffsMontgomery(sk, rfp.tmp1, share.RefreshShareDecrypt)

	contextQ.InvNTT(share.RefreshShareDecrypt, share.RefreshShareDecrypt)

	// h0 = s*ct[1]*P
	contextQ.MulScalarBigint(share.RefreshShareDecrypt, contextP.ModulusBigint, share.RefreshShareDecrypt)

	// h0 = s*ct[1]*P + e
	sampler.Sample(rfp.tmp1)
	contextQ.Add(share.RefreshShareDecrypt, rfp.tmp1, share.RefreshShareDecrypt)

	for x, i := 0, uint64(len(contextQ.Modulus)); i < uint64(len(rfp.context.contextKeys.Modulus)); x, i = x+1, i+1 {
		tmphP := rfp.hP.Coeffs[x]
		tmp1 := rfp.tmp1.Coeffs[i]
		for j := uint64(0); j < contextQ.N; j++ {
			tmphP[j] += tmp1[j]
		}
	}

	// h0 = (s*ct[1]*P + e)/P
	rfp.baseconverter.ModDownSplited(contextQ, contextP, rfp.context.rescaleParamsKeys, level, share.RefreshShareDecrypt, rfp.hP, share.RefreshShareDecrypt, rfp.tmp1)

	// h1 = -s*a
	contextKeys.NTT(crs, rfp.tmp1)
	contextKeys.MulCoeffsMontgomery(sk, rfp.tmp1, rfp.tmp2)
	contextKeys.Neg(rfp.tmp2, rfp.tmp2)
	contextKeys.InvNTT(rfp.tmp2, rfp.tmp2)

	// h1 = s*a + e'
	sampler.SampleAndAdd(rfp.tmp2)

	// h1 = (-s*a + e')/P
	rfp.baseconverter.ModDown(contextKeys, rfp.context.rescaleParamsKeys, level, rfp.tmp2, share.RefreshShareRecrypt, rfp.tmp1)

	// mask = (uniform plaintext in [0, T-1]) * floor(Q/T)
	coeffs := contextT.NewUniformPoly()
	lift(coeffs, rfp.tmp1, rfp.context)

	// h0 = (s*ct[1]*P + e)/P + mask
	contextQ.Add(share.RefreshShareDecrypt, rfp.tmp1, share.RefreshShareDecrypt)

	// h1 = (-s*a + e')/P - mask
	contextQ.Sub(share.RefreshShareRecrypt, rfp.tmp1, share.RefreshShareRecrypt)
}

// Aggregate sums share1 and share2 on shareOut.
func (rfp *RefreshProtocol) Aggregate(share1, share2, shareOut RefreshShare) {
	rfp.context.contextQ.Add(share1.RefreshShareDecrypt, share2.RefreshShareDecrypt, shareOut.RefreshShareDecrypt)
	rfp.context.contextQ.Add(share1.RefreshShareRecrypt, share2.RefreshShareRecrypt, shareOut.RefreshShareRecrypt)
}

// Decrypt operates a masked decryption on the input ciphertext using the provided decryption shares.
func (rfp *RefreshProtocol) Decrypt(ciphertext *bfv.Ciphertext, shareDecrypt RefreshShareDecrypt, sharePlaintext *ring.Poly) {
	rfp.context.contextQ.Add(ciphertext.Value()[0], shareDecrypt, sharePlaintext)
}

// Recode decodes and re-encode (removing the error) the masked decrypted ciphertext.
func (rfp *RefreshProtocol) Recode(sharePlaintext *ring.Poly, sharePlaintextOut *ring.Poly) {
	scaler := ring.NewSimpleScaler(rfp.context.t, rfp.context.contextQ)

	scaler.Scale(sharePlaintext, sharePlaintextOut)
	lift(sharePlaintextOut, sharePlaintextOut, rfp.context)
}

// Recrypt recrypts the input masked decrypted ciphertext with the recryption shares.
func (rfp *RefreshProtocol) Recrypt(sharePlaintext *ring.Poly, crs *ring.Poly, shareRecrypt RefreshShareRecrypt, ciphertextOut *bfv.Ciphertext) {

	// ciphertext[0] = (-crs*s + e')/P + m
	rfp.context.contextQ.Add(sharePlaintext, shareRecrypt, ciphertextOut.Value()[0])

	// ciphertext[1] = crs/P
	rfp.baseconverter.ModDown(rfp.context.contextKeys, rfp.context.rescaleParamsKeys, uint64(len(ciphertextOut.Value()[1].Coeffs)-1), crs, ciphertextOut.Value()[1], rfp.tmp1)

}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *RefreshProtocol) Finalize(ciphertext *bfv.Ciphertext, crs *ring.Poly, share RefreshShare, ciphertextOut *bfv.Ciphertext) {
	rfp.Decrypt(ciphertext, share.RefreshShareDecrypt, rfp.tmp1)
	rfp.Recode(rfp.tmp1, rfp.tmp1)
	rfp.Recrypt(rfp.tmp1, crs, share.RefreshShareRecrypt, ciphertextOut)
}

func lift(p0, p1 *ring.Poly, context *refreshProtocolContext) {
	for j := uint64(0); j < context.n; j++ {
		for i := len(context.contextQ.Modulus) - 1; i >= 0; i-- {
			p1.Coeffs[i][j] = ring.MRed(p0.Coeffs[0][j], context.deltaMont[i], context.contextQ.Modulus[i], context.contextQ.GetMredParams()[i])
		}
	}
}
