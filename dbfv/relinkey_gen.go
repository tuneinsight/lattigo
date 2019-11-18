package dbfv

import (
	"encoding/binary"
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"math"
)

// RKGProtocol is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type RKGProtocol struct {
	ringContext     *ring.Context
	bfvContext      *bfv.BfvContext
	keyswitchprimes []uint64
	alpha           uint64
	beta            uint64
	gaussianSampler *ring.KYSampler
	tmpPoly1        *ring.Poly
	tmpPoly2        *ring.Poly
	polypool        *ring.Poly
}

//todo should be
//type RKGShareRoundOne []*ring.Poly
//type RKGShareRoundTwo [][2]*ring.Poly
//type RKGShareRoundThree []*ring.Poly
type RKGShareRoundOne struct {
	//todo maybe optimize and not have to store modulus and bitLog in the struct
	modulus uint64
	bitLog  uint64
	share   [][]*ring.Poly
}
type RKGShareRoundTwo struct {
	modulus uint64
	bitLog  uint64
	share   [][][2]*ring.Poly
}
type RKGShareRoundThree struct {
	modulus uint64
	bitLog  uint64
	share   [][]*ring.Poly
}

func (share *RKGShareRoundOne) MarshalBinary() ([]byte, error) {
	//todo ask if Modulus , bitLog should be written on 1 byte or more.

	//we have modulus * bitLog * Len of 1 ring rings
	data := make([]byte, 2*8+(share.modulus*share.bitLog)*share.share[0][0].GetDataLen(true))

	//share.modulus = data[0]
	binary.LittleEndian.PutUint64(data[0:8], share.modulus)
	//share.bitLog = data[1]
	binary.LittleEndian.PutUint64(data[8:16], share.bitLog)
	//write all the polys
	ptr := uint64(16)
	for i := 0; i < int(share.modulus); i++ {
		for j := 0; j < int(share.bitLog); j++ {
			r := share.share[i][j]
			n, err := r.WriteTo(data[ptr : ptr+r.GetDataLen(true)])
			if err != nil {
				return []byte{}, err
			}

			ptr += n
		}
	}

	return data, nil
}

func (share *RKGShareRoundOne) UnmarshalBinary(data []byte) error {
	//share.modulus = data[0]
	share.modulus = (binary.LittleEndian.Uint64(data[0:8]))
	//share.bitLog = data[1]
	share.bitLog = (binary.LittleEndian.Uint64(data[8:16]))
	share.share = make([][]*ring.Poly, share.modulus)
	ptr := 16
	length := (len(data) - 16) / int(share.modulus) / int(share.bitLog)
	for i := 0; i < int(share.modulus); i++ {

		share.share[i] = make([]*ring.Poly, share.bitLog)
		for j := 0; j < int(share.bitLog); j++ {
			//put the poly
			share.share[i][j] = new(ring.Poly)
			err := share.share[i][j].UnmarshalBinary(data[ptr : ptr+length])
			if err != nil {
				return err
			}

			ptr += length
		}

	}

	return nil
}

func (share *RKGShareRoundTwo) MarshalBinary() ([]byte, error) {
	//we have modulus * bitLog * Len of 1 ring rings
	data := make([]byte, 2*8+2*(share.modulus*share.bitLog)*share.share[0][0][0].GetDataLen(true))

	//share.modulus = data[0]
	binary.LittleEndian.PutUint64(data[0:8], share.modulus)
	//share.bitLog = data[1]
	binary.LittleEndian.PutUint64(data[8:16], share.bitLog)

	//write all of our rings in the data.
	//write all the polys
	ptr := uint64(16)
	for i := 0; i < int(share.modulus); i++ {
		for j := 0; j < int(share.bitLog); j++ {
			r0 := share.share[i][j][0]
			r1 := share.share[i][j][1]
			//write first ring
			n, err := r0.WriteTo(data[ptr : ptr+r0.GetDataLen(true)])
			if err != nil {
				return []byte{}, err
			}

			ptr += (n)
			//write second ring
			n, err = r1.WriteTo(data[ptr : ptr+r1.GetDataLen(true)])
			if err != nil {
				return []byte{}, err
			}

			ptr += (n)
		}
	}

	return data, nil
}

func (share *RKGShareRoundTwo) UnmarshalBinary(data []byte) error {
	share.modulus = (binary.LittleEndian.Uint64(data[0:8]))
	//share.bitLog = data[1]
	share.bitLog = (binary.LittleEndian.Uint64(data[8:16]))
	share.share = make([][][2]*ring.Poly, share.modulus)
	ptr := 16
	//lenght of a single ring
	length := (len(data) - 16) / int(share.modulus) / int(share.bitLog) / 2
	//now retrieve all the rings.

	for i := 0; i < int(share.modulus); i++ {

		share.share[i] = make([][2]*ring.Poly, share.bitLog)
		for j := 0; j < int(share.bitLog); j++ {
			//put the poly
			share.share[i][j][0] = new(ring.Poly)
			share.share[i][j][1] = new(ring.Poly)
			err := share.share[i][j][0].UnmarshalBinary(data[ptr : ptr+length])

			if err != nil {
				return err
			}

			ptr += length

			//read second ring
			err = share.share[i][j][1].UnmarshalBinary(data[ptr : ptr+length])

			if err != nil {
				return err
			}

			ptr += length
		}

	}

	return nil
}

func (share *RKGShareRoundThree) MarshalBinary() ([]byte, error) {
	//we have modulus * bitLog * Len of 1 ring rings
	data := make([]byte, 2*8+(share.modulus*share.bitLog)*share.share[0][0].GetDataLen(true))

	//share.modulus = data[0]
	binary.LittleEndian.PutUint64(data[0:8], share.modulus)
	//share.bitLog = data[1]
	binary.LittleEndian.PutUint64(data[8:16], share.bitLog)
	//write all the polys
	ptr := uint64(16)
	for i := 0; i < int(share.modulus); i++ {
		for j := 0; j < int(share.bitLog); j++ {
			r := share.share[i][j]
			n, err := r.WriteTo(data[ptr : ptr+r.GetDataLen(true)])
			if err != nil {
				return []byte{}, err
			}

			ptr += (n)
		}
	}

	return data, nil
}

func (share *RKGShareRoundThree) UnmarshalBinary(data []byte) error {
	//share.modulus = data[0]
	share.modulus = (binary.LittleEndian.Uint64(data[0:8]))
	//share.bitLog = data[1]
	share.bitLog = (binary.LittleEndian.Uint64(data[8:16]))
	share.share = make([][]*ring.Poly, share.modulus)
	ptr := 16
	length := (len(data) - 16) / int(share.modulus) / int(share.bitLog)
	for i := 0; i < int(share.modulus); i++ {

		share.share[i] = make([]*ring.Poly, share.bitLog)
		for j := 0; j < int(share.bitLog); j++ {
			//put the poly
			share.share[i][j] = new(ring.Poly)
			err := share.share[i][j].UnmarshalBinary(data[ptr : ptr+length])
			if err != nil {
				return err
			}

			ptr += length
		}

	}

	return nil
}

func (ekg *RKGProtocol) AllocateShares() (r1 RKGShareRoundOne, r2 RKGShareRoundTwo, r3 RKGShareRoundThree) {
	r1 = make([]*ring.Poly, ekg.beta)
	r2 = make([][2]*ring.Poly, ekg.beta)
	r3 = make([]*ring.Poly, ekg.beta)
	for i := uint64(0); i < ekg.beta; i++ {
		r1[i] = ekg.ringContext.NewPoly()
		r2[i][0] = ekg.ringContext.NewPoly()
		r2[i][1] = ekg.ringContext.NewPoly()
		r3[i] = ekg.ringContext.NewPoly()
	}
	//todo find how to rm this
	r1.bitLog = ekg.bitLog
	r1.modulus = uint64(len(ekg.ringContext.Modulus))
	r2.bitLog = ekg.bitLog
	r2.modulus = uint64(len(ekg.ringContext.Modulus))
	r3.bitLog = ekg.bitLog
	r3.modulus = uint64(len(ekg.ringContext.Modulus))

	return
}

// NewEkgProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewEkgProtocol(context *bfv.BfvContext) *RKGProtocol {

	ekg := new(RKGProtocol)
	ekg.ringContext = context.ContextKeys()
	ekg.bfvContext = context

	ekg.keyswitchprimes = make([]uint64, len(context.KeySwitchPrimes()))
	for i, pi := range context.KeySwitchPrimes() {
		ekg.keyswitchprimes[i] = pi
	}

	ekg.alpha = uint64(len(ekg.keyswitchprimes))
	ekg.beta = uint64(math.Ceil(float64(len(ekg.ringContext.Modulus)-len(ekg.keyswitchprimes)) / float64(ekg.alpha)))

	ekg.gaussianSampler = context.GaussianSampler()

	ekg.tmpPoly1 = ekg.ringContext.NewPoly()
	ekg.tmpPoly2 = ekg.ringContext.NewPoly()
	ekg.polypool = ekg.ringContext.NewPoly()

	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first round).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *RKGProtocol) NewEphemeralKey(p float64) (ephemeralKey *ring.Poly, err error) {
	if ephemeralKey, err = ekg.ringContext.SampleTernaryMontgomeryNTTNew(p); err != nil {
		return nil, err
	}
	return
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(u, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShareRoundOne) {

	var index uint64

	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i]
	// where a_i = crp_i

	ekg.polypool.Copy(sk)

	ekg.ringContext.MulScalarBigint(ekg.polypool, ekg.bfvContext.ContextPKeys().ModulusBigint, ekg.polypool)

	ekg.ringContext.InvMForm(ekg.polypool, ekg.polypool)

	for i := uint64(0); i < ekg.beta; i++ {

		// h = e
		ekg.gaussianSampler.SampleNTT(shareOut[i])

		// h = sk*CrtBaseDecompQi + e
		for j := uint64(0); j < ekg.alpha; j++ {

			index = i*ekg.alpha + j
			qi := ekg.ringContext.Modulus[index]
			tmp0 := ekg.polypool.Coeffs[index]
			tmp1 := shareOut[i].Coeffs[index]

			for w := uint64(0); w < ekg.ringContext.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index == uint64(len(ekg.ringContext.Modulus)-len(ekg.keyswitchprimes)-1) {
				break
			}
		}

		// h = sk*CrtBaseDecompQi + -u*a + e
		ekg.ringContext.MulCoeffsMontgomeryAndSub(u, crp[i], shareOut[i])
	}

	ekg.polypool.Zero()

	return
}

func (ekg *RKGProtocol) AggregateShareRoundOne(share1, share2, shareOut RKGShareRoundOne) {

	for i := uint64(0); i < ekg.beta; i++ {
		ekg.ringContext.Add(share1[i], share2[i], shareOut[i])
	}

}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Uppon received the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(round1 RKGShareRoundOne, sk *ring.Poly, crp []*ring.Poly, shareOut RKGShareRoundTwo) {

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := uint64(0); i < ekg.beta; i++ {

		// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

		// (AggregateShareRoundTwo samples) * sk
		ekg.ringContext.MulCoeffsMontgomery(round1[i], sk, shareOut[i][0])

		// (AggregateShareRoundTwo samples) * sk + e_1i
		ekg.gaussianSampler.SampleNTT(ekg.tmpPoly1)
		ekg.ringContext.Add(shareOut[i][0], ekg.tmpPoly1, shareOut[i][0])

		// Second Element
		// e_2i
		ekg.gaussianSampler.SampleNTT(shareOut[i][1])
		// s*a + e_2i
		ekg.ringContext.MulCoeffsMontgomeryAndAdd(sk, crp[i], shareOut[i][1])
	}

}

// AggregateShareRoundTwo is the first part of the third and last round of the RKGProtocol protocol. Uppon receiving the j-1 elements, each party
// computues :
//
// [sum(s_j * (-u*a + s*w + e) + e_j1), sum(s_j*a + e_j2)]
//
// = [s * (-u*a + s*w + e) + e_1, s*a + e_2].
func (ekg *RKGProtocol) AggregateShareRoundTwo(share1, share2, shareOut RKGShareRoundTwo) {

	for i := uint64(0); i < ekg.beta; i++ {
		ekg.ringContext.Add(share1[i][0], share2[i][0], shareOut[i][0])
		ekg.ringContext.Add(share1[i][1], share2[i][1], shareOut[i][1])
	}

}

// GenShareRoundThree is the second pard of the third and last round of the RKGProtocol protocol. Each party operates a key-switch on [s*a + e_2],
// by computing :
//
// [(u_i - s_i)*(s*a + e_2)]
//
// and broadcasts the result the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundThree(round2 RKGShareRoundTwo, u, sk *ring.Poly, shareOut RKGShareRoundThree) {

	// (u_i - s_i)
	ekg.ringContext.Sub(u, sk, ekg.tmpPoly1)

	for i := uint64(0); i < ekg.beta; i++ {

		// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
		ekg.gaussianSampler.SampleNTT(shareOut[i])
		ekg.ringContext.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1, round2[i][1], shareOut[i])
	}
}

func (ekg *RKGProtocol) AggregateShareRoundThree(share1, share2, shareOut RKGShareRoundThree) {
	for i := uint64(0); i < ekg.beta; i++ {
		ekg.ringContext.Add(share1[i], share2[i], shareOut[i])
	}
}

func (ekg *RKGProtocol) GenRelinearizationKey(round2 RKGShareRoundTwo, round3 RKGShareRoundThree, evalKeyOut *bfv.EvaluationKey) {

	key := evalKeyOut.Get()[0].Get()
	for i := uint64(0); i < ekg.beta; i++ {

		ekg.ringContext.Add(round2[i][0], round3[i], key[i][0])
		key[i][1].Copy(round2[i][1])

		ekg.ringContext.MForm(key[i][0], key[i][0])
		ekg.ringContext.MForm(key[i][1], key[i][1])

	}
}

//todo find how to do it otherwise
func (ekg *RKGProtocol) AllocateEvaluationKey(ctx bfv.BfvContext) *bfv.EvaluationKey {
	return ctx.NewRelinKey(uint64(len(ekg.ringContext.Modulus)), ekg.bitDecomp)
}
