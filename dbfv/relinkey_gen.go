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
	ternarySampler  *ring.TernarySampler
	gaussianSampler *ring.KYSampler
	bitDecomp       uint64
	bitLog          uint64
	tmpPoly1        *ring.Poly
	tmpPoly2        *ring.Poly
}

type RKGShareRoundOne struct {
	//todo maybe optimize and not have to store modulus and bitLog in the struct
	//todo uint or int ?
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
	data := make([]byte, 2*8+int(share.modulus*share.bitLog)*share.share[0][0].GetDataLen())

	//share.modulus = data[0]
	binary.LittleEndian.PutUint64(data[0:8], share.modulus)
	//share.bitLog = data[1]
	binary.LittleEndian.PutUint64(data[8:16], share.bitLog)
	//write all the polys
	ptr := 16
	for i := 0; i < int(share.modulus); i++ {
		for j := 0; j < int(share.bitLog); j++ {
			r := share.share[i][j]
			n, err := r.WriteTo(data[ptr : ptr+r.GetDataLen()])
			if err != nil {
				return []byte{}, err
			}

			ptr += int(n)
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
	data := make([]byte, 2*8+2*int(share.modulus*share.bitLog)*share.share[0][0][0].GetDataLen())

	//share.modulus = data[0]
	binary.LittleEndian.PutUint64(data[0:8], share.modulus)
	//share.bitLog = data[1]
	binary.LittleEndian.PutUint64(data[8:16], share.bitLog)

	//write all of our rings in the data.
	//write all the polys
	ptr := 16
	for i := 0; i < int(share.modulus); i++ {
		for j := 0; j < int(share.bitLog); j++ {
			r0 := share.share[i][j][0]
			r1 := share.share[i][j][1]
			//write first ring
			n, err := r0.WriteTo(data[ptr : ptr+r0.GetDataLen()])
			if err != nil {
				return []byte{}, err
			}

			ptr += int(n)
			//write second ring
			n, err = r1.WriteTo(data[ptr : ptr+r1.GetDataLen()])
			if err != nil {
				return []byte{}, err
			}

			ptr += int(n)
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
	data := make([]byte, 2*8+int(share.modulus*share.bitLog)*share.share[0][0].GetDataLen())

	//share.modulus = data[0]
	binary.LittleEndian.PutUint64(data[0:8], share.modulus)
	//share.bitLog = data[1]
	binary.LittleEndian.PutUint64(data[8:16], share.bitLog)
	//write all the polys
	ptr := 16
	for i := 0; i < int(share.modulus); i++ {
		for j := 0; j < int(share.bitLog); j++ {
			r := share.share[i][j]
			n, err := r.WriteTo(data[ptr : ptr+r.GetDataLen()])
			if err != nil {
				return []byte{}, err
			}

			ptr += int(n)
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
	r1.share = make([][]*ring.Poly, len(ekg.ringContext.Modulus))
	r2.share = make([][][2]*ring.Poly, len(ekg.ringContext.Modulus))
	r3.share = make([][]*ring.Poly, len(ekg.ringContext.Modulus))
	for i := range ekg.ringContext.Modulus {
		r1.share[i] = make([]*ring.Poly, ekg.bitLog)
		r2.share[i] = make([][2]*ring.Poly, ekg.bitLog)
		r3.share[i] = make([]*ring.Poly, ekg.bitLog)
		for w := uint64(0); w < ekg.bitLog; w++ {
			r1.share[i][w] = ekg.ringContext.NewPoly()
			r2.share[i][w][0] = ekg.ringContext.NewPoly()
			r2.share[i][w][1] = ekg.ringContext.NewPoly()
			r3.share[i][w] = ekg.ringContext.NewPoly()
		}
	}

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
func NewEkgProtocol(context *bfv.BfvContext, bitDecomp uint64) *RKGProtocol {
	ekg := new(RKGProtocol)
	ekg.ringContext = context.ContextQ()
	ekg.ternarySampler = context.TernarySampler()
	ekg.gaussianSampler = context.GaussianSampler()
	ekg.bitDecomp = bitDecomp
	ekg.bitLog = uint64(math.Ceil(float64(60) / float64(bitDecomp)))
	ekg.tmpPoly1 = ekg.ringContext.NewPoly()
	ekg.tmpPoly2 = ekg.ringContext.NewPoly()
	return ekg
}

// NewEphemeralKey generates a new Ephemeral Key u_i (needs to be stored for the 3 first round).
// Each party is required to pre-compute a secret additional ephemeral key in addition to its share
// of the collective secret-key.
func (ekg *RKGProtocol) NewEphemeralKey(p float64) (ephemeralKey *ring.Poly, err error) {
	if ephemeralKey, err = ekg.ternarySampler.SampleMontgomeryNTTNew(p); err != nil {
		return nil, err
	}
	return
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
func (ekg *RKGProtocol) GenShareRoundOne(u, sk *ring.Poly, crp [][]*ring.Poly, shareOut RKGShareRoundOne) {

	mredParams := ekg.ringContext.GetMredParams()

	// Given a base decomposition w (here the CRT decomposition)
	// computes [-u_i*a + s_i*w + e_i]
	// where a = crp
	for i, qi := range ekg.ringContext.Modulus {

		for w := uint64(0); w < ekg.bitLog; w++ {

			// h = e
			ekg.gaussianSampler.SampleNTT(shareOut.share[i][w])

			// h = sk*CrtBaseDecompQi + e
			for j := uint64(0); j < ekg.ringContext.N; j++ {
				shareOut.share[i][w].Coeffs[i][j] += ring.PowerOf2(sk.Coeffs[i][j], ekg.bitDecomp*w, qi, mredParams[i])
			}

			// h = sk*CrtBaseDecompQi + -u*a + e
			ekg.ringContext.MulCoeffsMontgomeryAndSub(u, crp[i][w], shareOut.share[i][w])
		}
	}

	return
}

func (ekg *RKGProtocol) AggregateShareRoundOne(share1, share2, shareOut RKGShareRoundOne) {

	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			ekg.ringContext.Add(share1.share[i][w], share2.share[i][w], shareOut.share[i][w])
		}
	}
}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Uppon received the j-1 shares, each party computes :
//
// [s_i * sum([-u_j*a + s_j*w + e_j]) + e_i1, s_i*a + e_i2]
//
// = [s_i * (-u*a + s*w + e) + e_i1, s_i*a + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(round1 RKGShareRoundOne, sk *ring.Poly, crp [][]*ring.Poly, shareOut RKGShareRoundTwo) {

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i :
	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {

			// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

			// (AggregateShareRoundTwo samples) * sk
			ekg.ringContext.MulCoeffsMontgomery(round1.share[i][w], sk, shareOut.share[i][w][0])

			// (AggregateShareRoundTwo samples) * sk + e_1i
			ekg.gaussianSampler.SampleNTT(ekg.tmpPoly1)
			ekg.ringContext.Add(shareOut.share[i][w][0], ekg.tmpPoly1, shareOut.share[i][w][0])

			// Second Element
			// e_2i
			ekg.gaussianSampler.SampleNTT(shareOut.share[i][w][1])
			// s*a + e_2i
			ekg.ringContext.MulCoeffsMontgomeryAndAdd(sk, crp[i][w], shareOut.share[i][w][1])
		}
	}
}

// AggregateShareRoundTwo is the first part of the third and last round of the RKGProtocol protocol. Uppon receiving the j-1 elements, each party
// computues :
//
// [sum(s_j * (-u*a + s*w + e) + e_j1), sum(s_j*a + e_j2)]
//
// = [s * (-u*a + s*w + e) + e_1, s*a + e_2].
func (ekg *RKGProtocol) AggregateShareRoundTwo(share1, share2, shareOut RKGShareRoundTwo) {

	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			ekg.ringContext.Add(share1.share[i][w][0], share2.share[i][w][0], shareOut.share[i][w][0])
			ekg.ringContext.Add(share1.share[i][w][1], share2.share[i][w][1], shareOut.share[i][w][1])
		}
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

	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			// (u - s) * (sum [x][s*a_i + e_2i]) + e3i
			ekg.gaussianSampler.SampleNTT(shareOut.share[i][w])
			ekg.ringContext.MulCoeffsMontgomeryAndAdd(ekg.tmpPoly1, round2.share[i][w][1], shareOut.share[i][w])
		}
	}
}

func (ekg *RKGProtocol) AggregateShareRoundThree(share1, share2, shareOut RKGShareRoundThree) {
	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			ekg.ringContext.Add(share1.share[i][w], share2.share[i][w], shareOut.share[i][w])
		}
	}
}

func (ekg *RKGProtocol) GenRelinearizationKey(round2 RKGShareRoundTwo, round3 RKGShareRoundThree, evalKeyOut *bfv.EvaluationKey) {

	key := evalKeyOut.Get()[0].Get()
	for i := range ekg.ringContext.Modulus {
		for w := uint64(0); w < ekg.bitLog; w++ {
			ekg.ringContext.Add(round2.share[i][w][0], round3.share[i][w], key[i][w][0])
			key[i][w][1].Copy(round2.share[i][w][1])

			ekg.ringContext.MForm(key[i][w][0], key[i][w][0])
			ekg.ringContext.MForm(key[i][w][1], key[i][w][1])
		}
	}
}
