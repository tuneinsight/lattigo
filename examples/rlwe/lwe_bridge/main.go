package main

import (
	"fmt"
	"math"
	"math/bits"
	"math/rand"
	"time"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// This example implements an oblivious shuffling of the plaintext slots of an RLWE encryption.
// In such a scenario, an evaluator performs an arbitrary permutation (in this example, a random shuffle)
// over the ciphertext's slots, homomorphically.
//
// The circuit uses the LWE <-> RLWE conversion "Efficient Homomorphic Conversion Between (Ring) LWE Ciphertexts"
// of Hao Chen and Wei Dai and Miran Kim and Yongsoo Song (https://eprint.iacr.org/2020/015) to extract each message
// slot as one LWE ciphertext, shuffles the obtained LWE ciphertexts and repacks all the LWE ciphertexts' messages
// into a single a single RLWE ciphertext.
func main() {

	LogN := 12
	LogSlots := 8 // can go up to LogN

	RLWEParams := rlwe.ParametersLiteral{
		LogN:     LogN,
		Q:        []uint64{0x80000000080001},
		P:        []uint64{0x1fffffffffe00001},
		Sigma:    rlwe.DefaultSigma,
		RingType: ring.Standard,
	}
	scale := float64(1 << 40)

	params, _ := rlwe.NewParametersFromLiteral(RLWEParams)
	ringQ := params.RingQ()
	ks := rlwe.NewKeySwitcher(params)
	kgen := rlwe.NewKeyGenerator(params)

	NTimesLogSlotsInv := ring.ModExp((1 << (2*LogN - LogSlots)), ringQ.Modulus[0]-2, ringQ.Modulus[0])
	NTimesLogSlotsInv = ring.MForm(NTimesLogSlotsInv, ringQ.Modulus[0], ringQ.BredParams[0])

	sk := kgen.GenSecretKey()
	encryptor := rlwe.NewEncryptor(params, sk)
	decryptor := rlwe.NewDecryptor(params, sk)

	// Rotation Keys
	rotations := []int{}
	for i := 1; i < params.N(); i <<= 1 {
		rotations = append(rotations, i)
	}

	rtks := kgen.GenRotationKeysForRotations(rotations, true, sk)

	permuteNTTIndex := make(map[uint64][]uint64, len(rtks.Keys))
	for galEl := range rtks.Keys {
		permuteNTTIndex[galEl] = ringQ.PermuteNTTIndex(galEl)
	}

	// Plaintext generation & Encryption
	plaintext := rlwe.NewPlaintext(params, params.MaxLevel())
	plaintext.Value.IsNTT = true
	gap := 1 << (LogN - LogSlots)
	for i, j := 0, 0; i < ringQ.N; i, j = i+gap, j+1 {
		plaintext.Value.Coeffs[0][i] = uint64(float64(j%2) * scale)
	}

	ringQ.NTT(plaintext.Value, plaintext.Value)
	ciphertext := rlwe.NewCiphertextNTT(params, 1, plaintext.Level())
	encryptor.Encrypt(plaintext, ciphertext)

	fmt.Println("Plaintext before slot-shuffling:\nM(X) =")
	DecryptAndPrint(decryptor, LogSlots, ringQ, ciphertext, plaintext, scale)
	fmt.Println()

	//RLWE to LWEs: extracts each coefficient of enc(M(X)) = RLWE into a separate LWE ciphertext
	// such that dec(RLWE)[i] = dec(LWE[i])
	now := time.Now()
	fmt.Printf("Extracting RLWE  -> LWEs")
	LWE := RLWEToLWE(ciphertext, LogSlots, params)
	fmt.Printf("  Done : %s\n", time.Since(now))

	fmt.Printf("Shuffling the LWE samples")
	now = time.Now()
	rand.Shuffle(len(LWE), func(i, j int) { // for the sake of the example, this is using an insecure, fixed-seed RNG
		lwei := LWE[i]
		LWE[i] = LWE[j]
		LWE[j] = lwei
	})
	fmt.Printf(" Done : %s\n", time.Since(now))

	// LWEs to RLWEs: switches each individual LWE ciphertext into RLWE ciphertexts such that
	// dec(RLWEs[i])[0] = dec(LWE[i])
	now = time.Now()
	fmt.Printf("Switching  LWEs   -> RLWEs")
	ciphertexts := LWEToRLWE(LWE, ciphertext.Level(), params)
	fmt.Printf(" Done : %s\n", time.Since(now))

	// RLWEs to RLWE: repacks all the RLWEs into a single RLWE such that dec(RLWE) = M(X)
	fmt.Printf("Repacking  RLWEs -> RLWE")
	XPow := make(map[int]*ring.Poly)
	for i := 1; i < LogSlots+1; i++ {
		XpowNoverL := ringQ.NewPoly()
		XpowNoverL.Coeffs[0][ringQ.N/(1<<i)] = ring.MForm(1, ringQ.Modulus[0], ringQ.BredParams[0])
		ringQ.NTT(XpowNoverL, XpowNoverL)
		XPow[i] = XpowNoverL
	}
	now = time.Now()
	ciphertext = PackLWEs(ciphertexts, ks, rtks, XPow, NTimesLogSlotsInv, permuteNTTIndex, params)

	// Trace
	tmp := rlwe.NewCiphertextNTT(params, 1, plaintext.Level())
	for i := LogSlots - 1; i < LogN-1; i++ {
		Rotate(ciphertext, params.GaloisElementForColumnRotationBy(1<<i), permuteNTTIndex, params, ks, rtks, tmp)
		ringQ.Add(ciphertext.Value[0], tmp.Value[0], ciphertext.Value[0])
		ringQ.Add(ciphertext.Value[1], tmp.Value[1], ciphertext.Value[1])
	}
	fmt.Printf("  Done : %s\n", time.Since(now))

	fmt.Println("\nPlaintext after slot-shuffling:\nM'(X) =")
	DecryptAndPrint(decryptor, LogSlots, ringQ, ciphertext, plaintext, scale)
}

// LWESample is a struct for RNS LWE samples
type LWESample struct {
	b uint64
	a []uint64
}

// RLWEToLWE extracts each slot of a R-LWE ciphertext as a LWE sample.
func RLWEToLWE(ciphertext *rlwe.Ciphertext, LogSlots int, params rlwe.Parameters) (LWE []LWESample) {

	ringQ := params.RingQ()

	lvl := ciphertext.Level()
	a := ringQ.NewPolyLvl(lvl)
	b := ringQ.NewPolyLvl(lvl)

	ringQ.InvNTT(ciphertext.Value[0], b)
	ringQ.InvNTT(ciphertext.Value[1], a)

	LWE = make([]LWESample, 1<<LogSlots)

	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	acc := ringQ.NewPoly()
	for i, qi := range ringQ.Modulus {
		tmp0 := acc.Coeffs[i]
		tmp1 := a.Coeffs[i]
		tmp0[0] = tmp1[0]
		for j := 1; j < ringQ.N; j++ {
			tmp0[j] = qi - tmp1[ringQ.N-j]
		}
	}

	pol := b

	gap := 1 << (params.LogN() - LogSlots)

	// Real values
	for i, j := 0, 0; i < ringQ.N; i, j = i+gap, j+1 {

		LWE[j].b = pol.Coeffs[0][i]
		LWE[j].a = make([]uint64, ringQ.N)
		copy(LWE[j].a, acc.Coeffs[0])

		// Multiplies the accumulator by X
		MulBySmallMonomial(ringQ, acc, gap)
	}

	return
}

// LWEToRLWE transforms a set of LWE samples into their respective RLWE ciphertext such that decrypt(ct)[0] = decrypt(LWE)
func LWEToRLWE(lwe []LWESample, level int, params rlwe.Parameters) (ciphertexts []*rlwe.Ciphertext) {
	ringQ := params.RingQ()
	acc := ringQ.NewPoly()
	gap := params.N() / len(lwe)
	ciphertexts = make([]*rlwe.Ciphertext, params.N())
	for i, j := 0, 0; i < params.N(); i, j = i+gap, j+1 {

		// Alocates ciphertext
		ciphertexts[i] = rlwe.NewCiphertextNTT(params, 1, level)
		ciphertexts[i].Value[0].Coeffs[0][0] = lwe[j].b

		// Copies coefficients multiplied by X^{N-1} in reverse order:
		// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
		tmp0, tmp1 := acc.Coeffs[0], lwe[j].a
		tmp0[0] = tmp1[0]
		for k := 1; k < ringQ.N; k++ {
			tmp0[k] = ringQ.Modulus[0] - tmp1[ringQ.N-k]
		}

		copy(ciphertexts[i].Value[1].Coeffs[0], acc.Coeffs[0])

		// Switches to NTT domain
		ringQ.NTT(ciphertexts[i].Value[0], ciphertexts[i].Value[0])
		ringQ.NTT(ciphertexts[i].Value[1], ciphertexts[i].Value[1])
	}
	return
}

// PackLWEs repacks LWE ciphertexts into a RLWE ciphertext
func PackLWEs(ciphertexts []*rlwe.Ciphertext, ks *rlwe.KeySwitcher, rtks *rlwe.RotationKeySet, XPow map[int]*ring.Poly, NTimesLogSlotsInv uint64, permuteNTTIndex map[uint64][]uint64, params rlwe.Parameters) *rlwe.Ciphertext {

	ringQ := params.RingQ()

	L := bits.Len64(uint64(len(ciphertexts))) - 1

	if L == 0 {
		// Multiplies by (Slots * N) ^-1 mod Q
		if ciphertexts[0] != nil {
			ring.MulScalarMontgomeryVec(ciphertexts[0].Value[0].Coeffs[0], ciphertexts[0].Value[0].Coeffs[0], NTimesLogSlotsInv, ringQ.Modulus[0], ringQ.MredParams[0])
			ring.MulScalarMontgomeryVec(ciphertexts[0].Value[1].Coeffs[0], ciphertexts[0].Value[1].Coeffs[0], NTimesLogSlotsInv, ringQ.Modulus[0], ringQ.MredParams[0])
		}
		return ciphertexts[0]
	}

	odd := make([]*rlwe.Ciphertext, len(ciphertexts)>>1)
	even := make([]*rlwe.Ciphertext, len(ciphertexts)>>1)

	for i := 0; i < len(ciphertexts)>>1; i++ {
		odd[i] = ciphertexts[2*i]
		even[i] = ciphertexts[2*i+1]
	}

	ctEven := PackLWEs(odd, ks, rtks, XPow, NTimesLogSlotsInv, permuteNTTIndex, params)
	ctOdd := PackLWEs(even, ks, rtks, XPow, NTimesLogSlotsInv, permuteNTTIndex, params)

	if ctEven == nil && ctOdd == nil {
		return nil
	}

	var tmpEven *rlwe.Ciphertext
	if ctEven != nil {
		tmpEven = ctEven.CopyNew()
	}

	// ctOdd * X^(N/2^L)
	if ctOdd != nil {
		//X^(N/2^L)
		ringQ.MulCoeffsMontgomery(ctOdd.Value[0], XPow[L], ctOdd.Value[0])
		ringQ.MulCoeffsMontgomery(ctOdd.Value[1], XPow[L], ctOdd.Value[1])

		// ctEven + ctOdd * X^(N/2^L)
		ringQ.Add(ctEven.Value[0], ctOdd.Value[0], ctEven.Value[0])
		ringQ.Add(ctEven.Value[1], ctOdd.Value[1], ctEven.Value[1])

		// phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.Sub(tmpEven.Value[0], ctOdd.Value[0], tmpEven.Value[0])
		ringQ.Sub(tmpEven.Value[1], ctOdd.Value[1], tmpEven.Value[1])
	}

	if ctEven != nil {

		// if L-2 == -1, then gal = 2N-1
		if L == 1 {
			Rotate(tmpEven, uint64(2*ringQ.N-1), permuteNTTIndex, params, ks, rtks, tmpEven)
		} else {
			Rotate(tmpEven, params.GaloisElementForColumnRotationBy(1<<(L-2)), permuteNTTIndex, params, ks, rtks, tmpEven)
		}

		// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.Add(ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
		ringQ.Add(ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])
	}

	return ctEven
}

//MulBySmallMonomial multiplies pol by x^n
func MulBySmallMonomial(ringQ *ring.Ring, pol *ring.Poly, n int) {
	for i, qi := range ringQ.Modulus[:pol.Level()+1] {
		pol.Coeffs[i] = append(pol.Coeffs[i][ringQ.N-n:], pol.Coeffs[i][:ringQ.N-n]...)
		tmp := pol.Coeffs[i]
		for j := 0; j < n; j++ {
			tmp[j] = qi - tmp[j]
		}
	}
}

// Rotate rotates a ciphertext
func Rotate(ctIn *rlwe.Ciphertext, galEl uint64, permuteNTTindex map[uint64][]uint64, params rlwe.Parameters, ks *rlwe.KeySwitcher, rtks *rlwe.RotationKeySet, ctOut *rlwe.Ciphertext) {
	ringQ := params.RingQ()
	rtk, _ := rtks.GetRotationKey(galEl)
	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	index := permuteNTTindex[galEl]
	ks.SwitchKeysInPlace(level, ctIn.Value[1], rtk, ks.Pool[1].Q, ks.Pool[2].Q)
	ringQ.AddLvl(level, ks.Pool[1].Q, ctIn.Value[0], ks.Pool[1].Q)
	ringQ.PermuteNTTWithIndexLvl(level, ks.Pool[1].Q, index, ctOut.Value[0])
	ringQ.PermuteNTTWithIndexLvl(level, ks.Pool[2].Q, index, ctOut.Value[1])
}

// DecryptAndPrint decrypts and prints the first N values.
func DecryptAndPrint(decryptor rlwe.Decryptor, LogSlots int, ringQ *ring.Ring, ciphertext *rlwe.Ciphertext, plaintext *rlwe.Plaintext, scale float64) {
	decryptor.Decrypt(ciphertext, plaintext)
	ringQ.InvNTT(plaintext.Value, plaintext.Value)

	v := make([]float64, 1<<LogSlots)

	gap := ringQ.N / (1 << LogSlots)
	for i, j := 0, 0; i < 1<<LogSlots; i, j = i+1, j+gap {
		if plaintext.Value.Coeffs[0][j] >= ringQ.Modulus[0]>>1 {
			v[i] = -float64(ringQ.Modulus[0] - plaintext.Value.Coeffs[0][j])
		} else {
			v[i] = float64(plaintext.Value.Coeffs[0][j])
		}

		v[i] /= scale
	}

	fmt.Printf("   ")
	for i := 0; i < 1<<LogSlots; i++ {
		if i != 0 && i&31 == 0 {
			fmt.Printf("\n   ")
		}
		fmt.Printf("%d ", uint64(math.Round(v[i])))

	}
	fmt.Printf("\n")
}
