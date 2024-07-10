package sign

import (
	"bytes"
	"encoding/binary"
	"log"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
	"github.com/zeebo/blake3"
)

const keySize = 32

// PRNGKey generates a key for PRNG using the secret key share
func PRNGKey(skShare structs.Vector[ring.Poly]) []byte {
	hasher := blake3.New()
	buf := new(bytes.Buffer)
	skShare.WriteTo(buf)
	hasher.Write(buf.Bytes())

	skHash := hasher.Sum(nil)
	return skHash[:keySize]
}

// GenerateMAC generates a MAC for a given TildeD matrix and mask
func GenerateMAC(TildeD structs.Matrix[ring.Poly], MACKey []byte, partyID int, sid int, T []int, otherParty int, verify bool) []byte {
	hasher := blake3.New()
	buf := new(bytes.Buffer)

	if verify {
		binary.Write(buf, binary.BigEndian, int64(otherParty))
	} else {
		binary.Write(buf, binary.BigEndian, int64(partyID))
	}

	binary.Write(buf, binary.BigEndian, MACKey)
	TildeD.WriteTo(buf)
	binary.Write(buf, binary.BigEndian, int64(sid))
	binary.Write(buf, binary.BigEndian, T)

	hasher.Write(buf.Bytes())
	MAC := hasher.Sum(nil)
	return MAC[:keySize]
}

// Hashes parameters to a Gaussian distribution
func GaussianHash(r *ring.Ring, hash []byte, mu string, sigmaU float64, boundU float64, length int) structs.Vector[ring.Poly] {
	hasher := blake3.New()
	buf := new(bytes.Buffer)

	binary.Write(buf, binary.BigEndian, hash)
	buf.WriteString(mu)

	hasher.Write(buf.Bytes())
	hashOutput := hasher.Sum(nil)

	prng, _ := sampling.NewKeyedPRNG(hashOutput[:keySize])
	gaussianParams := ring.DiscreteGaussian{Sigma: sigmaU, Bound: boundU}
	hashGaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	return SamplePolyVector(r, length, hashGaussianSampler, true, true)
}

// PRF generates pseudorandom ring elements
func PRF(r *ring.Ring, sd_ij []byte, PRFKey []byte, mu string, hash []byte, n int) structs.Vector[ring.Poly] {
	hasher := blake3.New()
	buf := new(bytes.Buffer)

	binary.Write(buf, binary.BigEndian, PRFKey)
	binary.Write(buf, binary.BigEndian, sd_ij)
	binary.Write(buf, binary.BigEndian, hash)
	buf.WriteString(mu)

	hasher.Write(buf.Bytes())
	hashOutput := hasher.Sum(nil)

	prng, _ := sampling.NewKeyedPRNG(hashOutput[:keySize])
	PRFUniformSampler := ring.NewUniformSampler(prng, r)
	mask := SamplePolyVector(r, n, PRFUniformSampler, true, true)
	return mask
}

// Hashes precomputable values
func Hash(A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], D map[int]structs.Matrix[ring.Poly], sid int, T []int) []byte {
	hasher := blake3.New()
	buf := new(bytes.Buffer)

	if _, err := A.WriteTo(buf); err != nil {
		log.Fatalf("Error writing matrix A: %v\n", err)
	}

	if _, err := b.WriteTo(buf); err != nil {
		log.Fatalf("Error writing vector b: %v\n", err)
	}

	binary.Write(buf, binary.BigEndian, int64(sid))
	binary.Write(buf, binary.BigEndian, T)

	for i := 0; i < len(D); i++ {
		if _, err := D[i].WriteTo(buf); err != nil {
			log.Fatalf("Error writing matrix D_i: %v\n", err)
		}
	}

	hasher.Write(buf.Bytes())
	hashOutput := hasher.Sum(nil)
	return hashOutput[:keySize]
}

// Hashes to low norm ring elements
func LowNormHash(r *ring.Ring, A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], h structs.Vector[ring.Poly], mu string, kappa int) ring.Poly {
	hasher := blake3.New()
	buf := new(bytes.Buffer)

	if _, err := A.WriteTo(buf); err != nil {
		log.Fatalf("Error writing matrix A: %v\n", err)
	}

	if _, err := b.WriteTo(buf); err != nil {
		log.Fatalf("Error writing vector b: %v\n", err)
	}

	if _, err := h.WriteTo(buf); err != nil {
		log.Fatalf("Error writing vector h: %v\n", err)
	}

	binary.Write(buf, binary.BigEndian, []byte(mu))

	hasher.Write(buf.Bytes())
	hashOutput := hasher.Sum(nil)

	prng, _ := sampling.NewKeyedPRNG(hashOutput[:keySize])
	ternaryParams := ring.Ternary{H: kappa}
	ternarySampler, err := ring.NewTernarySampler(prng, r, ternaryParams, false)
	if err != nil {
		log.Fatalf("Error creating ternary sampler: %v", err)
	}
	c := ternarySampler.ReadNew()
	r.NTT(c, c)
	r.MForm(c, c)

	return c
}

// GenerateRandomSeed generates a random seed of length ell
func GenerateRandomSeed() []byte {
	return GetRandomBytes(keySize)
}
