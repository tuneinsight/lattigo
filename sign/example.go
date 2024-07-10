package sign

import (
	"fmt"
	"log"
	"math/big"
	"time"

	"github.com/montanaflynn/stats"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
)

var K int
var Threshold int

// main function orchestrates the threshold signature protocol
func LocalRun(x int) {
	var totalGenDuration, totalFinalizeDuration, totalVerifyDuration time.Duration

	// Create maps to collect durations across all runs
	totalSignRound1Durations := make(map[int]float64)
	totalSignRound2PreprocessDurations := make(map[int]float64)
	totalSignRound2Durations := make(map[int]float64)

	for run := 0; run < x; run++ {
		log.Println("RUN:", run)
		var genDuration, finalizeDuration, verifyDuration time.Duration

		randomKey := make([]byte, KeySize)

		r, _ := ring.NewRing(1<<LogN, []uint64{Q})
		r_xi, _ := ring.NewRing(1<<LogN, []uint64{QXi})
		r_nu, _ := ring.NewRing(1<<LogN, []uint64{QNu})

		prng, _ := sampling.NewKeyedPRNG(randomKey)
		uniformSampler := ring.NewUniformSampler(prng, r)
		trustedDealerKey := randomKey

		parties := make([]*Party, K)
		for i := range parties {
			prng, _ := sampling.NewKeyedPRNG(randomKey)
			uniformSampler := ring.NewUniformSampler(prng, r)
			parties[i] = NewParty(i, r, r_xi, r_nu, uniformSampler)
		}

		// GEN: Generate secret shares, seeds, and MAC keys
		T := make([]int, K) // Active parties
		for i := 0; i < K; i++ {
			T[i] = i
		}
		lagrangeCoeffs := ComputeLagrangeCoefficients(r, T, big.NewInt(int64(Q)))
		log.Println("Gen")

		start := time.Now()
		A, skShares, seeds, MACKeys, b := Gen(r, r_xi, uniformSampler, trustedDealerKey, lagrangeCoeffs)
		genDuration = time.Since(start)
		log.Println("Gen Duration:", genDuration)
		for partyID := 0; partyID < K; partyID++ {
			parties[partyID].SkShare = skShares[partyID]
			parties[partyID].Seed = seeds
			parties[partyID].MACKeys = MACKeys[partyID]
		}

		// Create maps to collect durations for this run
		signRound1Durations := make(map[int]time.Duration)
		signRound2PreprocessDurations := make(map[int]time.Duration)
		signRound2Durations := make(map[int]time.Duration)

		// SIGNATURE ROUND 1
		mu := "Message"
		sid := 1
		PRFKey := GenerateRandomSeed()
		log.Println("Generating lagrange coefficients...")
		D := make(map[int]structs.Matrix[ring.Poly])
		MACs := make(map[int]map[int][]byte)
		for _, partyID := range T {
			r.NTT(lagrangeCoeffs[partyID], lagrangeCoeffs[partyID])
			r.MForm(lagrangeCoeffs[partyID], lagrangeCoeffs[partyID])
			parties[partyID].Lambda = lagrangeCoeffs[partyID]
			parties[partyID].Seed = seeds
			log.Println("Sign Round 1, party", partyID)
			start = time.Now()
			D[partyID], MACs[partyID] = parties[partyID].SignRound1(A, sid, []byte(PRFKey), T)
			signRound1Durations[partyID] = time.Since(start)
		}

		// SIGNATURE ROUND 2
		z := make(map[int]structs.Vector[ring.Poly])

		for _, partyID := range T {
			log.Println("Sign Round 2 preprocess, party", partyID)
			start = time.Now()
			valid, DSum, hash := parties[partyID].SignRound2Preprocess(A, b, D, MACs, sid, T)
			if !valid {
				log.Fatalf("MAC verification failed for party %d", partyID)
			}
			signRound2PreprocessDurations[partyID] = time.Since(start)

			log.Println("Sign round 2 party", partyID)
			start = time.Now()
			z[partyID] = parties[partyID].SignRound2(A, b, DSum, sid, mu, T, []byte(PRFKey), hash)
			signRound2Durations[partyID] = time.Since(start)
		}

		// SIGNATURE FINALIZE
		log.Println("finalizing...")
		finalParty := parties[0]
		start = time.Now()
		_, sig, Delta := finalParty.SignFinalize(z, A, b)
		finalizeDuration = time.Since(start)

		// Verify the signature
		start = time.Now()
		valid := Verify(r, r_xi, r_nu, sig, A, mu, b, finalParty.C, Delta)
		verifyDuration = time.Since(start)
		fmt.Printf("Signature Verification Result: %v\n", valid)

		// Accumulate durations
		totalGenDuration += genDuration
		totalFinalizeDuration += finalizeDuration
		totalVerifyDuration += verifyDuration

		// Accumulate phase durations
		for partyID, duration := range signRound1Durations {
			totalSignRound1Durations[partyID] += float64(duration.Nanoseconds()) / 1e6
		}
		for partyID, duration := range signRound2PreprocessDurations {
			totalSignRound2PreprocessDurations[partyID] += float64(duration.Nanoseconds()) / 1e6
		}
		for partyID, duration := range signRound2Durations {
			totalSignRound2Durations[partyID] += float64(duration.Nanoseconds()) / 1e6
		}
	}

	// Print averaged durations
	fmt.Println("Averaged durations over", x, "runs:")
	fmt.Printf("Gen duration: %.3f ms\n", float64(totalGenDuration.Nanoseconds())/1e6/float64(x))
	fmt.Printf("Finalize duration: %.3f ms\n", float64(totalFinalizeDuration.Nanoseconds())/1e6/float64(x))
	fmt.Printf("Verify duration: %.3f ms\n", float64(totalVerifyDuration.Nanoseconds())/1e6/float64(x))

	// Calculate and print averaged statistics for each phase
	printAveragedStats("Signature Round 1", totalSignRound1Durations, x)
	printAveragedStats("Signature Round 2 Preprocess", totalSignRound2PreprocessDurations, x)
	printAveragedStats("Signature Round 2", totalSignRound2Durations, x)

	// Calculate and print total signing and offline durations
	totalSigningDurations := make(map[int]float64)
	for partyID := range totalSignRound1Durations {
		totalSigningDurations[partyID] = totalSignRound1Durations[partyID] + totalSignRound2PreprocessDurations[partyID] + totalSignRound2Durations[partyID]
	}
	printAveragedStats("Total Signing", totalSigningDurations, x)
}

// printAveragedStats prints the mean, median, and standard deviation for a map of durations averaged over x runs
func printAveragedStats(phaseName string, totalDurations map[int]float64, x int) {
	var values []float64
	for _, totalDuration := range totalDurations {
		values = append(values, totalDuration/float64(x))
	}
	mean, _ := stats.Mean(values)
	median, _ := stats.Median(values)
	stddev, _ := stats.StandardDeviation(values)

	fmt.Printf("%s averaged duration stats over %d runs:\n", phaseName, x)
	fmt.Printf("  Mean: %.3f ms\n", mean)
	fmt.Printf("  Median: %.3f ms\n", median)
	fmt.Printf("  Standard Deviation: %.3f ms\n", stddev)
}
