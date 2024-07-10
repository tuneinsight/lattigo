package sign

import (
	"crypto/rand"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
)

// ShamirSecretSharing shares each coefficient of a vector of ring.Poly across k parties using (t, k)-threshold Shamir secret sharing.
func ShamirSecretSharingGeneral(r *ring.Ring, s []ring.Poly, t, k int) map[int]structs.Vector[ring.Poly] {

	degree := r.N() // Number of coefficients in each ring.Poly
	q := r.Modulus()

	// Initialize shares for each party
	shares := make(map[int]structs.Vector[ring.Poly], k)

	for i := 0; i < k; i++ {
		shares[i] = make([]ring.Poly, len(s))
		for j := range shares[i] {
			shares[i][j] = r.NewPoly()
		}
	}

	for polyIndex, poly := range s {
		coeffs := make([]*big.Int, degree)
		r.PolyToBigint(poly, 1, coeffs)

		for coeffIndex, secret := range coeffs {
			polyCoeffs := make([]*big.Int, t)
			polyCoeffs[0] = secret
			for i := 1; i < t; i++ {
				randomCoeff, _ := rand.Int(rand.Reader, q)
				polyCoeffs[i] = randomCoeff
			}

			for i := 1; i <= k; i++ {
				x := big.NewInt(int64(i))
				shareValue := big.NewInt(0)
				xPow := big.NewInt(1)

				for _, coeff := range polyCoeffs {
					term := new(big.Int).Mul(coeff, xPow)
					shareValue.Add(shareValue, term)
					shareValue.Mod(shareValue, q)
					xPow.Mul(xPow, x)
				}

				if shares[i-1][polyIndex].Coeffs[0] == nil {
					shares[i-1][polyIndex].Coeffs[0] = make([]uint64, degree)
				}

				shares[i-1][polyIndex].Coeffs[0][coeffIndex] = shareValue.Uint64()
			}
		}
	}

	return shares
}

// ShamirSecretSharing shares each coefficient of a vector of ring.Poly across k parties using (t, k)-threshold Shamir secret sharing. This optimized implementation only works when t = k.
func ShamirSecretSharing(r *ring.Ring, s []ring.Poly, k int, lambdas []ring.Poly) map[int]structs.Vector[ring.Poly] {

	degree := r.N() // Number of coefficients in each ring.Poly
	q := r.Modulus()

	// Initialize shares for each party
	shares := make(map[int]structs.Vector[ring.Poly], k)

	for i := 0; i < k; i++ {
		shares[i] = make([]ring.Poly, len(s))
		for j := range shares[i] {
			shares[i][j] = r.NewPoly()
		}
	}

	for polyIndex, poly := range s {
		coeffs := make([]*big.Int, degree)
		r.PolyToBigint(poly, 1, coeffs)

		for coeffIndex, secret := range coeffs {
			randomShares := make([]*big.Int, k-1)
			for i := 0; i < k-1; i++ {
				randomShares[i] = GetRandomInt(q)
			}

			sum := big.NewInt(0)
			for i := 0; i < k-1; i++ {
				lambdaCoeff := new(big.Int).SetUint64(lambdas[i].Coeffs[0][0])
				shareTerm := new(big.Int).Mul(randomShares[i], lambdaCoeff)
				sum.Add(sum, shareTerm)
				sum.Mod(sum, q)
			}

			lastShare := new(big.Int).Sub(secret, sum)
			lastShareLambdaCoeff := new(big.Int).SetUint64(lambdas[k-1].Coeffs[0][0])
			lastShare.Mul(lastShare, new(big.Int).ModInverse(lastShareLambdaCoeff, q))
			lastShare.Mod(lastShare, q)

			for i := 0; i < k-1; i++ {
				if shares[i][polyIndex].Coeffs[0] == nil {
					shares[i][polyIndex].Coeffs[0] = make([]uint64, degree)
				}
				shares[i][polyIndex].Coeffs[0][coeffIndex] = randomShares[i].Uint64()
			}

			if shares[k-1][polyIndex].Coeffs[0] == nil {
				shares[k-1][polyIndex].Coeffs[0] = make([]uint64, degree)
			}
			shares[k-1][polyIndex].Coeffs[0][coeffIndex] = lastShare.Uint64()
		}
	}

	return shares
}

// ComputeLagrangeCoefficients computes the Lagrange coefficients for interpolation based on the indices of available shares.
func ComputeLagrangeCoefficients(r *ring.Ring, T []int, modulus *big.Int) []ring.Poly {
	lagrangeCoefficients := make([]ring.Poly, len(T))
	for i := 0; i < len(T); i++ {
		xi := big.NewInt(int64(T[i] + 1))
		numerator := big.NewInt(1)
		denominator := big.NewInt(1)
		for j := 0; j < len(T); j++ {
			if i != j {
				xj := big.NewInt(int64(T[j] + 1))
				numerator.Mul(numerator, new(big.Int).Neg(xj))
				numerator.Mod(numerator, modulus)
				temp := new(big.Int).Sub(xi, xj)
				denominator.Mul(denominator, temp)
				denominator.Mod(denominator, modulus)
			}
		}
		denomInv := new(big.Int).ModInverse(denominator, modulus)
		coeff := new(big.Int).Mul(numerator, denomInv)
		coeff.Mod(coeff, modulus)
		lagrangePoly := r.NewPoly()
		r.SetCoefficientsBigint([]*big.Int{coeff}, lagrangePoly)
		lagrangeCoefficients[i] = lagrangePoly
	}
	return lagrangeCoefficients
}
