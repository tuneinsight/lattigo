package sign

import (
	"fmt"
	"log"
	"math/big"
	"strings"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
	"github.com/zeebo/blake3"
)

// MatrixVectorMul performs matrix-vector multiplication.
func MatrixVectorMulNTT(r *ring.Ring, M structs.Matrix[ring.Poly], vec structs.Vector[ring.Poly], result structs.Vector[ring.Poly]) {
	// Convert all elements of the matrix and the vector to the NTT domain
	ConvertMatrixToNTT(r, M)
	ConvertVectorToNTT(r, vec)

	// Perform the multiplications coefficient-wise
	for i := range M {
		result[i] = r.NewPoly()
		for j := range (M)[i] {
			r.MulCoeffsMontgomeryThenAdd(M[i][j], vec[j], result[i])
		}
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertVectorFromNTT(r, result)
	ConvertMatrixFromNTT(r, M)
	ConvertVectorFromNTT(r, vec)
}

// MatrixMatrixMul performs matrix-matrix multiplication.
func MatrixMatrixMulNTT(r *ring.Ring, M1, M2 structs.Matrix[ring.Poly], result structs.Matrix[ring.Poly]) {
	if M1 == nil || M2 == nil || len(M1) == 0 || len(M2) == 0 || len((M1)[0]) != len(M2) {
		log.Fatalf("Matrix dimensions are not compatible for multiplication.")
		return
	}

	m := len(M1)
	p := len(M1[0]) // Assuming all rows in M1 are of the same length
	n := len(M2[0]) // Assuming all rows in M2 are of the same length

	// Convert all elements of M1 and M2 to the NTT domain
	ConvertMatrixToNTT(r, M1)
	ConvertMatrixToNTT(r, M2)

	// Initialize the result matrix with zeros
	for i := 0; i < m; i++ {
		result[i] = make([]ring.Poly, n)
		for j := 0; j < n; j++ {
			result[i][j] = r.NewPoly()
		}
	}

	temp := r.NewPoly()

	// Perform matrix multiplication coefficient-wise
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < p; k++ {
				// r.MForm(M1[i][k], temp)
				r.MulCoeffsMontgomeryThenAdd(temp, M2[k][j], result[i][j])
			}
		}
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertMatrixFromNTT(r, result)
	ConvertMatrixFromNTT(r, M1)
	ConvertMatrixFromNTT(r, M2)
}

// VectorPolyMul performs element-wise multiplication of a vector by a polynomial.
func VectorPolyMulNTT(r *ring.Ring, vec structs.Vector[ring.Poly], poly ring.Poly, result structs.Vector[ring.Poly]) {
	// Convert the polynomial to the NTT domain
	r.NTT(poly, poly)

	// Convert all elements of the vector to the NTT domain
	ConvertVectorToNTT(r, vec)
	ConvertVectorToNTT(r, result)
	temp := r.NewPoly()

	// Perform the multiplications coefficient-wise
	for i := range vec {
		r.MulCoeffsMontgomery(temp, poly, result[i])
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertVectorFromNTT(r, result)
	ConvertVectorFromNTT(r, vec)
	r.INTT(poly, poly)
}

// No included NTT

// MatrixVectorMul performs matrix-vector multiplication.
func MatrixVectorMul(r *ring.Ring, M structs.Matrix[ring.Poly], vec structs.Vector[ring.Poly], result structs.Vector[ring.Poly]) {
	for i := range M {
		for j := range M[i] {
			r.MulCoeffsMontgomeryThenAdd(M[i][j], vec[j], result[i])
		}
	}
}

// MatrixMatrixMul performs matrix-matrix multiplication.
func MatrixMatrixMul(r *ring.Ring, M1, M2 structs.Matrix[ring.Poly], result structs.Matrix[ring.Poly]) {
	if M1 == nil || M2 == nil || len(M1) == 0 || len(M2) == 0 || len((M1)[0]) != len(M2) {
		log.Fatalf("Matrix dimensions are not compatible for multiplication.")
		return
	}

	m := len(M1)
	p := len(M1[0])
	n := len(M2[0])

	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < p; k++ {
				r.MulCoeffsMontgomeryThenAdd(M1[i][k], M2[k][j], result[i][j])
			}
		}
	}
}

// VectorPolyMul performs element-wise multiplication of a vector by a polynomial.
func VectorPolyMul(r *ring.Ring, vec structs.Vector[ring.Poly], poly ring.Poly, result structs.Vector[ring.Poly]) {
	for i := range vec {
		r.MulCoeffsMontgomery(vec[i], poly, result[i])
	}
}

// MatrixAdd adds two matrices of ring.Poly element-wise and stores the result in a given result matrix.
func MatrixAdd(r *ring.Ring, M1, M2, result structs.Matrix[ring.Poly]) {
	if M1 == nil || M2 == nil || len(M1) == 0 || len(M2) == 0 || len(M1) != len(M2) || len((M1)[0]) != len((M2)[0]) {
		log.Fatalf("Matrix dimensions must match for element-wise addition.")
		return
	}

	m := len(M1)
	n := len(M1[0])

	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			r.Add(M1[i][j], M2[i][j], result[i][j])
		}
	}
}

// VectorAdd adds two vectors of ring.Poly element-wise and stores the result in a result vector.
func VectorAdd(r *ring.Ring, v1, v2, result structs.Vector[ring.Poly]) {
	for i := range v1 {
		r.Add(v1[i], v2[i], result[i])
	}
}

// VectorSub subtracts two vectors of ring.Poly element-wise and stores the result in a result vector.
func VectorSub(r *ring.Ring, v1, v2, result structs.Vector[ring.Poly]) {
	for i := range v1 {
		r.Sub(v1[i], v2[i], result[i])
	}
}

// SAMPLER HELPERS

// SamplePolyVector samples a vector of polynomials of a given length using the provided sampler.
func SamplePolyVector(r *ring.Ring, length int, sampler ring.Sampler, NTT bool, montgomery bool) structs.Vector[ring.Poly] {
	vector := structs.Vector[ring.Poly](make([]ring.Poly, length))
	for i := 0; i < length; i++ {
		vector[i] = sampler.ReadNew()
		if NTT {
			r.NTT(vector[i], vector[i])
		}
		if montgomery {
			r.MForm(vector[i], vector[i])
		}
	}
	return vector
}

// SamplePolyMatrix samples a matrix of polynomials with given dimensions (rows and cols) using the provided sampler.
func SamplePolyMatrix(r *ring.Ring, rows, cols int, sampler ring.Sampler, NTT bool, montgomery bool) structs.Matrix[ring.Poly] {
	matrix := structs.Matrix[ring.Poly](make([][]ring.Poly, rows))
	for i := 0; i < rows; i++ {
		matrix[i] = make([]ring.Poly, cols)
		for j := 0; j < cols; j++ {
			matrix[i][j] = sampler.ReadNew()
			if NTT {
				r.NTT(matrix[i][j], matrix[i][j])
			}
			if montgomery {
				r.MForm(matrix[i][j], matrix[i][j])
			}
		}
	}
	return matrix
}

// PRINT FUNCTIONS

func PrintMatrix(label string, matrix structs.Matrix[ring.Poly]) {
	log.Println(label)
	for i, row := range matrix {
		for j, poly := range row {
			log.Printf("[%d][%d]: %s\n", i, j, formatCoeffs(poly.Coeffs[0]))
		}
	}
}

func PrintVector(label string, vector structs.Vector[ring.Poly]) {
	log.Println(label)
	for i, poly := range vector {
		log.Printf("[%d]: %s\n", i, formatCoeffs(poly.Coeffs[0]))
	}
}

func PrintPolynomial(label string, poly ring.Poly) {
	log.Println(label)
	log.Printf("%s\n", formatCoeffs(poly.Coeffs[0]))
}

func formatCoeffs(coeffs []uint64) string {
	var coeffStr []string
	for _, coeff := range coeffs {
		coeffStr = append(coeffStr, fmt.Sprintf("%v", coeff))
	}
	return strings.Join(coeffStr, ", ")
}

// NTT CONVERSION

// ConvertMatrixToNTT converts a matrix of polynomials to the NTT domain.
func ConvertMatrixToNTT(r *ring.Ring, M structs.Matrix[ring.Poly]) {
	for i := range M {
		for j := range M[i] {
			r.NTT(M[i][j], M[i][j])
			r.MForm(M[i][j], M[i][j])
		}
	}
}

// ConvertMatrixFromNTT converts a matrix of polynomials from the NTT domain back to the standard domain.
func ConvertMatrixFromNTT(r *ring.Ring, M structs.Matrix[ring.Poly]) {
	for i := range M {
		for j := range M[i] {
			r.IMForm(M[i][j], M[i][j])
			r.INTT(M[i][j], M[i][j])
		}
	}
}

// ConvertVectorToNTT converts a vector of polynomials to the NTT domain.
func ConvertVectorToNTT(r *ring.Ring, vec structs.Vector[ring.Poly]) {
	for i := range vec {
		r.NTT(vec[i], vec[i])
		r.MForm(vec[i], vec[i])
	}
}

// ConvertVectorFromNTT converts a vector of polynomials from the NTT domain back to the standard domain.
func ConvertVectorFromNTT(r *ring.Ring, vec structs.Vector[ring.Poly]) {
	for i := range vec {
		r.IMForm(vec[i], vec[i])
		r.INTT(vec[i], vec[i])
	}
}

// INITIALIZE HELPERS

// InitializeVector creates and returns a vector of the given length, initializing each element as a new polynomial.
func InitializeVector(r *ring.Ring, length int) structs.Vector[ring.Poly] {
	vector := make(structs.Vector[ring.Poly], length)
	for i := range vector {
		vector[i] = r.NewPoly()
	}
	return vector
}

// InitializeMatrix creates and returns a matrix of the given dimensions, initializing each element as a new polynomial.
func InitializeMatrix(r *ring.Ring, rows, cols int) structs.Matrix[ring.Poly] {
	matrix := make(structs.Matrix[ring.Poly], rows)
	for i := range matrix {
		matrix[i] = make([]ring.Poly, cols)
		for j := range matrix[i] {
			matrix[i][j] = r.NewPoly()
		}
	}
	return matrix
}

// PrintBigIntVector prints a vector of slices of *big.Int elements
func PrintBigIntVector(label string, vector structs.Vector[[]*big.Int]) {
	log.Println(label)
	for i, slice := range vector {
		log.Printf("[%d]: %s\n", i, FormatBigIntSlice(slice))
	}
}

func FormatBigIntSlice(slice []*big.Int) string {
	var sliceStr []string
	for _, val := range slice {
		sliceStr = append(sliceStr, val.String())
	}
	return strings.Join(sliceStr, ", ")
}

// Copy map helpers

func CopyMatrixMap(original map[int]structs.Matrix[ring.Poly]) map[int]structs.Matrix[ring.Poly] {
	copy := make(map[int]structs.Matrix[ring.Poly])
	for key, value := range original {
		copy[key] = value
	}
	return copy
}

func CopyVectorMap(original map[int]structs.Vector[ring.Poly]) map[int]structs.Vector[ring.Poly] {
	copy := make(map[int]structs.Vector[ring.Poly])
	for key, value := range original {
		copy[key] = value
	}
	return copy
}

// CalculateBetaDelta computes ((B * p) - q) / (2 * q) as a big.Int
func CalculateBetaDelta(p uint64, B float64, q uint64) *big.Int {
	pInt := new(big.Int).SetUint64(p)
	BInt := new(big.Int)
	BInt.SetString(fmt.Sprintf("%.0f", B), 10)
	qInt := new(big.Int).SetUint64(q)

	Bp := new(big.Int).Mul(BInt, pInt)
	BpSubQ := new(big.Int).Sub(Bp, qInt)
	twoQ := new(big.Int).Mul(big.NewInt(2), qInt)
	betaDelta := new(big.Int).Div(BpSubQ, twoQ)
	log.Println("betaddelta", betaDelta)
	return betaDelta
}

func PrintSignRepresentation(r *ring.Ring, poly ring.Poly, modulus uint64) {
	qBig := new(big.Int).SetUint64(modulus)
	halfQ := new(big.Int).Div(qBig, big.NewInt(2))

	coeffs := make([]*big.Int, r.N())
	r.PolyToBigint(poly, 1, coeffs)

	for i := 0; i < len(coeffs); i++ {
		if coeffs[i].Cmp(halfQ) > 0 {
			coeffs[i].Sub(coeffs[i], qBig)
		}
	}

	log.Println("Coeffs", coeffs)
}

func PrintSignRepresentationVector(r *ring.Ring, vec structs.Vector[ring.Poly], modulus uint64) {
	for i := 0; i < len(vec); i++ {
		PrintSignRepresentation(r, vec[i], modulus)
	}
}

func PrintSignRepresentationMatrix(r *ring.Ring, matrix structs.Matrix[ring.Poly], modulus uint64) {
	for i := 0; i < len(matrix); i++ {
		PrintSignRepresentationVector(r, matrix[i], modulus)
	}
}

// Rounding

// RoundCoefficients rounds each coefficient of the polynomial as specified
func RoundCoefficients(r *ring.Ring, newRing *ring.Ring, poly ring.Poly, roundingVal uint) ring.Poly {
	roundedPoly := newRing.NewPoly()
	for i := range poly.Coeffs[0] {
		coeff := poly.Coeffs[0][i]
		roundedCoeff := (coeff + (1 << (roundingVal - 1))) >> roundingVal
		roundedPoly.Coeffs[0][i] = roundedCoeff
	}
	return roundedPoly
}

// RoundVector rounds each polynomial in the vector
func RoundVector(r *ring.Ring, newRing *ring.Ring, v structs.Vector[ring.Poly], roundingVal uint) structs.Vector[ring.Poly] {
	roundedVector := InitializeVector(newRing, len(v))
	for i := range v {
		roundedVector[i] = RoundCoefficients(r, newRing, v[i], roundingVal)
	}
	return roundedVector
}

// RestoreCoefficients multiplies each coefficient of the polynomial by 2^roundingVal
func RestoreCoefficients(r *ring.Ring, newRing *ring.Ring, poly ring.Poly, roundingVal uint) ring.Poly {
	restoredPoly := r.NewPoly()
	for i := range poly.Coeffs[0] {
		coeff := poly.Coeffs[0][i]
		restoredCoeff := coeff << roundingVal
		restoredPoly.Coeffs[0][i] = restoredCoeff
	}
	return restoredPoly
}

// RestoreVector multiplies each polynomial in the vector by 2^roundingVal
func RestoreVector(r *ring.Ring, newRing *ring.Ring, v structs.Vector[ring.Poly], roundingVal uint) structs.Vector[ring.Poly] {
	restoredVector := InitializeVector(r, len(v))
	for i := range v {
		restoredVector[i] = RestoreCoefficients(r, newRing, v[i], roundingVal)
	}
	return restoredVector
}

// Global variable to hold precomputed randomness
var PrecomputedRandomness []byte
var RandomnessIndex int

// PrecomputeRandomness precomputes all necessary randomness and stores it in the global variable
func PrecomputeRandomness(size int, key []byte) {
	RandomnessIndex = 0
	hasher := blake3.New()
	hasher.Write(key)
	digest := hasher.Digest()
	PrecomputedRandomness = make([]byte, size)
	digest.Read(PrecomputedRandomness)
}

// GetRandomBytes returns the next n bytes of precomputed randomness
func GetRandomBytes(n int) []byte {
	bytes := PrecomputedRandomness[RandomnessIndex : RandomnessIndex+n]
	RandomnessIndex += n
	return bytes
}

// GetRandomInt returns a random integer from the precomputed randomness
func GetRandomInt(q *big.Int) *big.Int {
	randBytes := GetRandomBytes(len(q.Bytes()))
	randInt := new(big.Int).SetBytes(randBytes)
	return randInt.Mod(randInt, q)
}

// Convert to signed representation
func SignedRepresentation(coeffs []*big.Int, Q uint64) {
	qBig := new(big.Int).SetUint64(Q)
	halfQ := new(big.Int).Div(qBig, big.NewInt(2))

	for i := 0; i < len(coeffs); i++ {
		if coeffs[i].Cmp(halfQ) > 0 {
			coeffs[i].Sub(coeffs[i], qBig)
		}
	}
}

// GaussianEliminationModQ performs Gaussian elimination mod q and checks if the matrix is full rank
func GaussianEliminationModQ(mat [][]*big.Int, q *big.Int) bool {
	m := copyMat(mat)
	rows := len(m)
	cols := len(m[0])
	for i := 0; i < rows; i++ {
		foundPivot := false
		var pivotIdx int
		for j := 0; j < cols; j++ {
			if m[i][j].Cmp(big.NewInt(0)) != 0 {
				foundPivot = true
				pivotIdx = j
				break
			}
		}
		if !foundPivot {
			return false
		}
		invPivot := new(big.Int).ModInverse(m[i][pivotIdx], q)
		for j := 0; j < cols; j++ {
			m[i][j].Mul(m[i][j], invPivot).Mod(m[i][j], q)
		}
		for k := i + 1; k < rows; k++ {
			if m[k][pivotIdx].Cmp(big.NewInt(0)) != 0 {
				factor := new(big.Int).Set(m[k][pivotIdx])
				for j := 0; j < cols; j++ {
					m[k][j].Sub(m[k][j], new(big.Int).Mul(factor, m[i][j])).Mod(m[k][j], q)
				}
			}
		}
	}
	return true
}

func copyMat(matBig [][]*big.Int) [][]*big.Int {
	rows := len(matBig)
	cols := len(matBig[0])
	data := make([][]*big.Int, rows)
	for i := 0; i < rows; i++ {
		data[i] = make([]*big.Int, cols)
		for j := 0; j < cols; j++ {
			data[i][j] = new(big.Int).Set(matBig[i][j])
		}
	}
	return data
}
