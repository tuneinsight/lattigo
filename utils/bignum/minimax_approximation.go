package bignum

import (
	"fmt"
	"math"
	"math/big"
	"sync"
)

// Remez implements the optimized multi-interval minimax approximation
// algorithm of Lee et al. (https://eprint.iacr.org/2020/552).
// This is an iterative algorithm that returns the minimax polynomial
// approximation of any function that is smooth over a set of interval
// [a0, b0] U [a1, b1] U ... U [ai, bi].
type Remez struct {
	RemezParameters
	Degree int

	extremePoints      []point
	localExtremePoints []point
	nbExtremePoints    int

	MaxErr, MinErr *big.Float

	Nodes  []point
	Matrix [][]*big.Float
	Vector []*big.Float
	Coeffs []*big.Float
}

type point struct {
	x, y      *big.Float
	slopesign int
}

// RemezParameters is a struct storing the parameters
// required to initialize the Remez algorithm.
type RemezParameters struct {
	// Function is the function to approximate.
	// It has to be smooth in the defined intervals.
	Function func(x *big.Float) (y *big.Float)

	// Basis is the basis to use.
	// Supported basis are: Monomial and Chebyshev
	Basis Basis

	// Intervals is the set of interval [ai, bi] on which to approximate
	// the function. Each interval also define the number of nodes (points)
	// that will be used to approximate the function inside this interval.
	// This allows the user to implement a separate algorithm that allocates
	// an optimal number of nodes per interval.
	Intervals []Interval

	// ScanStep is the size of the default step used to find the extreme points.
	// The smaller this value is, the lower the probability to miss an extreme point is
	// but the longer each iteration will be.
	// A good starting value is 2^{-10}.
	ScanStep *big.Float

	// Prec defines the bit precision of the overall computation.
	Prec uint

	// OptimalScanStep is a boolean to use a dynamic update of the scan step during each
	// iteration.
	OptimalScanStep bool
}

// NewRemez instantiates a new Remez algorithm from the provided parameters.
func NewRemez(p RemezParameters) (r *Remez) {

	r = &Remez{
		RemezParameters: p,
		MaxErr:          new(big.Float).SetPrec(p.Prec),
		MinErr:          new(big.Float).SetPrec(p.Prec),
	}

	for i := range r.Intervals {
		r.Degree += r.Intervals[i].Nodes
	}

	r.Degree -= 2

	r.Nodes = make([]point, r.Degree+2)

	r.Coeffs = make([]*big.Float, r.Degree+1)
	for i := range r.Coeffs {
		r.Coeffs[i] = new(big.Float)
	}

	r.extremePoints = make([]point, 3*r.Degree)

	for i := range r.extremePoints {
		r.extremePoints[i].x = new(big.Float)
		r.extremePoints[i].y = new(big.Float)
	}

	r.localExtremePoints = make([]point, 3*r.Degree)
	for i := range r.localExtremePoints {
		r.localExtremePoints[i].x = new(big.Float)
		r.localExtremePoints[i].y = new(big.Float)
	}

	r.Matrix = make([][]*big.Float, r.Degree+2)
	for i := range r.Matrix {
		r.Matrix[i] = make([]*big.Float, r.Degree+2)

		for j := range r.Matrix[i] {
			r.Matrix[i][j] = new(big.Float)
		}
	}

	r.Vector = make([]*big.Float, r.Degree+2)
	for i := range r.Vector {
		r.Vector[i] = new(big.Float)
	}

	return r
}

// Approximate starts the approximation process.
// maxIter: the maximum number of iterations before the approximation process is terminated.
// threshold: the minimum value that (maxErr-minErr)/minErr (the normalized absolute difference
// between the maximum and minimum approximation error over the defined intervals) must take
// before the approximation process is terminated.
func (r *Remez) Approximate(maxIter int, threshold float64) {

	decimals := int(-math.Log(threshold)/math.Log(10)+0.5) + 10

	r.initialize()

	for i := 0; i < maxIter; i++ {

		// Solves the linear system and gets the new set of coefficients
		r.getCoefficients()

		// Finds the extreme points of p(x) - f(x) (where the absolute error is max)
		r.findExtremePoints()

		// Choose the new nodes based on the set of extreme points
		r.chooseNewNodes()

		nErr := new(big.Float).Sub(r.MaxErr, r.MinErr)
		nErr.Quo(nErr, r.MinErr)

		fmt.Printf("Iteration: %2d - %.*f\n", i, decimals, nErr)

		if nErr.Cmp(new(big.Float).SetFloat64(threshold)) < 1 {
			break
		}

	}
}

// ShowCoeffs prints the coefficient of the approximate
// prec: the bit precision of the printed values.
func (r *Remez) ShowCoeffs(prec int) {
	fmt.Printf("{")
	for _, c := range r.Coeffs {
		fmt.Printf("%.*f, ", prec, c)
	}
	fmt.Println("}")
}

// ShowError prints the minimum and maximum error of the approximate
// prec: the bit precision of the printed values.
func (r *Remez) ShowError(prec int) {
	fmt.Printf("MaxErr: %.*f\n", prec, r.MaxErr)
	fmt.Printf("MinErr: %.*f\n", prec, r.MinErr)
}

func (r *Remez) initialize() {

	var idx int

	switch r.Basis {
	case Monomial:

		for _, inter := range r.Intervals {

			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			A := &inter.A
			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			B := &inter.B

			nodes := inter.Nodes

			for j := 0; j < nodes; j++ {

				x := new(big.Float).Sub(B, A)
				x.Mul(x, NewFloat(float64(j+1)/float64(nodes+1), r.Prec))
				x.Add(x, A)

				y := r.Function(x)

				r.Nodes[idx+j].x = x
				r.Nodes[idx+j].y = y
			}

			idx += nodes
		}

	case Chebyshev:

		for _, inter := range r.Intervals {

			nodes := chebyshevNodes(inter.Nodes, inter)

			for j := range nodes {
				r.Nodes[idx+j].x = nodes[j]
				r.Nodes[idx+j].y = r.Function(nodes[j])
			}

			idx += len(nodes)
		}
	}
}

func (r *Remez) getCoefficients() {

	// Constructs the linear system
	// | 1 x0 x0^2 x0^3 ...  1 | f(x0)
	// | 1 x1 x1^2 x1^3 ... -1 | f(x1)
	// | 1 x2 x2^2 x2^3 ...  1 | f(x2)
	// | 1 x3 x3^2 x3^3 ... -1 | f(x3)
	// |          .            |   .
	// |          .            |   .
	// |          .            |   .

	switch r.Basis {
	case Monomial:
		for i := 0; i < r.Degree+2; i++ {
			r.Matrix[i][0] = NewFloat(1, r.Prec)
			for j := 1; j < r.Degree+1; j++ {
				r.Matrix[i][j].Mul(r.Nodes[i].x, r.Matrix[i][j-1])
			}
		}
	case Chebyshev:
		for i := 0; i < r.Degree+2; i++ {
			chebyshevBasisInPlace(r.Degree+1, r.Nodes[i].x, Interval{A: r.Intervals[0].A, B: r.Intervals[len(r.Intervals)-1].B}, r.Matrix[i])
		}
	}

	for i := 0; i < r.Degree+2; i++ {
		if i&1 == 0 {
			r.Matrix[i][r.Degree+1] = NewFloat(-1, r.Prec)
		} else {
			r.Matrix[i][r.Degree+1] = NewFloat(1, r.Prec)
		}
	}

	/*
		for i := 0; i < r.Degree+2; i++{
			for j := 0; j < r.Degree+2; j++{
				fmt.Printf("%v\n", r.Matrix[i][j])
			}
			fmt.Println()
		}
		fmt.Println()
	*/

	for i := 0; i < r.Degree+2; i++ {
		r.Vector[i].Set(r.Nodes[i].y)
	}

	// Solves the linear system
	solveLinearSystemInPlace(r.Matrix, r.Vector)

	// Updates the new [x0, x1, ..., xi]
	for i := 0; i < r.Degree+1; i++ {
		r.Coeffs[i].Set(r.Vector[i])
	}
}

func (r *Remez) findExtremePoints() {

	r.nbExtremePoints = 0

	// e = p(x) - f(x) over [a, b]
	fErr := func(x *big.Float) (y *big.Float) {
		return new(big.Float).Sub(r.eval(x), r.Function(x))
	}

	for j := 0; j < len(r.Intervals); j++ {

		points := r.findLocalExtrempointsWithSlope(fErr, r.Intervals[j])

		for i, j := r.nbExtremePoints, 0; i < r.nbExtremePoints+len(points); i, j = i+1, j+1 {
			r.extremePoints[i].x.Set(points[j].x)
			r.extremePoints[i].y.Set(points[j].y)
			r.extremePoints[i].slopesign = points[j].slopesign
		}

		r.nbExtremePoints += len(points)
	}

	// show error message
	if r.nbExtremePoints < r.Degree+2 {
		panic("number of extrem points is smaller than deg + 2, some points have been missed, consider reducing the size of the initial scan step or the approximation degree")
	}
}

// ChooseNewNodes implements Algorithm 3 of High-Precision Bootstrapping
// of RNS-CKKS Homomorphic Encryption Using Optimal Minimax Polynomial
// Approximation and Inverse Sine Function (https://eprint.iacr.org/2020/552).
// This is an optimized Go reimplementation of Remez::choosemaxs at
// https://github.com/snu-ccl/FHE-MP-CNN/blob/main-3.6.6/cnn_ckks/common/Remez.cpp
func (r *Remez) chooseNewNodes() {

	// Allocates the list of new nodes
	newNodes := []point{}

	// Retrieve the list of extrem points
	extremePoints := r.extremePoints

	// Resets max and min error
	r.MaxErr.SetFloat64(0)
	r.MinErr.SetFloat64(1e15)

	//=========================
	//========= PART 1 ========
	//=========================

	// Line 1 to 8 of Algorithm 3

	// The first part of the algorithm is to remove
	// consecutive extreme points with the same slope sign,
	// which will ensure that new linear system has a
	// solution by the Haar condition.

	// Stores consecutive extreme points with the same slope sign
	// It is unlikely that more that two consecutive extreme points
	// will have the same slope sign.
	idxAdjSameSlope := []int{}

	// To find the maximum value between extreme points that have the
	// same slope sign.
	maxpoint := new(big.Float)

	// Tracks the total number of extreme points iterated on
	ind := 0
	for ind < r.nbExtremePoints {

		// If idxAdjSameSlope is empty then adds the next point
		if len(idxAdjSameSlope) == 0 {
			idxAdjSameSlope = append(idxAdjSameSlope, ind)
			ind++
		} else {

			// If the slope of two consecutive extreme points is not alternating in sign
			// then adds the point index to the temporary array
			if extremePoints[ind-1].slopesign*extremePoints[ind].slopesign == 1 {
				mid := new(big.Float).Add(extremePoints[ind-1].x, extremePoints[ind].x)
				mid.Quo(mid, new(big.Float).SetInt64(2))
				idxAdjSameSlope = append(idxAdjSameSlope, ind)
				ind++
			} else {

				maxpoint.SetFloat64(0)

				// If the next point has alternating sign, then iterates over all the index in the temporary array
				// with extreme points whose slope is of the same sign and looks for the one with the maximum
				// absolute value
				maxIdx := 0
				for i := range idxAdjSameSlope {
					if maxpoint.Cmp(new(big.Float).Abs(extremePoints[idxAdjSameSlope[i]].y)) == -1 {
						maxpoint.Abs(extremePoints[idxAdjSameSlope[i]].y)
						maxIdx = idxAdjSameSlope[i]
					}
				}

				// Adds to the new nodes the extreme points whose absolute value is the largest
				// between all consecutive extreme points with the same slope sign
				newNodes = append(newNodes, extremePoints[maxIdx])
				idxAdjSameSlope = []int{}
			}
		}
	}

	// The above loop might terminate without flushing the array of extreme points
	// with the same slope sign, the second part of the loop is called one last time.
	maxpoint.SetInt64(0)
	maxIdx := 0
	for i := range idxAdjSameSlope {
		if maxpoint.Cmp(new(big.Float).Abs(extremePoints[idxAdjSameSlope[i]].y)) == -1 {
			maxpoint.Abs(extremePoints[idxAdjSameSlope[i]].y)
			maxIdx = idxAdjSameSlope[i]
		}
	}

	newNodes = append(newNodes, extremePoints[maxIdx])

	if len(newNodes) < r.Degree+2 {
		panic("number of alternating extreme points is less than deg+2, some points have been missed, consider reducing the size of the initial scan step or the approximation degree")
	}

	//=========================
	//========= PART 2 ========
	//=========================

	// Lines 11 to 24 of Algorithm 3

	// Choosing the new nodes if the set of alternating extreme points
	// is larger than degree+2.

	minPair := new(big.Float)
	tmp := new(big.Float)

	// Loops run as long as the number of extreme points is not equal to deg+2 (the dimension of the linear system)
	var minIdx int
	for len(newNodes) > r.Degree+2 {

		minPair.SetFloat64(1e15)

		// If the number of remaining extreme points is one more than the number needed
		// then we can remove only one point
		if len(newNodes) == r.Degree+3 {

			// Removes the largest one between the first and the last
			if new(big.Float).Abs(newNodes[0].y).Cmp(new(big.Float).Abs(newNodes[len(newNodes)-1].y)) == 1 {
				newNodes = newNodes[:len(newNodes)-1]
			} else {
				newNodes = newNodes[1:]
			}

			// If the number of remaining extreme points is two more than the number needed
			// then we can remove two points.
		} else if len(newNodes) == r.Degree+4 {

			// Finds the minimum index of the sum of two adjacent points
			for i := range newNodes {
				tmp.Add(new(big.Float).Abs(newNodes[i].y), new(big.Float).Abs(newNodes[(i+1)%len(newNodes)].y))
				if minPair.Cmp(tmp) == 1 {
					minPair.Set(tmp)
					minIdx = i
				}
			}

			// If the index is the last, then remove the first and last points
			if minIdx == len(newNodes)-1 {
				newNodes = newNodes[1:]
				// Else remove the two consecutive points
			} else {
				newNodes = append(newNodes[:minIdx], newNodes[minIdx+2:]...)
			}

			// If the number of remaining extreme points is more four over the number needed
			// then remove up to two points, prioritizing the first and last points.
		} else {

			// Finds the minimum index of the sum of two adjacent points
			for i := range newNodes[:len(newNodes)-1] {

				tmp.Add(new(big.Float).Abs(newNodes[i].y), new(big.Float).Abs(newNodes[i+1].y))

				if minPair.Cmp(tmp) == 1 {
					minPair.Set(tmp)
					minIdx = i
				}
			}

			// If the first element is included in the smallest sum, then removes it
			if minIdx == 0 {
				newNodes = newNodes[1:]
				// If the last element is included in the smallest sum, then removes it
			} else if minIdx == len(newNodes)-2 {
				newNodes = newNodes[:len(newNodes)-1]
				// Else removes the two consecutive points adding to the smallest sum
			} else {
				newNodes = append(newNodes[:minIdx], newNodes[minIdx+2:]...)
			}
		}
	}

	// Assigns the new points to the nodes and computes the min and max error
	for i := 0; i < r.Degree+2; i++ {

		// Deep copy
		r.Nodes[i].x.Set(newNodes[i].x)
		r.Nodes[i].y = r.Function(r.Nodes[i].x)      // we must evaluate, because Y was the error Function)
		r.Nodes[i].slopesign = newNodes[i].slopesign // should have alternating sign

		if r.MaxErr.Cmp(new(big.Float).Abs(newNodes[i].y)) == -1 {
			r.MaxErr.Abs(newNodes[i].y)
		}

		if r.MinErr.Cmp(new(big.Float).Abs(newNodes[i].y)) == 1 {
			r.MinErr.Abs(newNodes[i].y)
		}
	}
}

// findLocalExtrempointsWithSlope finds local extrema/minima of a function.
// It starts by scanning the interval with a pre-defined window size, until it finds that the function is concave or convex
// in this window. Then it uses a binary search to find the local maximum/minimum in this window. The process is repeated
// until the entire interval has been scanned.
// This is an optimized Go re-implementation of the method find_extreme that can be found at
// https://github.com/snu-ccl/FHE-MP-CNN/blob/main-3.6.6/cnn_ckks/common/MinicompFunc.cpp
func (r *Remez) findLocalExtrempointsWithSlope(fErr func(*big.Float) (y *big.Float), interval Interval) []point {

	extrempoints := r.localExtremePoints
	prec := r.Prec
	scan := r.ScanStep

	var slopeLeft, slopeRight, s int

	scanMid := new(big.Float)
	scanRight := new(big.Float)
	scanLeft := new(big.Float)
	fErrLeft := new(big.Float)
	fErrRight := new(big.Float)

	nbextrempoints := 0
	extrempoints[nbextrempoints].x.Set(&interval.A)
	extrempoints[nbextrempoints].y.Set(fErr(&interval.A))
	extrempoints[nbextrempoints].slopesign = extrempoints[nbextrempoints].y.Cmp(new(big.Float))
	nbextrempoints++

	optScan := new(big.Float).Set(scan)

	if r.OptimalScanStep {
		s = 15
		optScan.Quo(scan, NewFloat(1e15, prec))
	} else {
		optScan.Set(scan)
	}

	scanMid.Set(&interval.A)
	scanRight.Add(&interval.A, optScan)
	fErrLeft.Set(fErr(scanMid))
	fErrRight.Set(fErr(scanRight))

	if slopeRight = fErrRight.Cmp(fErrLeft); slopeRight == 0 {
		panic("slope 0 occurred: consider increasing the precision")
	}

	for {

		if r.OptimalScanStep {

			for i := 0; i < s; i++ {

				pow10 := NewFloat(math.Pow(10, float64(i)), prec)

				// start + 10*scan/pow(10,i)
				a := new(big.Float).Mul(scan, NewFloat(10, prec))
				a.Quo(a, pow10)
				a.Add(&interval.A, a)

				// end - 10*scan/pow(10,i)
				b := new(big.Float).Mul(scan, NewFloat(10, prec))
				b.Quo(b, pow10)
				b.Sub(&interval.B, b)

				// a < scanRight && scanRight < b
				if a.Cmp(scanRight) == -1 && scanRight.Cmp(b) == -1 {
					optScan.Quo(scan, pow10)
					break
				}

				if i == s-1 {
					optScan.Quo(scan, pow10)
					optScan.Quo(optScan, NewFloat(10, prec))
				}
			}

		} else {
			optScan.Set(scan)
		}

		// Breaks when the scan window gets out of the interval
		if new(big.Float).Add(scanRight, optScan).Cmp(&interval.B) >= 0 {
			break
		}

		slopeLeft = slopeRight
		scanLeft.Set(scanMid)
		scanMid.Set(scanRight)
		scanRight.Add(scanMid, optScan)

		fErrLeft.Set(fErrRight)
		fErrRight.Set(fErr(scanRight))

		if slopeRight = fErrRight.Cmp(fErrLeft); slopeRight == 0 {
			panic("slope 0 occurred: consider increasing the precision")
		}

		// Positive and negative slope (concave)
		if slopeLeft == 1 && slopeRight == -1 {
			findLocalMaximum(fErr, scanLeft, scanRight, prec, &extrempoints[nbextrempoints])
			nbextrempoints++
			// Negative and positive slope (convex)
		} else if slopeLeft == -1 && slopeRight == 1 {
			findLocalMinimum(fErr, scanLeft, scanRight, prec, &extrempoints[nbextrempoints])
			nbextrempoints++
		}
	}

	extrempoints[nbextrempoints].x.Set(&interval.B)
	extrempoints[nbextrempoints].y.Set(fErr(&interval.B))
	extrempoints[nbextrempoints].slopesign = extrempoints[nbextrempoints].y.Cmp(new(big.Float))
	nbextrempoints++

	return extrempoints[:nbextrempoints]
}

// findLocalMaximum finds the local maximum of a function that is concave in a given window.
func findLocalMaximum(fErr func(x *big.Float) (y *big.Float), start, end *big.Float, prec uint, p *point) {

	windowStart := new(big.Float).Set(start)
	windowEnd := new(big.Float).Set(end)
	quarter := new(big.Float).Sub(windowEnd, windowStart)
	quarter.Quo(quarter, NewFloat(4, prec))

	for i := 0; i < int(prec); i++ {

		// Obtains the sign of the err Function in the interval (normalized and zeroed)
		// 0: [0.00, 0.25]
		// 1: [0.25, 0.50]
		// 2: [0.50, 0.75]
		// 3: [0.75, 1.00]

		slopeWin0, slopeWin1, slopeWin2, slopeWin3 := slopes(fErr, windowStart, windowEnd, quarter)

		// Look for a sign change between the 4 intervals.
		// Since we are here in a concave Function, we look
		// for the point in the interval where the sign of the
		// err Function changes.

		// Sign change occurs between [0, 0.5]
		if slopeWin0 == 1 && slopeWin1 == -1 {

			// Reduces the windowEnd from 1 to 0.5
			windowEnd.Sub(windowEnd, quarter)
			windowEnd.Sub(windowEnd, quarter)

			// Divides the scan step by half
			quarter.Quo(quarter, NewFloat(2.0, prec))

			// Sign change occurs between [0.25, 0.75]
		} else if slopeWin1 == 1 && slopeWin2 == -1 {

			// Increases windowStart from 0 to 0.25
			windowStart.Add(windowStart, quarter)

			// Decreases windowEnd from 1 to 0.75
			windowEnd.Sub(windowEnd, quarter)

			// Divides the scan step by half
			quarter.Quo(quarter, NewFloat(2.0, prec))

			// Sign change occurs between [0.5, 1.0]
		} else if slopeWin2 == 1 && slopeWin3 == -1 {

			// Increases windowStart fro 0 to 0.5
			windowStart.Add(windowStart, quarter)
			windowStart.Add(windowStart, quarter)

			// Divides the scan step by half
			quarter.Quo(quarter, NewFloat(2.0, prec))
		}
	}

	p.x.Quo(new(big.Float).Add(windowStart, windowEnd), NewFloat(2, prec))
	p.y.Set(fErr(p.x))
	p.slopesign = 1
}

// findLocalMaximum finds the local maximum of a function that is convex in a given window.
func findLocalMinimum(fErr func(x *big.Float) (y *big.Float), start, end *big.Float, prec uint, p *point) {

	windowStart := new(big.Float).Set(start)
	windowEnd := new(big.Float).Set(end)
	quarter := new(big.Float).Sub(windowEnd, windowStart)
	quarter.Quo(quarter, NewFloat(4, prec))

	for i := 0; i < int(prec); i++ {

		// Obtains the sign of the err Function in the interval (normalized and zeroed)
		// 0: [0.00, 0.25]
		// 1: [0.25, 0.50]
		// 2: [0.50, 0.75]
		// 3: [0.75, 1.00]
		slopeWin0, slopeWin1, slopeWin2, slopeWin3 := slopes(fErr, windowStart, windowEnd, quarter)

		// Look for a sign change between the 4 intervals.
		// Since we are here in a convex Function, we look
		// for the point in the interval where the sign of the
		// err Function changes.

		// Sign change occurs between [0, 0.5]
		if slopeWin0 == -1 && slopeWin1 == 1 {

			// Reduces the windowEnd from 1 to 0.5
			windowEnd.Sub(windowEnd, quarter)
			windowEnd.Sub(windowEnd, quarter)

			// Divides the scan step by half
			quarter.Quo(quarter, NewFloat(2.0, prec))

			// Sign change occurs between [0.25, 0.75]
		} else if slopeWin1 == -1 && slopeWin2 == 1 {

			// Increases windowStart from 0 to 0.25
			windowStart.Add(windowStart, quarter)

			// Decreases windowEnd from 1 to 0.75
			windowEnd.Sub(windowEnd, quarter)

			// Divides the scan step by half
			quarter.Quo(quarter, NewFloat(2.0, prec))

			// Sign change occurs between [0.5, 1.0]
		} else if slopeWin2 == -1 && slopeWin3 == 1 {

			// Increases windowStart fro 0 to 0.5
			windowStart.Add(windowStart, quarter)
			windowStart.Add(windowStart, quarter)

			// Divides the scan step by half
			quarter.Quo(quarter, NewFloat(2.0, prec))
		}

	}

	p.x.Quo(new(big.Float).Add(windowStart, windowEnd), NewFloat(2, prec))
	p.y.Set(fErr(p.x))
	p.slopesign = -1
}

// slopes takes a window, divides it into four intervals and computes the sign of the slope of the error function in each sub-interval.
func slopes(fErr func(x *big.Float) (y *big.Float), searchStart, searchEnd, searchquarter *big.Float) (searchslopeLeft, searchslopeRight, searchInc3, searchInc4 int) {

	slope := func(fErr func(x *big.Float) (y *big.Float), start, end *big.Float) (sign int) {

		if a := fErr(start); a.Cmp(fErr(end)) == -1 {
			return 1
		}

		return -1
	}

	var wg sync.WaitGroup
	wg.Add(4)
	go func() {

		// [start, start + sc]
		start := searchStart
		end := new(big.Float).Add(searchStart, searchquarter)

		searchslopeLeft = slope(fErr, start, end)
		wg.Done()
	}()

	go func() {

		//[start + sc, start + 2*sc]
		start := new(big.Float).Add(searchStart, searchquarter)
		end := new(big.Float).Add(searchStart, searchquarter)
		end.Add(end, searchquarter)

		searchslopeRight = slope(fErr, start, end)
		wg.Done()
	}()

	go func() {

		// [start + 2*sc, enc-sc]
		start := new(big.Float).Add(searchStart, searchquarter)
		start.Add(start, searchquarter)
		end := new(big.Float).Sub(searchEnd, searchquarter)

		searchInc3 = slope(fErr, start, end)
		wg.Done()
	}()

	go func() {

		// [end-sc, end]
		start := new(big.Float).Sub(searchEnd, searchquarter)
		end := searchEnd

		searchInc4 = slope(fErr, start, end)
		wg.Done()
	}()

	wg.Wait()

	return
}

func (r *Remez) eval(x *big.Float) (y *big.Float) {
	switch r.Basis {
	case Monomial:
		return MonomialEval(x, r.Coeffs)
	case Chebyshev:
		return ChebyshevEval(x, r.Coeffs, Interval{A: r.Intervals[0].A, B: r.Intervals[len(r.Intervals)-1].B, Nodes: r.Degree + 1})
	default:
		panic("invalid Basis")
	}
}

// solves for y the system matrix * y = vector using Gaussian elimination.
func solveLinearSystemInPlace(matrix [][]*big.Float, vector []*big.Float) {

	n, m := len(matrix), len(matrix[0])

	var tmp = new(big.Float)
	for i := 0; i < n; i++ {

		a := matrix[i][i]

		vector[i].Quo(vector[i], a)

		for j := m - 1; j >= i; j-- {
			b := matrix[i][j]
			b.Quo(b, a)
		}

		for j := i + 1; j < m; j++ {
			c := matrix[j][i]
			vector[j].Sub(vector[j], tmp.Mul(vector[i], c))
			for k := m - 1; k >= i; k-- {
				matrix[j][k].Sub(matrix[j][k], tmp.Mul(matrix[i][k], c))
			}
		}
	}

	for i := m - 1; i > 0; i-- {
		c := vector[i]
		for j := i - 1; j >= 0; j-- {
			vector[j].Sub(vector[j], tmp.Mul(matrix[j][i], c))
		}
	}
}
