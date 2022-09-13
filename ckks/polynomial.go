package ckks

import (
	"fmt"
	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

type Polynomial struct {
	rlwe.Polynomial
}

func NewPolynomial(polynomialBasis rlwe.PolynomialBasisType, coeffs interface{}, slotsIndex [][]int) (Polynomial, error) {
	pol, err := rlwe.NewPolynomial(polynomialBasis, coeffs, slotsIndex)
	if err != nil {
		return Polynomial{}, err
	}
	pol.Scale = &Scale{}
	return Polynomial{pol}, nil
}

func (p *Polynomial) MarshalBinary() (data []byte, err error) {
	return p.Polynomial.MarshalBinary()
}

func (p *Polynomial) UnmarshalBinary(data []byte) (err error) {
	p.Polynomial = rlwe.Polynomial{}
	p.Coefficients = rlwe.Coefficients{}
	p.Scale = &Scale{}
	return p.Polynomial.UnmarshalBinary(data)
}

func (p *Polynomial) Encode(ecd Encoder, level int, inputScale, outputScale rlwe.Scale) (Polynomial, error) {

	params := ecd.(*encoderComplex128).params

	ptPoly := rlwe.Polynomial{}

	switch p.Coefficients.Value.(type) {
	case []float64, []complex128, [][]float64, [][]complex128:

		ptPoly.Coefficients = rlwe.Coefficients{Value: [][]*rlwe.Plaintext{}}
		ptPoly.Basis = p.Basis
		ptPoly.Odd = p.Odd
		ptPoly.Even = p.Even
		ptPoly.SlotsIndex = p.SlotsIndex

		getScaledBSGSCoefficients(params, ecd, nil, level, inputScale, p.Polynomial, outputScale, &ptPoly)

	default:
		return Polynomial{}, fmt.Errorf("Polynomial.Encode(*): underlying polynomial is already encoded or encrypted")
	}

	return Polynomial{ptPoly}, nil
}

func (p *Polynomial) Encrypt(ecd Encoder, enc Encryptor, level int, inputScale, outputScale rlwe.Scale) (Polynomial, error) {

	params := ecd.(*encoderComplex128).params

	ctPoly := rlwe.Polynomial{}

	switch p.Coefficients.Value.(type) {
	case []float64, []complex128, [][]float64, [][]complex128:

		ctPoly.Coefficients = rlwe.Coefficients{Value: [][]*rlwe.Ciphertext{}}
		ctPoly.Basis = p.Basis
		ctPoly.Odd = p.Odd
		ctPoly.Even = p.Even
		ctPoly.SlotsIndex = p.SlotsIndex

		getScaledBSGSCoefficients(params, ecd, enc, level, inputScale, p.Polynomial, outputScale, &ctPoly)

	default:
		return Polynomial{}, fmt.Errorf("Polynomial.Encrypt(*): underlying polynomial is already encrypted")
	}

	return Polynomial{ctPoly}, nil
}

type dummyPolynomialEvaluator struct {
	dummyPolynomialBasis
	params Parameters
	ecd    Encoder
	enc    Encryptor
	giant  int
	baby   int
	poly   *rlwe.Polynomial
}

func getScaledBSGSCoefficients(params Parameters, ecd Encoder, enc Encryptor, level int, scale rlwe.Scale, polIn rlwe.Polynomial, targetScale rlwe.Scale, polOut *rlwe.Polynomial) {

	dummbpb := newDummyPolynomialBasis(params, &dummyCiphertext{level, scale.CopyNew()})

	giant, baby := polIn.BSGSSplit()

	isRingStandard := params.RingType() == ring.Standard

	odd, even := polIn.OddEven()

	for i := (1 << baby) - 1; i > 1; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			dummbpb.GenPower(i, isRingStandard, targetScale)
		}
	}

	for i := baby; i < giant; i++ {
		dummbpb.GenPower(1<<i, false, targetScale)
	}

	polyEval := &dummyPolynomialEvaluator{}
	polyEval.params = params
	polyEval.ecd = ecd
	polyEval.enc = enc
	polyEval.dummyPolynomialBasis = *dummbpb
	polyEval.giant = giant
	polyEval.baby = baby
	polyEval.poly = polOut

	polyEval.recurse(dummbpb.Value[1].Level-giant+1, targetScale, evalPoly{polIn, true, polIn.Degree()})
}

func (polyEval *dummyPolynomialEvaluator) recurse(targetLevel int, targetScale rlwe.Scale, polIn evalPoly) (res *dummyCiphertext) {

	params := polyEval.params

	baby := polyEval.baby

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if polIn.Degree() < (1 << baby) {

		if polIn.lead && baby > 1 && polIn.maxDeg%(1<<(baby+1)) > (1<<(baby-1)) {

			giant := int(bits.Len64(uint64(polIn.Degree())))
			baby := giant >> 1

			polyEvalBis := new(dummyPolynomialEvaluator)
			polyEvalBis.params = polyEval.params
			polyEvalBis.ecd = polyEval.ecd
			polyEvalBis.enc = polyEval.enc
			polyEvalBis.giant = giant
			polyEvalBis.baby = baby
			polyEvalBis.dummyPolynomialBasis = polyEval.dummyPolynomialBasis
			polyEvalBis.poly = polyEval.poly

			return polyEvalBis.recurse(targetLevel, targetScale, polIn)
		}

		if polIn.lead {
			targetScale.Mul(params.QiFloat64(targetLevel))
		}

		switch polyEval.poly.Coefficients.Value.(type) {
		case [][]*rlwe.Plaintext:
			return polyEval.evaluatePolyFromPolynomialBasisPlaintext(targetScale, targetLevel, &polIn.Polynomial, polyEval.poly)
		case [][]*rlwe.Ciphertext:
			return polyEval.evaluatePolyFromPolynomialBasisCiphertext(targetScale, targetLevel, &polIn.Polynomial, polyEval.poly)
		}
	}

	var nextPower = 1 << polyEval.baby
	for nextPower < (polIn.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := polIn.splitBSGS(nextPower)

	XPow := polyEval.dummyPolynomialBasis.Value[nextPower]

	level := targetLevel

	var currentQi float64
	if polIn.lead {
		currentQi = params.QiFloat64(level)
	} else {
		currentQi = params.QiFloat64(level + 1)
	}

	nextScale := targetScale.CopyNew()
	nextScale.Mul(currentQi)
	nextScale.Div(XPow.Scale)

	res = polyEval.recurse(targetLevel+1, nextScale, coeffsq)

	res.rescale(params, params.DefaultScale())

	res.Scale.Mul(XPow.Scale)

	tmp := polyEval.recurse(res.Level, res.Scale, coeffsr)

	res.Level = utils.MinInt(tmp.Level, res.Level)

	return
}

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasisPlaintext(targetScale rlwe.Scale, level int, polIn, polOut *rlwe.Polynomial) (res *dummyCiphertext) {

	params := polyEval.params
	ecd := polyEval.ecd

	coeffs := polIn.Coefficients.Value.([][]complex128)

	slotsIndex := polIn.SlotsIndex

	X := polyEval.dummyPolynomialBasis.Value

	degree := polIn.Degree()
	nbPoly := len(coeffs)

	//[#poly][coefficients]
	pt := make([]*rlwe.Plaintext, degree+1)

	values := make([]complex128, params.Slots())

	var toEncode bool

	for i := 0; i < nbPoly; i++ {

		c := coeffs[i][0] // [poly][coeff]

		if c != 0 {

			toEncode = true

			for _, j := range slotsIndex[i] {
				values[j] = c
			}
		}
	}

	if toEncode {

		pt[0] = ecd.EncodeNew(values, level, targetScale, params.LogSlots()).Plaintext

		for i := range values {
			values[i] = 0
		}

		toEncode = false
	}

	for i := 1; i < degree+1; i++ {

		for j := 0; j < nbPoly; j++ {

			c := coeffs[j][i] // [poly][coeff]

			if c != 0 {

				toEncode = true

				for _, k := range slotsIndex[j] {
					values[k] = c
				}
			}
		}

		if toEncode {

			scale := targetScale.CopyNew()
			scale.Div(X[i].Scale)

			pt[i] = ecd.EncodeNew(values, level, scale, params.LogSlots()).Plaintext

			for i := range values {
				values[i] = 0
			}

			toEncode = false
		}
	}

	polOut.Coefficients.Value = append([][]*rlwe.Plaintext{pt}, polOut.Coefficients.Value.([][]*rlwe.Plaintext)...)

	return &dummyCiphertext{level, targetScale.CopyNew()}
}

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasisCiphertext(targetScale rlwe.Scale, level int, polIn, polOut *rlwe.Polynomial) (res *dummyCiphertext) {

	params := polyEval.params
	ecd := polyEval.ecd
	enc := polyEval.enc

	coeffs := polIn.Coefficients.Value.([][]complex128)

	slotsIndex := polIn.SlotsIndex

	X := polyEval.dummyPolynomialBasis.Value

	degree := polIn.Degree()
	nbPoly := len(coeffs)

	//[#poly][coefficients]
	ct := make([]*rlwe.Ciphertext, degree+1)

	var values []complex128
	var pt *Plaintext
	if slotsIndex != nil {
		values = make([]complex128, params.Slots())
		pt = NewPlaintext(params, level, &Scale{})
	} else {
		values = make([]complex128, 1)
	}

	var toEncrypt bool

	for i := 0; i < nbPoly; i++ {

		c := coeffs[i][0] // [poly][coeff]

		if c != 0 {

			toEncrypt = true

			if slotsIndex != nil {
				for _, j := range slotsIndex[i] {
					values[j] = c
				}
			} else {
				values[0] = c
			}
		}
	}

	if toEncrypt {
		if slotsIndex != nil {
			pt.Scale().Set(targetScale)
			ecd.Encode(values, pt, params.LogSlots())
			ct[0] = enc.EncryptNew(pt).Ciphertext
		} else {
			tmp := enc.EncryptZeroNew(level, targetScale)
			addConst(params, tmp, values[0], tmp)
			ct[0] = tmp.Ciphertext
		}

		for i := range values {
			values[i] = 0
		}

		toEncrypt = false
	}

	for i := 1; i < degree+1; i++ {

		for j := 0; j < nbPoly; j++ {

			c := coeffs[j][i] // [poly][coeff]

			if c != 0 {

				toEncrypt = true

				if slotsIndex != nil {
					for _, k := range slotsIndex[j] {
						values[k] = c
					}
				} else {
					values[0] = c
				}

			}
		}

		if toEncrypt {

			scale := targetScale.CopyNew()
			scale.Div(X[i].Scale)

			if slotsIndex != nil {
				pt.Scale().Set(scale)
				ecd.Encode(values, pt, params.LogSlots())
				ct[i] = enc.EncryptNew(pt).Ciphertext
			} else {
				tmp := enc.EncryptZeroNew(level, scale)
				addConst(params, tmp, values[0], tmp)
				ct[i] = tmp.Ciphertext
			}

			for i := range values {
				values[i] = 0
			}

			toEncrypt = false
		}
	}

	polOut.Coefficients.Value = append([][]*rlwe.Ciphertext{ct}, polOut.Coefficients.Value.([][]*rlwe.Ciphertext)...)

	return &dummyCiphertext{level, targetScale}
}

type dummyCiphertext struct {
	Level int
	Scale rlwe.Scale
}

func (d *dummyCiphertext) rescale(params Parameters, minScale rlwe.Scale) {
	var nbRescales int
	for d.Level-nbRescales >= 0 && d.Scale.(*Scale).Value/float64(params.Q()[d.Level-nbRescales]) >= minScale.(*Scale).Value/2 {
		d.Scale.Div(float64(params.Q()[d.Level-nbRescales]))
		nbRescales++
	}

	d.Level -= nbRescales
}

type dummyPolynomialBasis struct {
	Value  map[int]*dummyCiphertext
	params Parameters
}

func newDummyPolynomialBasis(params Parameters, ct *dummyCiphertext) (p *dummyPolynomialBasis) {
	p = new(dummyPolynomialBasis)
	p.params = params
	p.Value = make(map[int]*dummyCiphertext)
	p.Value[1] = ct
	return
}

func (p *dummyPolynomialBasis) GenPower(n int, lazy bool, scale rlwe.Scale) {
	if p.Value[n] == nil {
		p.genPower(n, lazy, scale)
		p.Value[n].rescale(p.params, scale)
	}
}

func (p *dummyPolynomialBasis) genPower(n int, lazy bool, scale rlwe.Scale) (err error) {
	if p.Value[n] == nil {

		isPow2 := n&(n-1) == 0

		// Computes the index required to compute the asked ring evaluation
		var a, b int
		if isPow2 {
			a, b = n/2, n/2 //Necessary for optimal depth
		} else {
			// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
			// Maximize the number of odd terms of Chebyshev basis
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)
		}

		// Recurses on the given indexes
		p.genPower(a, lazy, scale)
		p.genPower(b, lazy, scale)

		p.Value[a].rescale(p.params, scale)
		p.Value[b].rescale(p.params, scale)

		// Computes C[n] = C[a]*C[b]
		p.Value[n] = new(dummyCiphertext)
		p.Value[n].Level = utils.MinInt(p.Value[a].Level, p.Value[b].Level)
		p.Value[n].Scale = p.Value[a].Scale.CopyNew()
		p.Value[n].Scale.Mul(p.Value[b].Scale)
		if !(lazy && !isPow2) {
			p.Value[n].rescale(p.params, scale)
		}
	}
	return
}
