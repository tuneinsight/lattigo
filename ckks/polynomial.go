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

func NewPolynomial(polynomialBasis rlwe.PolynomialBasisType, coeffs interface{}, slotsIndex [][]int) (*Polynomial, error) {
	pol, err := rlwe.NewPolynomial(polynomialBasis, coeffs, slotsIndex)
	if err != nil {
		return nil, err
	}
	pol.Scale = &Scale{}
	return &Polynomial{pol}, nil
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

func (p *Polynomial) Encode(ecd Encoder, level int, inputScale, outputScale rlwe.Scale) (*Polynomial, error) {

	params := ecd.(*encoderComplex128).params

	ptPoly := rlwe.Polynomial{}

	switch p.Coefficients.Value.(type) {
	case []float64, []complex128, [][]float64, [][]complex128:

		ptPoly.Coefficients = rlwe.Coefficients{Value: [][]rlwe.Operand{}}
		ptPoly.Basis = p.Basis
		ptPoly.Odd = p.Odd
		ptPoly.Even = p.Even
		ptPoly.SlotsIndex = p.SlotsIndex
		ptPoly.Scale = &Scale{}

		embed := func(values interface{}, level int, scale rlwe.Scale) (op rlwe.Operand) {

			v := values.([]complex128)

			var pt *rlwe.Plaintext
			if len(v) == 1 {
				pt = rlwe.NewPlaintext(params.Parameters, level)
				addConst(params, &Ciphertext{pt.El()}, v[0], &Ciphertext{pt.El()})
			} else {
				pt = ecd.EncodeNew(v, level, scale, params.LogSlots()).Plaintext
			}

			pt.Scale = scale.CopyNew()

			return pt
		}

		getScaledBSGSCoefficients(params, embed, level, inputScale, p.Polynomial, outputScale, &ptPoly)

	default:
		return nil, fmt.Errorf("Polynomial.Encode(*): underlying polynomial is already encoded or encrypted")
	}

	return &Polynomial{ptPoly}, nil
}

func (p *Polynomial) Encrypt(ecd Encoder, enc Encryptor, level int, inputScale, outputScale rlwe.Scale) (*Polynomial, error) {

	params := ecd.(*encoderComplex128).params

	ctPoly := rlwe.Polynomial{}

	switch p.Coefficients.Value.(type) {
	case []float64, []complex128, [][]float64, [][]complex128:

		ctPoly.Coefficients = rlwe.Coefficients{Value: [][]rlwe.Operand{}}
		ctPoly.Basis = p.Basis
		ctPoly.Odd = p.Odd
		ctPoly.Even = p.Even
		ctPoly.SlotsIndex = p.SlotsIndex
		ctPoly.Scale = &Scale{}

		buff := params.RingQ().NewPolyLvl(level)

		embed := func(values interface{}, level int, scale rlwe.Scale) (op rlwe.Operand) {

			v := values.([]complex128)

			var ct *Ciphertext
			if len(v) == 1 {
				ct = enc.EncryptZeroNew(level, scale)
				addConst(params, ct, v[0], ct)
			} else {
				pt := NewPlaintextAtLevelFromPoly(level, buff)
				pt.Plaintext.Scale = scale.CopyNew()
				ecd.Encode(values, pt, params.LogSlots())
				ct = enc.EncryptNew(pt)
			}

			return ct.Ciphertext
		}

		getScaledBSGSCoefficients(params, embed, level, inputScale, p.Polynomial, outputScale, &ctPoly)

	default:
		return nil, fmt.Errorf("Polynomial.Encrypt(*): underlying polynomial is already encrypted")
	}

	return &Polynomial{ctPoly}, nil
}

type dummyPolynomialEvaluator struct {
	dummyPolynomialBasis
	embed  func(values interface{}, level int, scale rlwe.Scale) (op rlwe.Operand)
	params Parameters
	giant  int
	baby   int
	poly   *rlwe.Polynomial
}

type evalPoly struct {
	rlwe.Polynomial
	lead   bool
	maxDeg int
}

func (p *evalPoly) splitBSGS(split int) (polyq, polyr evalPoly) {

	polyq = evalPoly{}
	polyr = evalPoly{}

	polyq.Polynomial, polyr.Polynomial = p.Polynomial.SplitBSGS(split)

	polyq.lead = p.lead
	polyq.maxDeg = p.maxDeg

	if p.maxDeg == p.Degree() {
		polyr.maxDeg = split - 1
	} else {
		polyr.maxDeg = p.maxDeg - (p.Degree() - split + 1)
	}

	return
}

func getScaledBSGSCoefficients(params Parameters, embed func(values interface{}, level int, scale rlwe.Scale) (op rlwe.Operand), level int, scale rlwe.Scale, polIn rlwe.Polynomial, targetScale rlwe.Scale, polOut *rlwe.Polynomial) {

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
	polyEval.dummyPolynomialBasis = *dummbpb
	polyEval.giant = giant
	polyEval.baby = baby
	polyEval.poly = polOut
	polyEval.embed = embed

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
			polyEvalBis.giant = giant
			polyEvalBis.baby = baby
			polyEvalBis.dummyPolynomialBasis = polyEval.dummyPolynomialBasis
			polyEvalBis.poly = polyEval.poly
			polyEvalBis.embed = polyEval.embed

			return polyEvalBis.recurse(targetLevel, targetScale, polIn)
		}

		if polIn.lead {
			targetScale.Mul(params.QiFloat64(targetLevel))
		}

		return polyEval.embedPoly(targetScale, targetLevel, &polIn.Polynomial, polyEval.poly)
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

func (polyEval *dummyPolynomialEvaluator) embedPoly(targetScale rlwe.Scale, level int, polIn, polOut *rlwe.Polynomial) (res *dummyCiphertext) {

	params := polyEval.params
	embed := polyEval.embed

	coeffs := polIn.Coefficients.Value.([][]complex128)

	slotsIndex := polIn.SlotsIndex

	X := polyEval.dummyPolynomialBasis.Value

	degree := polIn.Degree()
	nbPoly := len(coeffs)

	//[#poly][coefficients]
	ops := make([]rlwe.Operand, degree+1)

	var values []complex128
	if slotsIndex != nil {
		values = make([]complex128, params.Slots())
	} else {
		values = make([]complex128, 1)
	}

	var toEmbed bool

	for i := 0; i < nbPoly; i++ {

		c := coeffs[i][0] // [poly][coeff]

		if c != 0 {

			toEmbed = true

			if slotsIndex != nil {
				for _, j := range slotsIndex[i] {
					values[j] = c
				}
			} else {
				values[0] = c
			}
		}
	}

	if toEmbed {

		ops[0] = embed(values, level, targetScale)

		for i := range values {
			values[i] = 0
		}

		toEmbed = false
	}

	for i := 1; i < degree+1; i++ {

		for j := 0; j < nbPoly; j++ {

			c := coeffs[j][i] // [poly][coeff]

			if c != 0 {

				toEmbed = true

				if slotsIndex != nil {
					for _, k := range slotsIndex[j] {
						values[k] = c
					}
				} else {
					values[0] = c
				}
			}
		}

		if toEmbed {

			scale := targetScale.CopyNew()
			scale.Div(X[i].Scale)

			ops[i] = embed(values, level, scale)

			for i := range values {
				values[i] = 0
			}

			toEmbed = false
		}
	}

	polOut.Coefficients.Value = append([][]rlwe.Operand{ops}, polOut.Coefficients.Value.([][]rlwe.Operand)...)

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
