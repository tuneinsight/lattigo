package ring

// DivFloorByLastModulusNTT divides (floored) the polynomial by its last modulus.
// The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivFloorByLastModulusNTT(p0, buff, p1 *Poly) {

	level := r.level

	r.InvNTTSingleLazy(r.Tables[level], p0.Coeffs[level], buff.Coeffs[0])
	for i, table := range r.Tables[:level] {
		r.NTTSingleLazy(table, buff.Coeffs[0], buff.Coeffs[1])
		// (-x[i] + x[-1]) * -InvQ
		SubVecAndMulScalarMontgomeryTwoQiVec(buff.Coeffs[1], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[r.level-1][i], table.Modulus, table.MRedParams)
	}
}

// DivFloorByLastModulus divides (floored) the polynomial by its last modulus.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivFloorByLastModulus(p0, p1 *Poly) {

	level := r.level

	for i := 0; i < level; i++ {
		SubVecAndMulScalarMontgomeryTwoQiVec(p0.Coeffs[level], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], r.Tables[i].Modulus, r.Tables[i].MRedParams)
	}
}

// DivFloorByLastModulusManyNTT divides (floored) sequentially nbRescales times the polynomial by its last modulus. Input must be in the NTT domain.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivFloorByLastModulusManyNTT(nbRescales int, p0, buff, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			copy(p1.Buff, p0.Buff)
		}

	} else {

		rCpy := r.AtLevel(r.Level())

		rCpy.InvNTT(p0, buff)

		for i := 0; i < nbRescales; i++ {
			rCpy.DivFloorByLastModulus(buff, buff)
			rCpy = rCpy.AtLevel(rCpy.Level() - 1)
		}

		rCpy.NTT(buff, p1)
	}
}

// DivFloorByLastModulusMany divides (floored) sequentially nbRescales times the polynomial by its last modulus.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivFloorByLastModulusMany(nbRescales int, p0, buff, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			copy(p1.Buff, p0.Buff)
		}

	} else {

		if nbRescales > 1 {

			rCpy := r.AtLevel(r.Level())

			rCpy.DivFloorByLastModulus(p0, buff)
			rCpy = rCpy.AtLevel(rCpy.Level() - 1)

			for i := 1; i < nbRescales; i++ {

				if i == nbRescales-1 {
					rCpy.DivFloorByLastModulus(buff, p1)
				} else {
					rCpy.DivFloorByLastModulus(buff, buff)
				}

				rCpy = rCpy.AtLevel(rCpy.Level() - 1)
			}

		} else {
			r.DivFloorByLastModulus(p0, p1)
		}
	}
}

// DivRoundByLastModulusNTT divides (rounded) the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivRoundByLastModulusNTT(p0, buff, p1 *Poly) {

	level := r.level

	r.InvNTTSingleLazy(r.Tables[level], p0.Coeffs[level], buff.Coeffs[level])

	// Center by (p-1)/2
	pj := r.Tables[level].Modulus
	pHalf := (pj - 1) >> 1

	AddScalarVec(buff.Coeffs[level], buff.Coeffs[level], pHalf, pj)

	for i, table := range r.Tables[:level] {
		qi := table.Modulus
		AddScalarNoModVec(buff.Coeffs[level], buff.Coeffs[i], qi-BRedAdd(pHalf, qi, table.BRedParams))
		r.NTTSingleLazy(table, buff.Coeffs[i], buff.Coeffs[i])
		SubVecAndMulScalarMontgomeryTwoQiVec(buff.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], qi, table.MRedParams)
	}
}

// DivRoundByLastModulus divides (rounded) the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivRoundByLastModulus(p0, p1 *Poly) {

	level := r.level

	// Center by (p-1)/2
	pj := r.Tables[level].Modulus
	pHalf := (pj - 1) >> 1

	AddScalarVec(p0.Coeffs[level], p0.Coeffs[level], pHalf, pj)

	for i := 0; i < level; i++ {
		table := r.Tables[i]
		qi := table.Modulus
		AddScalarNoModAndNegTwoQiNoModVec(p0.Coeffs[i], p0.Coeffs[i], qi-BRedAdd(pHalf, qi, table.BRedParams), qi)
		AddVecNoModAndMulScalarMontgomeryVec(p0.Coeffs[level], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], qi, table.MRedParams)
	}
}

// DivRoundByLastModulusManyNTT divides (rounded) sequentially nbRescales times the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivRoundByLastModulusManyNTT(level, nbRescales int, p0, buff, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			copy(p1.Buff, p0.Buff)
		}

	} else {

		if nbRescales > 1 {

			rCpy := r.AtLevel(r.Level())

			rCpy.InvNTT(p0, buff)
			for i := 0; i < nbRescales; i++ {
				rCpy.DivRoundByLastModulus(buff, buff)
				rCpy = rCpy.AtLevel(rCpy.Level() - 1)
			}

			rCpy.NTT(buff, p1)

		} else {
			r.DivRoundByLastModulusNTT(p0, buff, p1)
		}
	}
}

// DivRoundByLastModulusMany divides (rounded) sequentially nbRescales times the polynomial by its last modulus.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivRoundByLastModulusMany(nbRescales int, p0, buff, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			copy(p1.Buff, p0.Buff)
		}

	} else {

		if nbRescales > 1 {

			rCpy := r.AtLevel(r.Level())

			rCpy.DivRoundByLastModulus(p0, buff)
			rCpy = rCpy.AtLevel(rCpy.Level() - 1)

			for i := 1; i < nbRescales; i++ {

				if i == nbRescales-1 {
					rCpy.DivRoundByLastModulus(buff, p1)
				} else {
					rCpy.DivRoundByLastModulus(buff, buff)
				}

				rCpy = rCpy.AtLevel(rCpy.Level() - 1)
			}

		} else {
			r.DivRoundByLastModulus(p0, p1)
		}
	}
}
