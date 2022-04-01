package ring

// DivFloorByLastModulusNTTLvl divides (floored) the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivFloorByLastModulusNTTLvl(level int, p0, pool, p1 *Poly) {
	r.InvNTTSingleLazy(level, p0.Coeffs[level], pool.Coeffs[0])
	for i := 0; i < level; i++ {
		r.NTTSingleLazy(i, pool.Coeffs[0], pool.Coeffs[1])
		// (-x[i] + x[-1]) * -InvQ
		SubVecAndMulScalarMontgomeryTwoQiVec(pool.Coeffs[1], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], r.Modulus[i], r.MredParams[i])
	}
}

// DivFloorByLastModulusLvl divides (floored) the polynomial by its last modulus.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivFloorByLastModulusLvl(level int, p0, p1 *Poly) {
	for i := 0; i < level; i++ {
		SubVecAndMulScalarMontgomeryTwoQiVec(p0.Coeffs[level], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], r.Modulus[i], r.MredParams[i])
	}
}

// DivFloorByLastModulusManyNTTLvl divides (floored) sequentially nbRescales times the polynomial by its last modulus. Input must be in the NTT domain.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivFloorByLastModulusManyNTTLvl(level, nbRescales int, p0, pool, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			CopyValuesLvl(level, p0, p1)
		}

	} else {

		r.InvNTTLvl(level, p0, pool)

		for i := 0; i < nbRescales; i++ {
			r.DivFloorByLastModulusLvl(level-i, pool, pool)
		}

		r.NTTLvl(level-nbRescales, pool, p1)
	}
}

// DivFloorByLastModulusManyLvl divides (floored) sequentially nbRescales times the polynomial by its last modulus.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivFloorByLastModulusManyLvl(level, nbRescales int, p0, pool, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			CopyValuesLvl(level, p0, p1)
		}

	} else {

		if nbRescales > 1 {
			r.DivFloorByLastModulusLvl(level, p0, pool)

			for i := 1; i < nbRescales; i++ {

				if i == nbRescales-1 {
					r.DivFloorByLastModulusLvl(level-i, pool, p1)
				} else {
					r.DivFloorByLastModulusLvl(level-i, pool, pool)
				}
			}

		} else {
			r.DivFloorByLastModulusLvl(level, p0, p1)
		}
	}
}

// DivRoundByLastModulusNTTLvl divides (rounded) the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivRoundByLastModulusNTTLvl(level int, p0, pool, p1 *Poly) {

	r.InvNTTSingleLazy(level, p0.Coeffs[level], pool.Coeffs[level])

	// Center by (p-1)/2
	pj := r.Modulus[level]
	pHalf := (pj - 1) >> 1

	AddScalarVec(pool.Coeffs[level], pool.Coeffs[level], pHalf, pj)

	for i := 0; i < level; i++ {
		qi := r.Modulus[i]
		AddScalarNoModVec(pool.Coeffs[level], pool.Coeffs[i], qi-BRedAdd(pHalf, qi, r.BredParams[i]))
		r.NTTSingleLazy(i, pool.Coeffs[i], pool.Coeffs[i])
		SubVecAndMulScalarMontgomeryTwoQiVec(pool.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], qi, r.MredParams[i])
	}
}

// DivRoundByLastModulusLvl divides (rounded) the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivRoundByLastModulusLvl(level int, p0, p1 *Poly) {

	// Center by (p-1)/2
	pj := r.Modulus[level]
	pHalf := (pj - 1) >> 1

	AddScalarVec(p0.Coeffs[level], p0.Coeffs[level], pHalf, pj)

	for i := 0; i < level; i++ {
		qi := r.Modulus[i]
		AddScalarNoModAndNegTwoQiNoModVec(p0.Coeffs[i], p0.Coeffs[i], qi-BRedAdd(pHalf, qi, r.BredParams[i]), qi)
		AddVecNoModAndMulScalarMontgomeryVec(p0.Coeffs[level], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], qi, r.MredParams[i])
	}
}

// DivRoundByLastModulusManyNTTLvl divides (rounded) sequentially nbRescales times the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivRoundByLastModulusManyNTTLvl(level, nbRescales int, p0, pool, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			CopyValuesLvl(level, p0, p1)
		}

	} else {

		if nbRescales > 1 {

			r.InvNTTLvl(level, p0, pool)
			for i := 0; i < nbRescales; i++ {
				r.DivRoundByLastModulusLvl(level-i, pool, pool)
			}
			r.NTTLvl(level-nbRescales, pool, p1)

		} else {
			r.DivRoundByLastModulusNTTLvl(level, p0, pool, p1)
		}
	}
}

// DivRoundByLastModulusManyLvl divides (rounded) sequentially nbRescales times the polynomial by its last modulus.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivRoundByLastModulusManyLvl(level, nbRescales int, p0, pool, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			CopyValuesLvl(level, p0, p1)
		}

	} else {

		if nbRescales > 1 {

			r.DivRoundByLastModulusLvl(level, p0, pool)

			for i := 1; i < nbRescales; i++ {

				if i == nbRescales-1 {
					r.DivRoundByLastModulusLvl(level-i, pool, p1)
				} else {
					r.DivRoundByLastModulusLvl(level-i, pool, pool)
				}
			}

		} else {
			r.DivRoundByLastModulusLvl(level, p0, p1)
		}
	}
}
