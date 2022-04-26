package gadget

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// AddPolyToCiphertext takes a plaintext polynomial and a list of Ciphertexts and adds the
// plaintext times the RNS and BIT decomposition to the i-th element of the i-th Ciphertexts.
func AddPolyToCiphertext(pt *ring.Poly, ct []Ciphertext, ringQP ringqp.Ring, logbase2 int, buff *ring.Poly) {

	ringQ := ringQP.RingQ

	levelQ := ct[0].LevelQ()
	levelP := ct[0].LevelP()

	if levelP != -1 {
		ringQ.MulScalarBigintLvl(levelQ, pt, ringQP.RingP.ModulusAtLevel[levelP], buff) // P * pt
	} else {
		levelP = 0
		if pt != buff {
			ring.CopyLvl(levelQ, pt, buff) // 1 * pt
		}
	}

	RNSDecomp := len(ct[0].Value)
	BITDecomp := len(ct[0].Value[0])

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e + (m * P * w^2j) * (q_star * q_tild) mod QP
			//
			// q_prod = prod(q[i*alpha+j])
			// q_star = Q/qprod
			// q_tild = q_star^-1 mod q_prod
			//
			// Therefore : (pt * P * w^2j) * (q_star * q_tild) = pt*P*w^2j mod q[i*alpha+j], else 0
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handle cases where #pj does not divide #qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.Modulus[index]
				p0tmp := buff.Coeffs[index]

				for u, el := range ct {
					p1tmp := el.Value[i][j][u].Q.Coeffs[index]
					for w := 0; w < ringQ.N; w++ {
						p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
					}
				}
			}
		}

		// w^2j
		ringQ.MulScalar(buff, 1<<logbase2, buff)
	}
}

// AddPolyToGadgetMatrix takes a plaintext polynomial and a list of ringqp.Poly adds the
// plaintext times the RNS and BIT decomposition to the list of ringqp.Poly.
func AddPolyToGadgetMatrix(pt *ring.Poly, gm [][]ringqp.Poly, ringQP ringqp.Ring, logbase2 int, buff *ring.Poly) {

	ringQ := ringQP.RingQ

	levelQ := gm[0][0].LevelQ()
	levelP := gm[0][0].LevelP()

	if levelP != -1 {
		ringQ.MulScalarBigintLvl(levelQ, pt, ringQP.RingP.ModulusAtLevel[levelP], buff) // P * pt
	} else {
		levelP = 0
		if pt != buff {
			ring.CopyLvl(levelQ, pt, buff) // 1 * pt
		}
	}

	RNSDecomp := len(gm)
	BITDecomp := len(gm[0])

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e + (m * P * w^2j) * (q_star * q_tild) mod QP
			//
			// q_prod = prod(q[i*alpha+j])
			// q_star = Q/qprod
			// q_tild = q_star^-1 mod q_prod
			//
			// Therefore : (pt * P * w^2j) * (q_star * q_tild) = pt*P*w^2j mod q[i*alpha+j], else 0
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handle cases where #pj does not divide #qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.Modulus[index]
				p0tmp := buff.Coeffs[index]

				p1tmp := gm[i][j].Q.Coeffs[index]
				for w := 0; w < ringQ.N; w++ {
					p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
				}
			}
		}

		// w^2j
		ringQ.MulScalar(buff, 1<<logbase2, buff)
	}
}
