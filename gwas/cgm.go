package gwas

import (
	"fmt"
	"math"

	"go.dedis.ch/onet/v3/log"

	"github.com/hhcho/sfgwas-private/crypto"

	// "github.com/hhcho/sfgwas-private/mpc"

	"github.com/ldsec/lattigo/v2/ckks"
)

// masks then rotates a one-element ciphertext such that it becomes a singleton vector whose only nonzero entry is `nrot` places to the right of 0
func maskAndRotate(cps *crypto.CryptoParams, ct *ckks.Ciphertext, nrot int) *ckks.Ciphertext {
	// since the ciphertext is one-element, we know its one element of relevance is at slot 0
	return crypto.RotateRight(cps, crypto.Mask(cps, ct, 0, false), nrot)
}

// TODO the below functions have equivalents in crypto/basics.go, but I have my own dummy implementations pending understanding how to choose the Intervals/the possibility such functions will be rewritten

// inverts a ciphertext representing a single scalar, TODO currently a dummy implementation
func cInvert(cps *crypto.CryptoParams, x *ckks.Ciphertext) *ckks.Ciphertext {
	// TODO pending figuring out how to determine "reasonable" values of A/B/Degree for the Chebyshev approximation
	// intv := crypto.IntervalApprox{
	// 	A: ,
	// 	B: ,
	// 	Degree: ,
	// 	Iter: ,
	// 	InverseNew: ,
	// }
	// return crypto.CInverse(cps, crypto.CipherVector{x}, intv)[0]
	return crypto.EncryptFloat(cps, 1/crypto.DecryptFloat(cps, x))
}

// returns the square root of a ciphertext representing a single scalar, TODO currently a dummy implementation
func cSqrt(cps *crypto.CryptoParams, x *ckks.Ciphertext) *ckks.Ciphertext {
	return crypto.EncryptFloat(cps, math.Sqrt(crypto.DecryptFloat(cps, x)))
}

// takes in a scalar represented as just a slot-0 Ciphertext and returns a Ciphertext in which the all slots are that Ciphertext
func floodScalar(cps *crypto.CryptoParams, k *ckks.Ciphertext) *ckks.Ciphertext {
	return crypto.Rebalance(cps, crypto.CMultConstRescale(cps, crypto.CipherVector{k}, float64(cps.GetSlots()), false)[0])
}

// approximates a solution `x` to `Qx = eta` through `n_iterations` of the Conjugated Gradient Method, where Q = (N/sigma2)(sum(LD_shares) + Psi^-1)
func FederatedCGM(cps *crypto.CryptoParams, eta []float64, LDshare crypto.PlainMatrix, Psi_inv crypto.CipherVector, sigma2 *ckks.Ciphertext, N, block, n_iterations int) crypto.CipherVector {
	// TODO bugfix this part for the comparisons
	// TODO eta will have to be a CipherVector actually, as I think the MVN sample added to Qmu is private
	// mpcObj := (*mpcObjs)[0]
	// pid := mpcObj.GetPid()
	// hub_pid := mpcObj.GetHubPid()
	pid := 0
	hub_pid := 0

	n := len(eta)
	n_slots := cps.GetSlots()
	// n_ciphertexts := (n % n_slots) + 1
	log.LLvl1("Value of eta:", eta)
	log.LLvl1("Value of LDshare:", crypto.DecodeFloatVector(cps, LDshare[0])[:2])
	log.LLvl1("Value of LDshare:", crypto.DecodeFloatVector(cps, LDshare[1])[:2])
	log.LLvl1("Value of Psi_inv:", crypto.DecryptFloatVector(cps, Psi_inv, n))

	
	if n != len(LDshare) {
		// well, technically we're just making an assumption about the CipherVector widths, as we can't directly check them
		panic(fmt.Sprintf("LD Share has dimension %dx%d but eta has length %d", len(LDshare), len(LDshare), n))
	}

	// x := crypto.CZeros(cps, 1) // TODO a lot of bugs stem from this — the int parameter is the # of ciphertexts, not the number of encrypted items
	x, _ := crypto.EncryptFloatVector(cps, []float64{2.0, 1.0}) // TODO a lot of bugs stem from this — the int parameter is the # of ciphertexts, not the number of encrypted items
	y, _ := crypto.EncryptFloatVector(cps, eta)

	// initialization aggregation of values
	Ax_plus := crypto.CZeros(cps, 1)
	// recall that the LD matrix is definitionally symmetric, so we can treat the rows as the columns:
	for index := 0; index < n; index++ {
		row_product := crypto.CPMult(cps, x, LDshare[index])
		// TODO if we change to InnerSumAll we can just have all the slots have the same value and expose the only relevant one
		coordinate := crypto.InnerSum(cps, row_product, n)
		// we must now put the scalar coordinate into the right spot, which crucially also means we must pick out the right Ciphertext to work with
		i := 0
		Ax_plus[i] = crypto.Add(cps, Ax_plus[i], maskAndRotate(cps, coordinate, index % n_slots))
	}
	// for i := 0; i < n_ciphertexts; i++ {
	// 	for j := 0; j < n_slots; j++ {
	// 		// TODO there's the N_i/N constant thing, but let's just ignore it for now
	// 		index := (i * n_slots) + j
	// 		row_product := crypto.CPMult(cps, y, LDshare[index])
	// 		coordinate := crypto.InnerSum(cps, row_product, n)
	// 		// we must now put the scalar coordinate into the right spot, which crucially also means we must pick out the right Ciphertext to work with
	// 		z_plus[i] = crypto.Add(cps, z_plus[i], maskAndRotate(cps, coordinate, j))
	// 	}
	// }
	// the hub will compute psi and do all the other necessary calculations
	// TODO maybe we will need to invert Psi ourselves, but that should be easy to add
	Ax_star := crypto.CMult(cps, Psi_inv, x)
	log.LLvl1("Value of Ax_plus:", crypto.DecryptFloatVector(cps, Ax_plus, n), "level", Ax_plus[0].Level())
	log.LLvl1("Value of Ax_star:", crypto.DecryptFloatVector(cps, Ax_star, n), "level", Ax_star[0].Level())
	// TODO somehow send out and aggregate the LD shares
	Ax := crypto.CAdd(cps, Ax_plus, Ax_star)
	log.LLvl1("Value of Ax:", crypto.DecryptFloatVector(cps, Ax, n), "level", Ax[0].Level())
	// crypto.CMultConst(cps, u, N, true)
	// log.LLvl1("Value of u:", crypto.DecryptFloatVector(cps, u, n))
	// crypto.CMultConst(cps, u, cInvert(cps, sigma2), true)
	// log.LLvl1("Value of u:", crypto.DecryptFloatVector(cps, u, n))
	y = crypto.CSub(cps, y, Ax)
	r := y
	delta_new := crypto.InnerProd(cps, r, r)

	log.LLvl1("Value of y:", crypto.DecryptFloatVector(cps, y, n), "level", y[0].Level())
	log.LLvl1("Value of r:", crypto.DecryptFloatVector(cps, r, n), "level", r[0].Level())
	log.LLvl1("Value of x:", crypto.DecryptFloatVector(cps, x, n), "level", x[0].Level())
	log.LLvl1("Value of delta_new:", crypto.DecryptFloat(cps, delta_new), "level", delta_new.Level())
	for t := 0; t < n_iterations; t++ {
		log.LLvl1("=== iteration", t, "===")

		// TODO bootstrap everything as necessary
		x = crypto.LevelTest(x, cps, 5, "foo", "bar")
		x = crypto.CRescale(cps, x)
		y = crypto.LevelTest(y, cps, 5, "foo", "bar")
		y = crypto.CRescale(cps, y)
		r = crypto.LevelTest(r, cps, 5, "foo", "bar")
		r = crypto.CRescale(cps, r)
		delta_new = crypto.LevelTest(crypto.CipherVector{delta_new}, cps, 5, "foo", "bar")[0]
		delta_new = crypto.CRescale(cps, crypto.CipherVector{delta_new})[0]

		// z_i := make(crypto.CipherVector, 0)
		z_plus := crypto.CZeros(cps, 1)
		// recall that the LD matrix is definitionally symmetric, so we can treat the rows as the columns:
		// TODO optimize number of additions to be log(n) by adding intermediate ciphertexts?
		for index := 0; index < n; index++ {
			row_product := crypto.CPMult(cps, y, LDshare[index])
			// TODO if we change to InnerSumAll we can just have all the slots have the same value and expose the only relevant one
			coordinate := crypto.InnerSum(cps, row_product, n)
			// we must now put the scalar coordinate into the right spot, which crucially also means we must pick out the right Ciphertext to work with
			i := 0
			z_plus[i] = crypto.Add(cps, z_plus[i], maskAndRotate(cps, coordinate, index % n_slots))
		}
		// for i := 0; i < n_ciphertexts; i++ {
		// 	for j := 0; j < n_slots; j++ {
		// 		// TODO there's the N_i/N constant thing, but let's just ignore it for now
		// 		index := (i * n_slots) + j
		// 		row_product := crypto.CPMult(cps, y, LDshare[index])
		// 		coordinate := crypto.InnerSum(cps, row_product, n)
		// 		// we must now put the scalar coordinate into the right spot, which crucially also means we must pick out the right Ciphertext to work with
		// 		z_plus[i] = crypto.Add(cps, z_plus[i], maskAndRotate(cps, coordinate, j))
		// 	}
		// }
		
		if pid != hub_pid {
			// TODO send over this party's z_plus and continue to the next iteration
			continue
		}
		// the hub will compute psi and do all the other necessary calculations
		// TODO maybe we will need to invert Psi ourselves, but that should be easy to add
		z_plus = crypto.CRescale(cps, z_plus)

		z_star := crypto.CMult(cps, Psi_inv, y)
		z_star = crypto.CRescale(cps, z_star)
		log.LLvl1("Value of z_plus:", crypto.DecryptFloatVector(cps, z_plus, n), "level", z_plus[0].Level())
		log.LLvl1("Value of z_star:", crypto.DecryptFloatVector(cps, z_star, n), "level", z_star[0].Level())
		// TODO somehow send out and aggregate the LD shares
		u := crypto.CAdd(cps, z_plus, z_star)
		u = crypto.CRescale(cps, u)
		log.LLvl1("Value of u:", crypto.DecryptFloatVector(cps, u, n), "level", u[0].Level())
		// crypto.CMultConst(cps, u, N, true)
		// log.LLvl1("Value of u:", crypto.DecryptFloatVector(cps, u, n))
		// crypto.CMultConst(cps, u, cInvert(cps, sigma2), true)
		// log.LLvl1("Value of u:", crypto.DecryptFloatVector(cps, u, n))

		denom := crypto.InnerProd(cps, y, u)
		denom = crypto.CRescale(cps, crypto.CipherVector{denom})[0]

		alpha := crypto.Mult(cps, delta_new, cInvert(cps, denom))
		alpha = floodScalar(cps, alpha)
		alpha = crypto.CRescale(cps, crypto.CipherVector{alpha})[0]
		log.LLvl1("Value of alpha:", crypto.DecryptFloat(cps, alpha), "level", alpha.Level())
		// tmp1 := crypto.CMultConst(cps, y, alpha, false)
		alpha_y := crypto.CMultScalar(cps, y, alpha)
		alpha_y = crypto.CRescale(cps, alpha_y)
		log.LLvl1("Value of alpha_y:", crypto.DecryptFloatVector(cps, alpha_y, n), "level", alpha_y[0].Level())
		// alpha_y1 := crypto.CRescale(cps, alpha_y)
		// log.LLvl1("Value of alpha_y1:", crypto.DecryptFloatVector(cps, alpha_y1, n), "level", alpha_y1[0].Level())
		// alpha_y11, _ := crypto.EncryptFloatVector(cps, crypto.DecryptFloatVector(cps, alpha_y, n))
		// log.LLvl1("Value of alpha_y11:", crypto.DecryptFloatVector(cps, alpha_y11, n), "level", alpha_y11[0].Level())

		// foo, _ := crypto.EncryptFloatVector(cps, []float64{0.25, 0.25})
		// tmp4 := crypto.CMultScalar(cps, y, foo[0]) // debug stuff
		x = crypto.CAdd(cps, x, alpha_y)
		x = crypto.CRescale(cps, x)
		log.LLvl1("Value of x:", crypto.DecryptFloatVector(cps, x, n), "level", x[0].Level())
		// tmp2 := crypto.CMultConst(cps, u, alpha, false)
		alpha_u := crypto.CMultScalar(cps, u, alpha)
		alpha_u = crypto.CRescale(cps, alpha_u)
		log.LLvl1("Value of alpha_u:", crypto.DecryptFloatVector(cps, alpha_u, n), "level", alpha_u[0].Level())
		r = crypto.CSub(cps, r, alpha_u)
		r = crypto.CRescale(cps, r)
		log.LLvl1("Value of r:", crypto.DecryptFloatVector(cps, r, n), "level", r[0].Level())
		delta_old := delta_new
		log.LLvl1("Value of delta_old:", crypto.DecryptFloat(cps, delta_old), "level", delta_old.Level())

		// delta_new = crypto.CRescale(cps, crypto.CipherVector{crypto.InnerProd(cps, r, r)})[0]
		delta_new = crypto.InnerProd(cps, r, r)
		// I think this is the one point where bootstrapping is necessary?
		delta_new = crypto.LevelTest(crypto.CipherVector{delta_new}, cps, 4, "foo", "bar")[0]
		delta_new = crypto.CRescale(cps, crypto.CipherVector{delta_new})[0]
		log.LLvl1("Value of delta_new:", crypto.DecryptFloat(cps, delta_new), "level", delta_new.Level())

		beta := crypto.Mult(cps, delta_new, cInvert(cps, delta_old))
		beta = floodScalar(cps, beta)
		beta = crypto.CRescale(cps, crypto.CipherVector{beta})[0]
		// TODO bugfix line below:
		beta = crypto.EncryptFloat(cps, 0.0088)
		log.LLvl1("Value of beta:", crypto.DecryptFloat(cps, beta), "level", beta.Level())
		// tmp3 := crypto.CMultConst(cps, y, beta, false)
		beta_y := crypto.CMultScalar(cps, y, beta)
		beta_y = crypto.CRescale(cps, beta_y)
		log.LLvl1("Value of beta_y:", crypto.DecryptFloatVector(cps, beta_y, n), "level", beta_y[0].Level())
		y = crypto.CAdd(cps, r, beta_y)
		y = crypto.CRescale(cps, y)
		log.LLvl1("Value of y:", crypto.DecryptFloatVector(cps, y, n), "level", y[0].Level())

		// TODO the hub party needs to send the new values of x, y, r, delta_new out to all the other parties too
	}
	return x
}

// TODO implement
// func PreconditionedFederatedCGM
