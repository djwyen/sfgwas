package gwas

import (
	"fmt"
	"math"

	"go.dedis.ch/onet/v3/log"

	"strconv"
	"time"

	"github.com/hhcho/sfgwas-private/crypto"

	mpc_core "github.com/hhcho/mpc-core"
	"gonum.org/v1/gonum/mat"

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

// approximates a solution `x` to `Qx = eta` through `n_iterations` of the Conjugated Gradient Method, where Q = (N/sigma2)(sum(LD_shares) + Psi^-1)
func FederatedCGM(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, eta []float64, LDshare crypto.PlainMatrix, Psi_inv crypto.CipherVector, sigma2 *ckks.Ciphertext, N, block, n_iterations int) crypto.CipherVector {
	// TODO eta will have to be a CipherVector actually, as I think the MVN sample added to Qmu is private
	mpcObj := (*mpcObjs)[0]
	pid := mpcObj.GetPid()
	hub_pid := mpcObj.GetHubPid()

	n := len(eta)
	n_slots := cps.GetSlots()
	n_ciphertexts := (n % n_slots) + 1

	
	if n != len(LDshare) {
		// well, technically we're just making an assumption about the CipherVector widths, as we can't directly check them
		panic(fmt.Sprintf("LD Share has dimension %dx%d but eta has length %d", len(LDshare), len(LDshare), n))
	}

	x := crypto.CZeros(cps, n)
	y, _ := crypto.EncryptFloatVector(cps, eta)
	r := y
	delta_new := crypto.InnerProd(cps, r, r)

	for t := 0; t < n_iterations; t++ {
		// z_i := make(crypto.CipherVector, 0)
		z_i := crypto.CZeros(cps, n) // TODO maybe choose a different index than i for the parties, as this is slightly misleading since we use i here for another meaning too
		// recall that the LD matrix is definitionally symmetric, so we can treat the rows as the columns:
		// TODO optimize number of additions to be log(n) by adding intermediate ciphertexts?
		for i := 0; i < n_ciphertexts; i++ {
			for j := 0; j < n_slots; j++ {
				// TODO there's the N_i/N constant thing, but let's just ignore it for now
				index := (i * n_slots) + j
				row_product := crypto.CPMult(cps, y, LDshare[index])
				coordinate := crypto.InnerSum(cps, row_product, n)
				// we must now put the scalar coordinate into the right spot, which crucially also means we must pick out the right Ciphertext to work with
				z_i[i] = crypto.Add(cps, z_i[i], maskAndRotate(cps, coordinate, j))
			}
		}
		
		if pid != hub_pid {
			// TODO send over this party's z_i and continue to the next iteration
			continue
		}
		// the hub will compute psi and do all the other necessary calculations
		// TODO maybe we will need to invert Psi ourselves, but that should be easy to add
		z_star := crypto.CMult(cps, Psi_inv, y)
		// TODO somehow send out and aggregate the LD shares
		u := crypto.CAdd(cps, z_i, z_star)
		crypto.CMultConst(cps, u, N, true)
		crypto.CMultConst(cps, u, cInvert(cps, sigma2), true)

		alpha := crypto.Mult(cps, delta_new, cInvert(cps, crypto.InnerProd(cps, y, u)))
		x = crypto.CAdd(cps, x, crypto.CMultConst(cps, y, alpha, false))
		r = crypto.CSub(cps, r, crypto.CMultConst(cps, u, alpha, false))
		delta_old := delta_new
		delta_new = crypto.InnerProd(cps, r, r)
		beta := crypto.Mult(cps, delta_new, cInvert(cps, delta_old))
		y = crypto.CAdd(cps, r, crypto.CMultConst(cps, y, beta, false))

		// TODO the hub party needs to send the new values of x, y, r, delta_new out to all the other parties too
	}
	return x
}

// TODO implement
// func PreconditionedFederatedCGM
