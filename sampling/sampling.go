package sampling

import (
	"math"

	"github.com/hhcho/sfgwas-private/mpc"
	"github.com/hhcho/sfgwas-private/crypto"

	"github.com/ldsec/lattigo/v2/ckks"

	"gonum.org/v1/gonum/stat/distuv"
	"go.dedis.ch/onet/v3/log"
)

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

// TODO in our application it looks like we will only ever sample from a normal with mean 0, so we can further optimize this
func CSampleNormal(cps *crypto.CryptoParams, mu, sigma2 *ckks.Ciphertext) *ckks.Ciphertext {
	stdNormalDist := distuv.Normal{Mu: 0, Sigma: 1}
	sample := stdNormalDist.Rand()
	rescaledSample := crypto.CMultConst(cps, crypto.CipherVector{cSqrt(cps, sigma2)}, sample, false)[0]
	return crypto.Add(cps, rescaledSample, mu)
}

// sample from the gamma distribution with public shape parameter `k` and private scale parameter `theta`
func CSampleGamma(cps *crypto.CryptoParams, k float64, theta *ckks.Ciphertext) *ckks.Ciphertext {
	// for our application, we apparently will only ever need to sample from gamma with public k and private theta, hence the specificity of this function
	gammaDist := distuv.Gamma{Alpha: k, Beta: 1} // this is a different parameterization than the one we use, but they coincide in the standard form, with rate = scale = 1
	sample := gammaDist.Rand()
	// invoke the property that X ~ G(k, 1) -> thetaX ~ G(k, theta)
	// inv_theta := cInvert(cps, theta)
	return crypto.CMultConst(cps, crypto.CipherVector{theta}, sample, false)[0]
}

// for k >= 1, generates a sample from the Gamma distribution with parameters `k, theta`.
func CSampleGammaMarsagliaTsang(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, k, theta *ckks.Ciphertext) *ckks.Ciphertext {
	mpcObj := (*mpcObjs)[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType().Zero()

	neg_one_third := crypto.EncryptFloat(cps, -1/3)
	d := crypto.Add(cps, k, neg_one_third)
	// addition is so much faster than even plaintext multiplication that we should just use it whenever possible
	two_d := crypto.Add(cps, d, d)
	four_d := crypto.Add(cps, two_d, two_d)
	eight_d := crypto.Add(cps, four_d, four_d)
	nine_d := crypto.Add(cps, eight_d, d)
	// c := crypto.InvSqrtApprox
	c := cInvert(cps, cSqrt(cps, nine_d))

	stdNormalDist := distuv.Normal{Mu: 1, Sigma: 0}
	z := 0.0
	v := crypto.Zero(cps)

	step2:
	for true {
		z = stdNormalDist.Rand()
		parenthetical_term := crypto.CMultConst(cps, crypto.CipherVector{c}, z, false)[0]
		parenthetical_term = crypto.AddConst(cps, parenthetical_term, 1)
		v = crypto.Mult(cps, parenthetical_term, parenthetical_term)
		v = crypto.Mult(cps, v, parenthetical_term) // v = (1 + cZ)^3

		v_secret_shared := mpcObj.CiphertextToSS(cps, mpcObj.GetRType(), v, pid, 1)
		outcome := mpcObj.LessThanPublic(v_secret_shared, rtype.Zero(), false)
		outcome = mpcObj.RevealSymVec(outcome)
		if outcome[0].Uint64() == 0 { // ie that v > 0 TODO I am eliding the difference between > and >= here, but should be fine since it's continuous?
			break step2
		}
	}

	stdUniformDist := distuv.Uniform{Min: 0, Max: 1}
	u := stdUniformDist.Rand()
	// TODO make sure that the variable scope on z in the for loop really still works. Also, the assigment operator might not work between loops?
	threshold_value_1 := 1 - (0.0331 * math.Pow(z, 4))
	// TODO add second threshold value. It involves d/v and therefore must be encrypted -- can we get away with just using the public threshold value condition and having a lower acceptance rate?
	if u < threshold_value_1 {
		return crypto.Mult(cps, d, v)
	}
	goto step2 // TODO rewrite this goto, this is terrible and considered harmful
}

func CSampleInverseGamma(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, k, theta *ckks.Ciphertext) *ckks.Ciphertext {
	inv_theta := cInvert(cps, theta)
	sample := CSampleGammaMarsagliaTsang(cps, mpcObjs, k, inv_theta)
	return cInvert(cps, sample)
}

func CSampleScaledInverseChiSquared(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, v, theta *ckks.Ciphertext) *ckks.Ciphertext {
	sample_k := crypto.CMultConst(cps, crypto.CipherVector{v}, 0.5, false)[0]
	vTheta := crypto.Mult(cps, v, theta)
	sample_theta := crypto.CMultConst(cps, crypto.CipherVector{vTheta}, 0.5, false)[0]
	return CSampleGammaMarsagliaTsang(cps, mpcObjs, sample_k, sample_theta) // Gamma(v/2, vTheta/2)
}

// for a list of `k_values` in which each element is interpreted as a separate scalar Ciphertext,
func CSampleDirichlet(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, k_values crypto.CipherVector) crypto.CipherVector {
	n := len(k_values)
	samples := make(crypto.CipherVector, n)
	one := crypto.EncryptFloat(cps, 1)

	Z := crypto.Zero(cps)
	for i, k_i := range k_values {
		sample := CSampleGammaMarsagliaTsang(cps, mpcObjs, k_i, one)
		samples[i] = sample
		Z = crypto.Add(cps, Z, sample)
	}
	
	Z_inv := cInvert(cps, Z)
	for i, X_i := range samples {
		samples[i] = crypto.Mult(cps, X_i, Z_inv)
	}
	return samples
}

func CSampleGeneralizedInverseGaussian(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, p float64, a, b *ckks.Ciphertext) *ckks.Ciphertext{
	if p == 0.5 {
		mu := crypto.Mult(cps, a, cInvert(cps, b))
		mu = cSqrt(cps, mu)
		sample := CSampleInverseGaussian(cps, mpcObjs, mu, a)
		return cInvert(cps, sample)
	}
	return nil
}

// sample from an inverse Gaussian distribution with parameters mu/lambda using the algorithm of Michael, Schucany, and Haas 1976
func CSampleInverseGaussian(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, mu, lambda *ckks.Ciphertext) *ckks.Ciphertext {
	// TODO bugfix, we're currently yielding impossibly massive samples
	// TODO reply: it's almost certainly due to exhausting Levels, so we need to throw in a bootstrap somewhere
	mpcObj := (*mpcObjs)[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType().Zero()
	
	stdNormalDist := distuv.Normal{Mu: 1, Sigma: 0}
	nu := stdNormalDist.Rand()
	y := math.Pow(nu, 2) // strictly speaking we could have sampled from a chi squared distribution

	mu_y := crypto.CMultConst(cps, crypto.CipherVector{mu}, y, false)[0]
	mu2_y2 := crypto.Mult(cps, mu_y, mu_y)
	// addition is so much cheaper than even plaintext constant multiplication that we should just use it whenever possible
	two_lambda := crypto.Add(cps, lambda, lambda)
	inv_two_lambda := cInvert(cps, two_lambda)
	four_lambda := crypto.Add(cps, two_lambda, two_lambda)
	four_lambda_mu_y := crypto.Mult(cps, four_lambda, mu_y)
	sqrt_term := cSqrt(cps, crypto.Add(cps, four_lambda_mu_y, mu2_y2))

	divided_term := crypto.CSub(cps, crypto.CipherVector{mu_y}, crypto.CipherVector{sqrt_term})[0]
	divided_term = crypto.Mult(cps, divided_term, inv_two_lambda)
	final_term := crypto.CAddConst(cps, crypto.CipherVector{divided_term}, 1)[0]
	x := crypto.Mult(cps, mu, final_term)

	log.LLvl1("Sampling finished, now checking boundary condition")

	threshold_value := crypto.Add(cps, x, mu)
	threshold_value = cInvert(cps, threshold_value)
	threshold_value = crypto.Mult(cps, mu, threshold_value)

	// TODO replace below with dummy decrypt comparison

	// we now need to return a certain value with probability (mu/(mu + x)), which cannot be done homomorphically (or at least not that easily)
	// so we must now convert the homomorphically encrypted value into a secret-shared one such that we can invoke an MPC protocol for evaluating less-than
	secret_shared := mpcObj.CiphertextToSS(cps, mpcObj.GetRType(), threshold_value, pid, 1)
	// TODO in examples I see source PID set to -1? Is this a meaningful convention?
	stdUniformDist := distuv.Uniform{Min: 0, Max: 1}
	z := stdUniformDist.Rand()
	z_ring := rtype.FromFloat64(z, mpcObj.GetFracBits())
	// as an aside, there might be a way to parallelize many samples if we can somehow do a LessThanPublic where it's a vector of public values to compare against
	outcome := mpcObj.LessThanPublic(secret_shared, z_ring, true)
	outcome = mpcObj.RevealSymVec(outcome)

	
	log.LLvl1("Continuing...")

	if outcome[0].Uint64() == 0 { // ie that z <= threshold_value
		return x
	}
	// otherwise we must manipulate further
	mu2 := crypto.Mult(cps, mu, mu) // potentially cheaper to divide mu2_y2 by plaintext y2?
	return crypto.Mult(cps, mu2, cInvert(cps, x))
}
