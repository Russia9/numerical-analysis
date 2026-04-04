package numericalanalysis

import "math"

// linesearch.go
// Line search methods

// LineSearchArmijo performs backtracking line search satisfying the Armijo
// (sufficient decrease) condition:
//
//	φ(α) ≤ φ(0) + c1·α·φ'(0)
//
// phi:   merit function φ(α)
// phi0:  φ(0)
// dphi0: φ'(0) — must be negative (descent direction)
// alpha0: initial step size
// c1:    sufficient decrease constant ∈ (0, 1)
// maxIter: maximum backtracking steps
func LineSearchArmijo(phi func(alpha float64) float64, phi0, dphi0, alpha0, c1 float64, maxIter int) (float64, error) {
	if dphi0 >= 0 || c1 <= 0 || c1 >= 1 || alpha0 <= 0 || maxIter <= 0 {
		return 0, ErrWrongInput
	}

	alpha := alpha0
	for range maxIter {
		if phi(alpha) <= phi0+c1*alpha*dphi0 {
			return alpha, nil
		}
		alpha *= 0.5
	}
	return alpha, ErrDidNotConverge
}

// LineSearchWolfe finds a step size satisfying the strong Wolfe conditions:
//
//  1. Sufficient decrease:    φ(α) ≤ φ(0) + c1·α·φ'(0)
//  2. Curvature condition:    |φ'(α)| ≤ c2·|φ'(0)|
//
// phi:   merit function φ(α)
// dphi:  derivative φ'(α)
// phi0, dphi0: values at α=0 (dphi0 must be negative)
// alpha0: initial trial step (should be 1.0 for Newton-like methods)
// c1 ∈ (0,1), c2 ∈ (c1,1)
// maxIter: maximum total function evaluations
func LineSearchWolfe(phi, dphi func(alpha float64) float64, phi0, dphi0, alpha0, c1, c2 float64, maxIter int) (float64, error) {
	if dphi0 >= 0 || c1 <= 0 || c1 >= c2 || c2 >= 1 || alpha0 <= 0 || maxIter <= 0 {
		return 0, ErrWrongInput
	}

	zoom := func(lo, hi, phiLo float64) (float64, error) {
		for range maxIter {
			// Bisect
			alpha := (lo + hi) / 2
			phiA := phi(alpha)
			if phiA > phi0+c1*alpha*dphi0 || phiA >= phiLo {
				hi = alpha
			} else {
				dphiA := dphi(alpha)
				if math.Abs(dphiA) <= -c2*dphi0 {
					return alpha, nil
				}
				if dphiA*(hi-lo) >= 0 {
					hi = lo
				}
				lo = alpha
				phiLo = phiA
			}
		}
		return (lo + hi) / 2, ErrDidNotConverge
	}

	prevAlpha := 0.0
	prevPhi := phi0
	alpha := alpha0

	for range maxIter {
		phiA := phi(alpha)
		if phiA > phi0+c1*alpha*dphi0 || (prevAlpha > 0 && phiA >= prevPhi) {
			return zoom(prevAlpha, alpha, prevPhi)
		}
		dphiA := dphi(alpha)
		if math.Abs(dphiA) <= -c2*dphi0 {
			return alpha, nil
		}
		if dphiA >= 0 {
			return zoom(alpha, prevAlpha, phiA)
		}
		prevAlpha = alpha
		prevPhi = phiA
		alpha *= 2
	}
	return alpha, ErrDidNotConverge
}
