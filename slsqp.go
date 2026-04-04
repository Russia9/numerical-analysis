package numericalanalysis

import "math"

// slsqp.go
// Sequential Least Squares Quadratic Programming

// gradient computes the gradient of f at x using central differences.
func gradient(f func([]float64) float64, x, deltaX []float64) []float64 {
	n := len(x)
	grad := make([]float64, n)
	for i := range n {
		xp := make([]float64, n)
		xm := make([]float64, n)
		copy(xp, x)
		copy(xm, x)
		xp[i] += deltaX[i]
		xm[i] -= deltaX[i]
		grad[i] = (f(xp) - f(xm)) / (2 * deltaX[i])
	}
	return grad
}

// SLSQP minimises f(x) subject to equality and inequality constraints using
// Sequential Least Squares Quadratic Programming.
//
//	min  f(x)
//	s.t. eqConstraints[i](x)   = 0   for all i
//	     ineqConstraints[j](x) >= 0  for all j
//
// Parameters:
//
//	f:               objective function
//	x0:              initial guess
//	deltaX:          finite-difference step for each variable (len == len(x0))
//	eqConstraints:   equality constraint functions g_i(x) = 0
//	ineqConstraints: inequality constraint functions h_j(x) >= 0
//	eps:             convergence tolerance (on the QP step norm and feasibility)
//	maxIter:         maximum number of outer iterations
func SLSQP(
	f func(x []float64) float64,
	x0 []float64,
	deltaX []float64,
	eqConstraints []func(x []float64) float64,
	ineqConstraints []func(x []float64) float64,
	eps float64,
	maxIter int,
) ([]float64, error) {
	n := len(x0)
	if n == 0 || len(deltaX) != n || eps <= 0 || maxIter <= 0 {
		return nil, ErrWrongInput
	}
	for i := range n {
		if deltaX[i] <= 0 {
			return nil, ErrWrongInput
		}
	}

	nEq := len(eqConstraints)
	nIneq := len(ineqConstraints)

	x := make([]float64, n)
	copy(x, x0)

	// BFGS approximation of the Lagrangian Hessian; start from identity
	B := IdentityMatrix(n)

	lagGrad := func(xk []float64, lambdaEq []float64, muIneq []float64) []float64 {
		g := gradient(f, xk, deltaX)
		for i, eq := range eqConstraints {
			ge := gradient(eq, xk, deltaX)
			for j := range n {
				g[j] += lambdaEq[i] * ge[j]
			}
		}
		for j, ineq := range ineqConstraints {
			gh := gradient(ineq, xk, deltaX)
			for k := range n {
				g[k] -= muIneq[j] * gh[k]
			}
		}
		return g
	}

	for range maxIter {
		gradf := gradient(f, x, deltaX)

		// Evaluate constraints and their gradients
		var AEq Matrix
		var bEq []float64
		gVals := make([]float64, nEq)
		if nEq > 0 {
			AEq = make(Matrix, nEq)
			bEq = make([]float64, nEq)
			for i, eq := range eqConstraints {
				gVals[i] = eq(x)
				AEq[i] = gradient(eq, x, deltaX)
				bEq[i] = -gVals[i] // linearised: ∇g·d = -g(x)
			}
		}

		var AIneq Matrix
		var bIneq []float64
		hVals := make([]float64, nIneq)
		if nIneq > 0 {
			AIneq = make(Matrix, nIneq)
			bIneq = make([]float64, nIneq)
			for j, ineq := range ineqConstraints {
				hVals[j] = ineq(x)
				AIneq[j] = gradient(ineq, x, deltaX)
				bIneq[j] = -hVals[j] // linearised: ∇h·d >= -h(x)
			}
		}

		// Solve QP subproblem for search direction d
		d, lambdaEq, muIneq, err := QPActiveSet(B, gradf, AEq, bEq, AIneq, bIneq)
		if err != nil {
			return nil, err
		}

		// Convergence: QP step is tiny and constraints are (nearly) satisfied
		eqFeas := 0.0
		for _, v := range gVals {
			eqFeas += math.Abs(v)
		}
		ineqFeas := 0.0
		for _, v := range hVals {
			if v < 0 {
				ineqFeas += -v
			}
		}
		if Norm(d) < eps && eqFeas < eps && ineqFeas < eps {
			return x, nil
		}

		// Penalty parameter: max of multiplier magnitudes + buffer
		rho := 1.0
		for _, v := range lambdaEq {
			if a := math.Abs(v); a > rho {
				rho = a
			}
		}
		for _, v := range muIneq {
			if v > rho {
				rho = v
			}
		}
		rho += 1.0

		// L1 exact merit function
		phi0 := f(x)
		for _, v := range gVals {
			phi0 += rho * math.Abs(v)
		}
		for _, v := range hVals {
			if v < 0 {
				phi0 -= rho * v
			}
		}

		phi := func(alpha float64) float64 {
			xNew := addVec(x, scaleVec(alpha, d))
			val := f(xNew)
			for _, eq := range eqConstraints {
				val += rho * math.Abs(eq(xNew))
			}
			for _, ineq := range ineqConstraints {
				if h := ineq(xNew); h < 0 {
					val -= rho * h
				}
			}
			return val
		}

		// Directional derivative of L1 merit at α=0 (subgradient)
		dphi0 := dot(gradf, d)
		for _, v := range gVals {
			dphi0 -= rho * math.Abs(v)
		}
		for _, v := range hVals {
			if v < 0 {
				dphi0 -= rho * v // adds rho*|violation|
			}
		}

		// Guarantee descent: if dphi0 >= 0 (can happen with steep constraints),
		// fall back to a small gradient step on f
		alpha := 1.0
		if dphi0 < 0 {
			alpha, _ = LineSearchArmijo(phi, phi0, dphi0, 1.0, 0.1, 50)
		} else {
			alpha = 0.1
		}

		s := scaleVec(alpha, d)
		xNew := addVec(x, s)

		// BFGS update: y = ∇L(x_new, λ) − ∇L(x_old, λ) using the same
		// multipliers for both points so that y ≈ ∇²L·s (proper secant condition).
		// Using different multipliers contaminates y with (λ_new−λ_old)·∇g and
		// causes ||y|| to explode, corrupting B.
		glNew := lagGrad(xNew, lambdaEq, muIneq)
		glOld := lagGrad(x, lambdaEq, muIneq)
		y := subVec(glNew, glOld)
		if Bnew, err := BFGSUpdate(B, s, y); err == nil {
			B = Bnew
		}

		x = xNew
	}

	return nil, ErrDidNotConverge
}
