package numericalanalysis

// qp.go
// Quadratic Programming solvers

// QPEqualityConstrained solves the equality-constrained QP:
//
//	min  ½ xᵀHx + cᵀx
//	s.t. Ax = b
//
// via the KKT system:
//
//	[H  Aᵀ] [x]   [-c]
//	[A  0 ] [λ] = [ b]
//
// Returns x (primal), lambda (dual/multipliers), error.
func QPEqualityConstrained(H Matrix, c []float64, A Matrix, b []float64) (x []float64, lambda []float64, err error) {
	n := len(c)
	m := len(b)

	if len(H) != n || (m > 0 && len(A) != m) {
		return nil, nil, ErrWrongInput
	}
	for i := range n {
		if len(H[i]) != n {
			return nil, nil, ErrWrongInput
		}
	}
	if m > 0 {
		for i := range m {
			if len(A[i]) != n {
				return nil, nil, ErrWrongInput
			}
		}
	}

	if m == 0 {
		// No constraints: solve Hx = -c
		sol, err := GaussElim(H, scaleVec(-1, c))
		if err != nil {
			return nil, nil, err
		}
		return sol, nil, nil
	}

	// Build KKT matrix (n+m) × (n+m)
	size := n + m
	KKT := make(Matrix, size)
	for i := range size {
		KKT[i] = make([]float64, size)
	}

	for i := range n {
		for j := range n {
			KKT[i][j] = H[i][j]
		}
	}
	for i := range n {
		for j := range m {
			KKT[i][n+j] = A[j][i] // Aᵀ
		}
	}
	for i := range m {
		for j := range n {
			KKT[n+i][j] = A[i][j] // A
		}
	}

	rhs := make([]float64, size)
	for i := range n {
		rhs[i] = -c[i]
	}
	for i := range m {
		rhs[n+i] = b[i]
	}

	sol, err := GaussElim(KKT, rhs)
	if err != nil {
		return nil, nil, err
	}
	return sol[:n], sol[n:], nil
}

// QPActiveSet solves the inequality-constrained QP:
//
//	min  ½ xᵀHx + cᵀx
//	s.t. A_eq x  = b_eq
//	     A_ineq x >= b_ineq
//
// using the active-set method. Returns x, equality multipliers, inequality
// multipliers (indexed by original inequality constraint).
// A_eq / b_eq may be nil for unconstrained equality side.
// A_ineq / b_ineq may be nil for unconstrained inequality side.
func QPActiveSet(H Matrix, c []float64, A_eq Matrix, b_eq []float64, A_ineq Matrix, b_ineq []float64) ([]float64, []float64, []float64, error) {
	n := len(c)
	nEq := len(b_eq)
	nIneq := len(b_ineq)

	if len(H) != n {
		return nil, nil, nil, ErrWrongInput
	}

	// No inequalities: delegate directly
	if nIneq == 0 {
		x, lam, err := QPEqualityConstrained(H, c, A_eq, b_eq)
		lamEq := lam
		if lamEq == nil {
			lamEq = make([]float64, nEq)
		}
		return x, lamEq, make([]float64, 0), err
	}

	// Initial x: satisfy equality constraints (min-norm), or zero if none
	x := make([]float64, n)
	if nEq > 0 {
		var err error
		x, err = minNormFeasible(A_eq, b_eq)
		if err != nil {
			x = make([]float64, n)
		}
	}

	// Initialize active set: treat any violated inequality as active
	active := make([]bool, nIneq)
	for j := range nIneq {
		if dot(A_ineq[j], x) < b_ineq[j] {
			active[j] = true
		}
	}

	const maxIter = 300

	for range maxIter {
		// Collect active inequality row indices
		activeRows := make([]int, 0, nIneq)
		for j := range nIneq {
			if active[j] {
				activeRows = append(activeRows, j)
			}
		}

		nA := nEq + len(activeRows)
		var AComb Matrix
		var bComb []float64
		if nA > 0 {
			AComb = make(Matrix, nA)
			bComb = make([]float64, nA)
			for i := range nEq {
				AComb[i] = A_eq[i]
				bComb[i] = b_eq[i]
			}
			for i, row := range activeRows {
				AComb[nEq+i] = A_ineq[row]
				bComb[nEq+i] = b_ineq[row]
			}
		}

		xStar, lam, err := QPEqualityConstrained(H, c, AComb, bComb)
		if err != nil {
			// Degenerate active set — drop the last added inequality and continue
			if len(activeRows) > 0 {
				active[activeRows[len(activeRows)-1]] = false
				continue
			}
			return nil, nil, nil, err
		}

		d := subVec(xStar, x)
		dNorm := Norm(d)

		if dNorm < 1e-10 {
			// d ≈ 0: check KKT multipliers for active inequalities.
			// QPEqualityConstrained returns λ from Hx + Aᵀλ = -c, so the
			// proper inequality KKT multiplier is μ = -λ (must be ≥ 0).
			maxLam := 0.0
			maxIdx := -1
			for i, row := range activeRows {
				lam_i := lam[nEq+i]
				if lam_i > maxLam {
					maxLam = lam_i
					maxIdx = row
				}
			}
			if maxIdx == -1 {
				// All λ ≤ 0  →  μ = -λ ≥ 0  →  KKT satisfied
				x = xStar
				break
			}
			// Drop the constraint with the most positive λ (most negative μ)
			active[maxIdx] = false
		} else {
			// Find the maximum step α before any inactive constraint is violated
			alpha := 1.0
			addConstraint := -1
			for j := range nIneq {
				if active[j] {
					continue
				}
				// Constraint: A_ineq[j]·(x + α·d) >= b_ineq[j]
				// A_ineq[j]·d * α >= b_ineq[j] - A_ineq[j]·x
				ad := dot(A_ineq[j], d)
				slack := dot(A_ineq[j], x) - b_ineq[j]
				if ad < -1e-12 {
					alphaMax := slack / (-ad)
					if alphaMax < alpha {
						alpha = alphaMax
						addConstraint = j
					}
				}
			}

			for i := range n {
				x[i] += alpha * d[i]
			}
			if addConstraint >= 0 && alpha < 1.0-1e-12 {
				active[addConstraint] = true
			}
		}
	}

	// Recompute multipliers at the final x
	activeRows := make([]int, 0, nIneq)
	for j := range nIneq {
		if active[j] {
			activeRows = append(activeRows, j)
		}
	}

	lambdaEq := make([]float64, nEq)
	muIneq := make([]float64, nIneq)

	nA := nEq + len(activeRows)
	if nA > 0 {
		AComb := make(Matrix, nA)
		bComb := make([]float64, nA)
		for i := range nEq {
			AComb[i] = A_eq[i]
			bComb[i] = b_eq[i]
		}
		for i, row := range activeRows {
			AComb[nEq+i] = A_ineq[row]
			bComb[nEq+i] = b_ineq[row]
		}
		_, lam, err := QPEqualityConstrained(H, c, AComb, bComb)
		if err == nil {
			for i := range nEq {
				lambdaEq[i] = lam[i]
			}
			for i, row := range activeRows {
				muIneq[row] = -lam[nEq+i] // μ = -λ ≥ 0 at optimum
			}
		}
	}

	return x, lambdaEq, muIneq, nil
}
