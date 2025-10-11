package numericalanalysis

import "math"

// Newton method for solving a system of nonlinear equations
// m: number of equations
// f[m]: system of nonlinear equations
// u0[m]: initial guess for the solution
// deltaU[m]: step size for each variable (for differential calculations)
// eps: tolerance for convergence
func Newton(f []func(u []float64) float64, u0 []float64, deltaU []float64, eps float64) ([]float64, error) {
	// Check input
	if len(f) == 0 || len(u0) == 0 || len(f) != len(u0) || eps <= 0 {
		return nil, ErrWrongInput
	}

	u := make([]float64, len(u0))
	copy(u, u0)

	for {
		// Calculate Jacobian matrix
		J := make(Matrix, len(f))
		for i := range f {
			J[i] = make([]float64, len(f))

			for j := range f {
				u_cur := make([]float64, len(u))
				copy(u_cur, u)

				u_cur[j] += deltaU[j]

				J[i][j] = (f[i](u_cur) - f[i](u)) / deltaU[j]
			}
		}

		// Calculate residual vector
		r := make([]float64, len(f))
		for i := range f {
			r[i] = f[i](u)
		}

		// Check convergence
		norm := 0.
		for i := range r {
			norm += r[i] * r[i]
		}
		norm = math.Sqrt(norm)
		if norm < eps {
			break
		}

		// Solve linear system J * du = -r
		du, err := Cramer(J, r)
		if err != nil {
			return nil, err
		}

		// Update solution
		for i := range u {
			u[i] -= du[i]
		}
	}

	return u, nil
}
