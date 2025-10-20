package numericalanalysis

import (
	"math"
)

// DampedNewtonExtremum method for finding an extremum of a function of many variables
// f: function to minimize
// x0: initial guess for the solution
// deltaX: step size for each variable (for differential calculations)
// alpha0: initial damping factor
// C1: coefficient, C1 in (0,1); C2=1/C1
// eps: tolerance
func DampedNewtonExtremum(f func(x []float64) []float64, x0 []float64, deltaX []float64, alpha0, C1 float64, eps float64) ([]float64, error) {
	// Check input
	if C1 <= 0 || C1 >= 1 || eps <= 0 {
		return nil, ErrWrongInput
	}

	C2 := 1 / C1
	x := x0
	alpha := alpha0

	for {
		y := f(x)

		// Calculate gradient
		grad := make([]float64, len(x))
		for i := range x {
			x_cur := make([]float64, len(x))
			copy(x_cur, x)

			x_cur[i] += deltaX[i]

			grad[i] = (f(x_cur)[i] - y[i]) / deltaX[i]
		}

		// Check grad G
		if Norm(grad) < eps {
			break
		}

		// Calculate Hessian matrix
		H := make(Matrix, len(x))
		for i := range x {
			H[i] = make([]float64, len(x))

			for j := range x {
				// A = g(..., x_i + Delta x_i, ..., x_j + Delta x_j, ...)
				// B = g(..., x_i + Delta x_i, ..., x_j
				// C = g(..., x_i, ..., x_j + Delta x_j, ...)
				// D = g(..., x_i, ..., x_j, ...)
				// H[i,j] = (A-B-C+D)/(Delta x_i * Delta x_j)
				x_A := make([]float64, len(x))
				x_B := make([]float64, len(x))
				x_C := make([]float64, len(x))
				x_D := make([]float64, len(x))

				copy(x_A, x)
				copy(x_B, x)
				copy(x_C, x)
				copy(x_D, x)

				x_A[i] += deltaX[i]
				x_A[j] += deltaX[j]
				x_B[i] += deltaX[i]
				x_C[j] += deltaX[j]

				A := f(x_A)[i]
				B := f(x_B)[i]
				C := f(x_C)[i]
				D := f(x_D)[i]

				H[i][j] = (A - B - C + D) / (deltaX[i] * deltaX[j])
			}
		}

		// X(i+1) = X(i) - (H(i)+alpha(i)I)^(-1)*grad(X(i))
		H1, err := H.Add(IdentityMatrix(len(x)).MulNumber(alpha))
		if err != nil {
			return nil, err
		}
		H2, err := H1.Inverse()
		if err != nil {
			return nil, err
		}
		H3, err := H2.Mul(Column(grad))
		if err != nil {
			return nil, err
		}
		H3 = H3.MulNumber(-1)
		H4, err := Column(x).Add(H3)
		if err != nil {
			return nil, err
		}

		// Convert back to float slice
		x1 := make([]float64, len(x))
		for i := range x {
			x1[i] = H4[i][0]
		}

		if Norm(f(x1)) < Norm(f(x)) {
			alpha = C1 * alpha
		} else {
			alpha = C2 * alpha
		}

		x = x1
	}

	return x, nil
}

// SENewton method for solving a system of nonlinear equations
// m: number of equations
// f[m]: system of nonlinear equations
// u0[m]: initial guess for the solution
// deltaU[m]: step size for each variable (for differential calculations)
// eps: tolerance for convergence
func SENewton(f []func(u []float64) float64, u0 []float64, deltaU []float64, eps float64) ([]float64, error) {
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
