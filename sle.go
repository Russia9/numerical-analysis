package numericalanalysis

import "math"

// sle.go
// Systems of linear equations solvers

// GaussElim solves the linear system Ax = b using Gaussian elimination with
// partial pivoting. O(n^3) — prefer this over Cramer for n > ~4.
func GaussElim(matrix Matrix, free []float64) ([]float64, error) {
	n := len(free)
	if len(matrix) != n {
		return nil, ErrWrongInput
	}
	for i := range n {
		if len(matrix[i]) != n {
			return nil, ErrWrongInput
		}
	}

	// Build augmented matrix [A | b]
	aug := make([][]float64, n)
	for i := range n {
		aug[i] = make([]float64, n+1)
		copy(aug[i], matrix[i])
		aug[i][n] = free[i]
	}

	// Forward elimination with partial pivoting
	for col := range n {
		// Find pivot row
		maxVal := math.Abs(aug[col][col])
		maxRow := col
		for row := col + 1; row < n; row++ {
			if v := math.Abs(aug[row][col]); v > maxVal {
				maxVal = v
				maxRow = row
			}
		}
		if maxVal < 1e-14 {
			return nil, ErrSingularMatrix
		}
		aug[col], aug[maxRow] = aug[maxRow], aug[col]

		for row := col + 1; row < n; row++ {
			factor := aug[row][col] / aug[col][col]
			for k := col; k <= n; k++ {
				aug[row][k] -= factor * aug[col][k]
			}
		}
	}

	// Back substitution
	result := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		result[i] = aug[i][n]
		for j := i + 1; j < n; j++ {
			result[i] -= aug[i][j] * result[j]
		}
		result[i] /= aug[i][i]
	}
	return result, nil
}

// minNormFeasible returns the minimum-norm x satisfying Ax = b,
// computed as x = Aᵀ (AAᵀ)⁻¹ b.
func minNormFeasible(A Matrix, b []float64) ([]float64, error) {
	m := len(b)
	n := len(A[0])
	// AAᵀ  (m×m)
	AAT := make(Matrix, m)
	for i := range m {
		AAT[i] = make([]float64, m)
		for j := range m {
			for k := range n {
				AAT[i][j] += A[i][k] * A[j][k]
			}
		}
	}
	z, err := GaussElim(AAT, b)
	if err != nil {
		return nil, err
	}
	// x = Aᵀ z
	x := make([]float64, n)
	for j := range n {
		for i := range m {
			x[j] += A[i][j] * z[i]
		}
	}
	return x, nil
}

func Cramer(matrix Matrix, free []float64) ([]float64, error) {
	// Check input
	if len(matrix) != len(free) {
		return nil, ErrWrongInput
	}

	// Find determinant of matrix
	delta, err := matrix.Det()
	if err != nil {
		return nil, err
	}
	if delta == 0 { // If determinant is zero, then there is no solution, need to use another method
		return nil, ErrNoSolution
	}

	// Find solution
	result := make([]float64, len(free))
	for i := range free {
		// Replace i-th column with free vector
		deltaIMatrix := make(Matrix, len(matrix))
		for j := range matrix {
			deltaIMatrix[j] = make([]float64, len(matrix[j]))
			copy(deltaIMatrix[j], matrix[j])
			deltaIMatrix[j][i] = free[j]
		}

		// Calculate determinant of new matrix and divide it by determinant of original matrix
		deltaI, err := deltaIMatrix.Det()
		if err != nil {
			return nil, err
		}
		result[i] = deltaI / delta
	}

	return result, nil
}
