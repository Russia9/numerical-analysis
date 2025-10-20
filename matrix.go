package numericalanalysis

import "math"

// matrix.go
// Matrix operations

func Norm(x []float64) float64 {
	var sum float64
	for _, v := range x {
		sum += v * v
	}
	return math.Sqrt(sum)
}

type Matrix [][]float64

func IdentityMatrix(n int) Matrix {
	result := make(Matrix, n)
	for i := range result {
		result[i] = make([]float64, n)
		result[i][i] = 1
	}
	return result
}

func Column(x []float64) Matrix {
	result := make(Matrix, len(x))
	for i := range result {
		result[i] = make([]float64, 1)
		result[i][0] = x[i]
	}
	return result
}

func (m Matrix) Equal(n Matrix) bool {
	if len(m) != len(n) || len(m[0]) != len(n[0]) {
		return false
	}
	for i := range m {
		for j := range m[i] {
			if m[i][j] != n[i][j] {
				return false
			}
		}
	}
	return true
}

func (m Matrix) Add(n Matrix) (Matrix, error) {
	if len(m) != len(n) || len(m[0]) != len(n[0]) {
		return nil, ErrWrongInput
	}
	result := make(Matrix, len(m))
	for i := range m {
		result[i] = make([]float64, len(m[i]))
		for j := range m[i] {
			result[i][j] = m[i][j] + n[i][j]
		}
	}
	return result, nil
}

func (m Matrix) Transpose() Matrix {
	result := make(Matrix, len(m[0]))
	for i := range m[0] {
		result[i] = make([]float64, len(m))
		for j := range m {
			result[i][j] = m[j][i]
		}
	}
	return result
}

func (m Matrix) Det() (float64, error) {
	// Check if the matrix is square
	if len(m) == 0 || len(m) != len(m[0]) {
		return 0, ErrWrongInput
	}

	// Base cases
	n := len(m)
	if n == 1 {
		return m[0][0], nil
	}
	if n == 2 {
		return m[0][0]*m[1][1] - m[0][1]*m[1][0], nil
	}

	// Recursive case
	var det float64
	for i := range n {
		// Create submatrix b by removing the first row and i-th column
		submatrix := make(Matrix, n-1)
		for j := range n - 1 {
			submatrix[j] = make([]float64, n-1)
			for k := range n - 1 {
				if k < i {
					submatrix[j][k] = m[j+1][k]
				} else {
					submatrix[j][k] = m[j+1][k+1]
				}
			}
		}

		// Recursively calculate the determinant
		submatrixDet, err := submatrix.Det()
		if err != nil {
			return 0, err
		}

		if i%2 == 0 {
			det += m[0][i] * submatrixDet
		} else {
			det -= m[0][i] * submatrixDet
		}
	}

	return det, nil
}

func (m Matrix) Mul(n Matrix) (Matrix, error) {
	if len(m[0]) != len(n) {
		return nil, ErrWrongInput
	}
	result := make(Matrix, len(m))
	for i := range m {
		result[i] = make([]float64, len(n[0]))
		for j := range n[0] {
			for k := range n {
				result[i][j] += m[i][k] * n[k][j]
			}
		}
	}
	return result, nil
}

func (m Matrix) MulNumber(a float64) Matrix {
	result := make(Matrix, len(m))
	for i := range m {
		result[i] = make([]float64, len(m[0]))
		for j := range m[0] {
			result[i][j] = m[i][j] * a
		}
	}
	return result
}

func (m Matrix) Inverse() (Matrix, error) {
	det, err := m.Det()
	if err != nil {
		return nil, err
	}
	if det == 0 {
		return nil, ErrSingularMatrix
	}
	if len(m) != len(m[0]) {
		return nil, ErrWrongInput
	}
	if len(m) == 1 {
		return Matrix{{1 / m[0][0]}}, nil
	}

	n := len(m)
	result := make(Matrix, n)
	for i := range n {
		result[i] = make([]float64, n)
	}

	for i := range n {
		for j := range n {
			// Create submatrix by removing the i-th row and j-th column
			submatrix := make(Matrix, n-1)
			for k := range n - 1 {
				submatrix[k] = make([]float64, n-1)
				for l := range n - 1 {
					srcRow := k
					if k >= i {
						srcRow = k + 1
					}
					srcCol := l
					if l >= j {
						srcCol = l + 1
					}
					submatrix[k][l] = m[srcRow][srcCol]
				}
			}

			// Recursively calculate the determinant
			submatrixDet, err := submatrix.Det()
			if err != nil {
				return nil, err
			}

			if (i+j)%2 == 0 {
				result[j][i] = submatrixDet / det
			} else {
				result[j][i] = -submatrixDet / det
			}
		}
	}

	return result, nil
}
