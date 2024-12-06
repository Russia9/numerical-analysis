package numericalanalysis

type Matrix [][]float64

func (m Matrix) Det() float64 {
	// Check if the matrix is square
	if len(m) != len(m[0]) {
		panic("matrix is not square")
	}

	// Base cases
	n := len(m)
	if n == 1 {
		return m[0][0]
	}
	if n == 2 {
		return m[0][0]*m[1][1] - m[0][1]*m[1][0]
	}

	// Recursive case
	var det float64
	for i := 0; i < n; i++ {
		// Create submatrix b by removing the first row and i-th column
		submatrix := make(Matrix, n-1)
		for j := 0; j < n-1; j++ {
			submatrix[j] = make([]float64, n-1)
			for k := 0; k < n-1; k++ {
				if k < i {
					submatrix[j][k] = m[j+1][k]
				} else {
					submatrix[j][k] = m[j+1][k+1]
				}
			}
		}

		// Recursively calculate the determinant
		if i%2 == 0 {
			det += m[0][i] * submatrix.Det()
		} else {
			det -= m[0][i] * submatrix.Det()
		}
	}

	return det
}

func (m Matrix) Add(n Matrix) Matrix {
	if len(m) != len(n) || len(m[0]) != len(n[0]) {
		panic("matrix size mismatch")
	}
	result := make(Matrix, len(m))
	for i := range m {
		result[i] = make([]float64, len(m[i]))
		for j := range m[i] {
			result[i][j] = m[i][j] + n[i][j]
		}
	}
	return result
}
