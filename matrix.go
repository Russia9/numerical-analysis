package numericalanalysis

type Matrix [][]float64

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
