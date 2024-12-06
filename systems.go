package numericalanalysis

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
	for i := 0; i < len(free); i++ {
		// Replace i-th column with free vector
		deltaIMatrix := make(Matrix, len(matrix))
		for j := 0; j < len(matrix); j++ {
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
