package numericalanalysis

func Cramer(matrix Matrix, free []float64) ([]float64, error) {
	// Find determinant of matrix
	det := matrix.Det()
	if det == 0 { // If determinant is zero, then there is no solution, need to use another method
		return nil, ErrNoSolution
	}

	// Find solution
	result := make([]float64, len(free))
	for i := 0; i < len(free); i++ {
		// Replace i-th column with free vector
		matrixCopy := make(Matrix, len(matrix))
		for j := 0; j < len(matrix); j++ {
			matrixCopy[j] = make([]float64, len(matrix[j]))
			copy(matrixCopy[j], matrix[j])
			matrixCopy[j][i] = free[j]
		}

		// Calculate determinant of new matrix and divide it by determinant of original matrix
		result[i] = matrixCopy.Det() / det
	}

	return result, nil
}
