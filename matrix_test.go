package numericalanalysis_test

import (
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestMatrix_Equal(t *testing.T) {
	t.Run("same size", func(t *testing.T) {
		matrix1 := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		matrix2 := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		if !matrix1.Equal(matrix2) {
			t.Errorf("matrix1 = %v, matrix2 = %v, want equal", matrix1, matrix2)
		}
	})

	t.Run("different size", func(t *testing.T) {
		matrix1 := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		matrix2 := numericalanalysis.Matrix{
			{1, 2},
			{4, 5},
		}

		if matrix1.Equal(matrix2) {
			t.Errorf("matrix1 = %v, matrix2 = %v, want not equal", matrix1, matrix2)
		}
	})

	t.Run("different values", func(t *testing.T) {
		matrix1 := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		matrix2 := numericalanalysis.Matrix{
			{7, 8, 9},
			{10, 11, 12},
		}

		if matrix1.Equal(matrix2) {
			t.Errorf("matrix1 = %v, matrix2 = %v, want not equal", matrix1, matrix2)
		}
	})
}

func TestMatrix_Add(t *testing.T) {
	t.Run("same size", func(t *testing.T) {
		matrix1 := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		matrix2 := numericalanalysis.Matrix{
			{7, 8, 9},
			{10, 11, 12},
		}

		result, err := matrix1.Add(matrix2)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		expected := numericalanalysis.Matrix{
			{8, 10, 12},
			{14, 16, 18},
		}
		if !result.Equal(expected) {
			t.Errorf("result = %v, want %v", result, expected)
		}
	})

	t.Run("different size", func(t *testing.T) {
		matrix1 := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		matrix2 := numericalanalysis.Matrix{
			{7, 8},
			{10, 11},
		}

		_, err := matrix1.Add(matrix2)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})
}

func TestMatrix_Transpose(t *testing.T) {
	t.Run("square", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9},
		}

		transposed := matrix.Transpose()
		expected := numericalanalysis.Matrix{
			{1, 4, 7},
			{2, 5, 8},
			{3, 6, 9},
		}
		if !transposed.Equal(expected) {
			t.Errorf("transposed = %v, want %v", transposed, expected)
		}
	})

	t.Run("non-square", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		transposed := matrix.Transpose()
		expected := numericalanalysis.Matrix{
			{1, 4},
			{2, 5},
			{3, 6},
		}
		if !transposed.Equal(expected) {
			t.Errorf("transposed = %v, want %v", transposed, expected)
		}
	})
}

func TestMatrix_Det(t *testing.T) {
	t.Run("zero", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9},
		}

		det, err := matrix.Det()
		if det != 0 {
			t.Errorf("det = %v, want 0", det)
		}
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
	})

	t.Run("nonzero", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 10},
		}

		det, err := matrix.Det()
		if det != -3 {
			t.Errorf("det = %v, want -3", det)
		}
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
	})

	t.Run("non-square", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		_, err := matrix.Det()
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}

		matrix = numericalanalysis.Matrix{}
		_, err = matrix.Det()
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})
}

func TestMatrix_Mul(t *testing.T) {
	t.Run("compatible matrices", func(t *testing.T) {
		matrix1 := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		matrix2 := numericalanalysis.Matrix{
			{7, 8},
			{9, 10},
			{11, 12},
		}

		result, err := matrix1.Mul(matrix2)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		expected := numericalanalysis.Matrix{
			{58, 64},
			{139, 154},
		}
		if !result.Equal(expected) {
			t.Errorf("result = %v, want %v", result, expected)
		}
	})

	t.Run("incompatible matrices", func(t *testing.T) {
		matrix1 := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		matrix2 := numericalanalysis.Matrix{
			{7, 8},
			{9, 10},
		}

		_, err := matrix1.Mul(matrix2)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("identity matrix", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9},
		}

		identity := numericalanalysis.Matrix{
			{1, 0, 0},
			{0, 1, 0},
			{0, 0, 1},
		}

		result, err := matrix.Mul(identity)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if !result.Equal(matrix) {
			t.Errorf("result = %v, want %v", result, matrix)
		}
	})

	t.Run("square matrices", func(t *testing.T) {
		matrix1 := numericalanalysis.Matrix{
			{1, 2},
			{3, 4},
		}

		matrix2 := numericalanalysis.Matrix{
			{5, 6},
			{7, 8},
		}

		result, err := matrix1.Mul(matrix2)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		expected := numericalanalysis.Matrix{
			{19, 22},
			{43, 50},
		}
		if !result.Equal(expected) {
			t.Errorf("result = %v, want %v", result, expected)
		}
	})
}

func TestMatrix_Inverse(t *testing.T) {
	t.Run("singular matrix", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2},
			{2, 4},
		}

		_, err := matrix.Inverse()
		if err != numericalanalysis.ErrSingularMatrix {
			t.Errorf("err = %v, want ErrSingularMatrix", err)
		}
	})

	t.Run("non-square matrix", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
		}

		_, err := matrix.Inverse()
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("identity matrix", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 0},
			{0, 1},
		}

		result, err := matrix.Inverse()
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if !result.Equal(matrix) {
			t.Errorf("result = %v, want %v", result, matrix)
		}
	})

	t.Run("square matrices", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 0},
			{0, -1, 2},
			{-1, 2, 0},
		}

		result, err := matrix.Inverse()
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		expected := numericalanalysis.Matrix{
			{0.5, 0, -0.5},
			{0.25, 0, 0.25},
			{0.125, 0.5, 0.125},
		}
		if !result.Equal(expected) {
			t.Errorf("result = %v, want %v", result, expected)
		}
	})
}
