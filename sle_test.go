package numericalanalysis_test

import (
	"math"
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestGaussElim(t *testing.T) {
	approxEq := func(a, b, tol float64) bool { return math.Abs(a-b) < tol }

	t.Run("2x2 system", func(t *testing.T) {
		// 2x + y = 5, x + 3y = 10  →  x=1, y=3
		matrix := numericalanalysis.Matrix{
			{2, 1},
			{1, 3},
		}
		free := []float64{5, 10}
		result, err := numericalanalysis.GaussElim(matrix, free)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approxEq(result[0], 1, 1e-10) || !approxEq(result[1], 3, 1e-10) {
			t.Errorf("got %v, want [1 3]", result)
		}
	})

	t.Run("3x3 system", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, -1},
			{2, 1, 3},
			{-1, 3, 2},
		}
		free := []float64{3, 13, 4}
		result, err := numericalanalysis.GaussElim(matrix, free)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		// Verify by substituting back
		for i := range 3 {
			sum := 0.0
			for j := range 3 {
				sum += matrix[i][j] * result[j]
			}
			if !approxEq(sum, free[i], 1e-10) {
				t.Errorf("row %d: got sum %v, want %v", i, sum, free[i])
			}
		}
	})

	t.Run("singular matrix", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2},
			{2, 4},
		}
		free := []float64{1, 2}
		_, err := numericalanalysis.GaussElim(matrix, free)
		if err != numericalanalysis.ErrSingularMatrix {
			t.Errorf("err = %v, want ErrSingularMatrix", err)
		}
	})

	t.Run("wrong input", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2},
		}
		free := []float64{1, 2}
		_, err := numericalanalysis.GaussElim(matrix, free)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})
}

func TestCramer(t *testing.T) {
	t.Run("ok", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{2, 3},
			{1, 1},
		}
		free := []float64{1, -1}

		result, err := numericalanalysis.Cramer(matrix, free)
		if err != nil {
			t.Errorf("unexpected error: %v", err)
		}

		expected := []float64{-4, 3}
		for i := range result {
			if result[i] != expected[i] {
				t.Errorf("unexpected result[%d]: got %v, want %v", i, result[i], expected[i])
			}
		}
	})

	t.Run("dimension mismatch", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2},
			{3, 4},
		}
		free := []float64{1, 2, 3} // wrong length
		_, err := numericalanalysis.Cramer(matrix, free)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("singular matrix", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2},
			{2, 4},
		}
		free := []float64{1, 2}
		_, err := numericalanalysis.Cramer(matrix, free)
		if err != numericalanalysis.ErrNoSolution {
			t.Errorf("err = %v, want ErrNoSolution", err)
		}
	})
}
