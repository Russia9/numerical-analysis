package numericalanalysis_test

import (
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

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
}
