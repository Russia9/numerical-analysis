package numericalanalysis_test

import (
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestMatrix_Det(t *testing.T) {
	t.Run("zero", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9},
		}

		det := matrix.Det()
		if det != 0 {
			t.Errorf("det = %v, want 0", det)
		}
	})

	t.Run("nonzero", func(t *testing.T) {
		matrix := numericalanalysis.Matrix{
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 10},
		}

		det := matrix.Det()
		if det != -3 {
			t.Errorf("det = %v, want -3", det)
		}
	})
}
