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
