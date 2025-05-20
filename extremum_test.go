package numericalanalysis_test

import (
	"math"
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestBisectionExtremum(t *testing.T) {
	tests := map[string]struct {
		f   func(float64) float64
		x   float64
		a   float64
		b   float64
		max bool
	}{
		"(x-2)^2": {
			f:   func(x float64) float64 { return (x - 2) * (x - 2) },
			x:   2,
			a:   0,
			b:   4,
			max: false,
		},
		"sin(x)": {
			f:   math.Sin,
			x:   math.Pi / 2,
			a:   0,
			b:   math.Pi,
			max: true,
		},
	}

	for name, test := range tests {
		x := numericalanalysis.BisectionExtremum(test.f, test.a, test.b, 1e-6, test.max)
		if math.Abs(x-test.x) > 1e-6 {
			t.Errorf("BisectionExtremum(%s) = %v, want %v", name, x, test.x)
		}
	}
}
