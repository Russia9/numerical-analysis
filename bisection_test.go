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
		t.Run(name, func(t *testing.T) {
			x, err := numericalanalysis.BisectionExtremum(test.f, test.a, test.b, 1e-6, test.max)
			if math.Abs(x-test.x) > 1e-6 {
				t.Errorf("BisectionExtremum(%s) = %v, want %v", name, x, test.x)
			}
			if err != nil {
				t.Errorf("BisectionExtremum(%s) returned error: %v", name, err)
			}
		})
	}
}

func TestBisectionValue(t *testing.T) {
	tests := map[string]struct {
		f        func(float64) float64
		a        float64
		b        float64
		value    float64
		expected float64
		backward bool
	}{
		"linear function": {
			f:        func(x float64) float64 { return 2*x + 1 },
			a:        0,
			b:        5,
			value:    7,
			expected: 3,
			backward: false,
		},
		"linear function backward": {
			f:        func(x float64) float64 { return -2*x + 1 },
			a:        -10,
			b:        10,
			value:    3,
			expected: -1,
			backward: true,
		},
		"quadratic function": {
			f:        func(x float64) float64 { return x * x },
			a:        1,
			b:        4,
			value:    9,
			expected: 3,
			backward: false,
		},
		"exponential function": {
			f:        func(x float64) float64 { return math.Exp(x) },
			a:        0,
			b:        2,
			value:    math.E,
			expected: 1,
			backward: false,
		},
		"exponential function backward": {
			f:        func(x float64) float64 { return -1 * math.Exp(2*x) },
			a:        0,
			b:        2,
			value:    -10,
			expected: 1.1512925465,
			backward: true,
		},
	}

	for name, test := range tests {
		t.Run(name, func(t *testing.T) {
			x, err := numericalanalysis.BisectionValue(test.f, test.a, test.b, test.value, 1e-6, test.backward)
			if math.Abs(x-test.expected) > 1e-5 {
				t.Errorf("BisectionValue(%s) = %v, want %v", name, x, test.expected)
			}
			if err != nil {
				t.Errorf("BisectionValue(%s) returned error: %v", name, err)
			}
		})
	}
}
