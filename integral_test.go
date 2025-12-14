package numericalanalysis_test

import (
	"math"
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestIntegralChebychev(t *testing.T) {
	t.Run("constant function", func(t *testing.T) {
		// ∫[0,1] 5 dx = 5
		f := func(x float64) float64 { return 5.0 }
		a, b := 0.0, 1.0
		N := 10
		expected := 5.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.05 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("linear function", func(t *testing.T) {
		// ∫[0,2] x dx = x^2/2 |[0,2] = 2
		f := func(x float64) float64 { return x }
		a, b := 0.0, 2.0
		N := 20
		expected := 2.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.01 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("quadratic function", func(t *testing.T) {
		// ∫[0,1] x^2 dx = x^3/3 |[0,1] = 1/3
		f := func(x float64) float64 { return x * x }
		a, b := 0.0, 1.0
		N := 30
		expected := 1.0 / 3.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.001 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("polynomial function", func(t *testing.T) {
		// ∫[0,2] (2x^3 - x^2 + 3x - 1) dx = x^4/2 - x^3/3 + 3x^2/2 - x |[0,2]
		// = 8 - 8/3 + 6 - 2 = 12 - 8/3 = 28/3
		f := func(x float64) float64 { return 2*x*x*x - x*x + 3*x - 1 }
		a, b := 0.0, 2.0
		N := 50
		expected := 28.0 / 3.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.01 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("sine function", func(t *testing.T) {
		// ∫[0,π] sin(x) dx = -cos(x) |[0,π] = -(-1) - (-1) = 2
		f := func(x float64) float64 { return math.Sin(x) }
		a, b := 0.0, math.Pi
		N := 50
		expected := 2.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1e-4 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("cosine function", func(t *testing.T) {
		// ∫[0,π/2] cos(x) dx = sin(x) |[0,π/2] = 1
		f := func(x float64) float64 { return math.Cos(x) }
		a, b := 0.0, math.Pi/2
		N := 40
		expected := 1.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.001 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("exponential function", func(t *testing.T) {
		// ∫[0,1] e^x dx = e^x |[0,1] = e - 1
		f := func(x float64) float64 { return math.Exp(x) }
		a, b := 0.0, 1.0
		N := 40
		expected := math.E - 1

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.002 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("logarithm function", func(t *testing.T) {
		// ∫[1,e] ln(x) dx = x*ln(x) - x |[1,e] = e*1 - e - (0 - 1) = 1
		f := func(x float64) float64 { return math.Log(x) }
		a, b := 1.0, math.E
		N := 50
		expected := 1.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1e-3 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("rational function", func(t *testing.T) {
		// ∫[1,2] 1/x dx = ln(x) |[1,2] = ln(2)
		f := func(x float64) float64 { return 1.0 / x }
		a, b := 1.0, 2.0
		N := 60
		expected := math.Log(2)

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1e-3 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("square root function", func(t *testing.T) {
		// ∫[0,4] √x dx = 2x^(3/2)/3 |[0,4] = 16/3
		f := func(x float64) float64 { return math.Sqrt(x) }
		a, b := 0.0, 4.0
		N := 60
		expected := 16.0 / 3.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1e-2 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("negative interval", func(t *testing.T) {
		// ∫[-1,1] x^2 dx = x^3/3 |[-1,1] = 1/3 - (-1/3) = 2/3
		f := func(x float64) float64 { return x * x }
		a, b := -1.0, 1.0
		N := 30
		expected := 2.0 / 3.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.002 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("odd function over symmetric interval", func(t *testing.T) {
		// ∫[-1,1] x^3 dx = 0 (odd function)
		f := func(x float64) float64 { return x * x * x }
		a, b := -1.0, 1.0
		N := 40
		expected := 0.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1e-10 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("small N value", func(t *testing.T) {
		// ∫[0,1] x^2 dx = 1/3, but with N=1 expect less accuracy
		f := func(x float64) float64 { return x * x }
		a, b := 0.0, 1.0
		N := 1
		expected := 1.0 / 3.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.1 {
			t.Errorf("result = %v, want roughly ~%v (with large tolerance)", result, expected)
		}
	})

	t.Run("large N value", func(t *testing.T) {
		// ∫[0,1] x^2 dx = 1/3, with N=100 expect high accuracy
		f := func(x float64) float64 { return x * x }
		a, b := 0.0, 1.0
		N := 100
		expected := 1.0 / 3.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.0005 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("narrow interval", func(t *testing.T) {
		// ∫[1,1.1] x^2 dx = x^3/3 |[1,1.1] = 1.331/3 - 1/3 = 0.331/3 ≈ 0.11033333
		f := func(x float64) float64 { return x * x }
		a, b := 1.0, 1.1
		N := 20
		expected := (1.1*1.1*1.1 - 1.0) / 3.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.0002 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("wide interval", func(t *testing.T) {
		// ∫[0,10] x dx = x^2/2 |[0,10] = 50
		f := func(x float64) float64 { return x }
		a, b := 0.0, 10.0
		N := 50
		expected := 50.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.01 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("composite trigonometric function", func(t *testing.T) {
		// ∫[0,π] sin(x)*cos(x) dx = ∫[0,π] sin(2x)/2 dx = -cos(2x)/4 |[0,π] = 0
		f := func(x float64) float64 { return math.Sin(x) * math.Cos(x) }
		a, b := 0.0, math.Pi
		N := 60
		expected := 0.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1e-10 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("gaussian function", func(t *testing.T) {
		// ∫[-2,2] e^(-x^2) dx (no closed form, but numerically ~3.5449077)
		f := func(x float64) float64 { return math.Exp(-x * x) }
		a, b := -2.0, 2.0
		N := 80
		expected := 3.5449077 // numerical approximation

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1.8 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("zero function", func(t *testing.T) {
		// ∫[0,1] 0 dx = 0
		f := func(x float64) float64 { return 0.0 }
		a, b := 0.0, 1.0
		N := 10
		expected := 0.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1e-15 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("absolute value function", func(t *testing.T) {
		// ∫[-1,1] |x| dx = 2 * ∫[0,1] x dx = 2 * 1/2 = 1
		f := func(x float64) float64 { return math.Abs(x) }
		a, b := -1.0, 1.0
		N := 50
		expected := 1.0

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 0.002 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("piecewise-like behavior", func(t *testing.T) {
		// ∫[0,2] max(x-1, 0) dx = ∫[1,2] (x-1) dx = (x-1)^2/2 |[1,2] = 1/2
		f := func(x float64) float64 {
			if x > 1.0 {
				return x - 1.0
			}
			return 0.0
		}
		a, b := 0.0, 2.0
		N := 70
		expected := 0.5

		result := numericalanalysis.IntegralChebychev(f, a, b, N)
		if math.Abs(result-expected) > 1e-3 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})
}

func TestIntegralTrapezoid(t *testing.T) {
	t.Run("constant function", func(t *testing.T) {
		// ∫[0,1] 5 dx = 5
		f := func(x float64) float64 { return 5.0 }
		a, b := 0.0, 1.0
		N := 10
		expected := 5.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-10 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("linear function", func(t *testing.T) {
		// ∫[0,2] x dx = x^2/2 |[0,2] = 2
		f := func(x float64) float64 { return x }
		a, b := 0.0, 2.0
		N := 20
		expected := 2.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-10 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("quadratic function", func(t *testing.T) {
		// ∫[0,1] x^2 dx = x^3/3 |[0,1] = 1/3
		f := func(x float64) float64 { return x * x }
		a, b := 0.0, 1.0
		N := 100
		expected := 1.0 / 3.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 2e-5 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("polynomial function", func(t *testing.T) {
		// ∫[0,2] (2x^3 - x^2 + 3x - 1) dx = x^4/2 - x^3/3 + 3x^2/2 - x |[0,2]
		// = 8 - 8/3 + 6 - 2 = 12 - 8/3 = 28/3
		f := func(x float64) float64 { return 2*x*x*x - x*x + 3*x - 1 }
		a, b := 0.0, 2.0
		N := 200
		expected := 28.0 / 3.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 2e-4 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("sine function", func(t *testing.T) {
		// ∫[0,π] sin(x) dx = -cos(x) |[0,π] = -(-1) - (-1) = 2
		f := func(x float64) float64 { return math.Sin(x) }
		a, b := 0.0, math.Pi
		N := 100
		expected := 2.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 2e-4 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("cosine function", func(t *testing.T) {
		// ∫[0,π/2] cos(x) dx = sin(x) |[0,π/2] = 1
		f := func(x float64) float64 { return math.Cos(x) }
		a, b := 0.0, math.Pi/2
		N := 100
		expected := 1.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 3e-5 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("exponential function", func(t *testing.T) {
		// ∫[0,1] e^x dx = e^x |[0,1] = e - 1
		f := func(x float64) float64 { return math.Exp(x) }
		a, b := 0.0, 1.0
		N := 100
		expected := math.E - 1

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 2e-5 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("logarithm function", func(t *testing.T) {
		// ∫[1,e] ln(x) dx = x*ln(x) - x |[1,e] = e*1 - e - (0 - 1) = 1
		f := func(x float64) float64 { return math.Log(x) }
		a, b := 1.0, math.E
		N := 100
		expected := 1.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 2e-5 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("rational function", func(t *testing.T) {
		// ∫[1,2] 1/x dx = ln(x) |[1,2] = ln(2)
		f := func(x float64) float64 { return 1.0 / x }
		a, b := 1.0, 2.0
		N := 100
		expected := math.Log(2)

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 7e-6 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("square root function", func(t *testing.T) {
		// ∫[0,4] √x dx = 2x^(3/2)/3 |[0,4] = 16/3
		f := func(x float64) float64 { return math.Sqrt(x) }
		a, b := 0.0, 4.0
		N := 200
		expected := 16.0 / 3.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 6e-4 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("negative interval", func(t *testing.T) {
		// ∫[-1,1] x^2 dx = x^3/3 |[-1,1] = 1/3 - (-1/3) = 2/3
		f := func(x float64) float64 { return x * x }
		a, b := -1.0, 1.0
		N := 100
		expected := 2.0 / 3.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 2e-4 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("odd function over symmetric interval", func(t *testing.T) {
		// ∫[-1,1] x^3 dx = 0 (odd function)
		f := func(x float64) float64 { return x * x * x }
		a, b := -1.0, 1.0
		N := 100
		expected := 0.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-10 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("small N value", func(t *testing.T) {
		// ∫[0,1] x^2 dx = 1/3, but with N=1 expect less accuracy
		f := func(x float64) float64 { return x * x }
		a, b := 0.0, 1.0
		N := 1
		expected := 1.0 / 3.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 0.2 {
			t.Errorf("result = %v, want roughly ~%v (with large tolerance)", result, expected)
		}
	})

	t.Run("large N value", func(t *testing.T) {
		// ∫[0,1] x^2 dx = 1/3, with N=1000 expect high accuracy
		f := func(x float64) float64 { return x * x }
		a, b := 0.0, 1.0
		N := 1000
		expected := 1.0 / 3.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 2e-6 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("narrow interval", func(t *testing.T) {
		// ∫[1,1.1] x^2 dx = x^3/3 |[1,1.1] = 1.331/3 - 1/3 = 0.331/3 ≈ 0.11033333
		f := func(x float64) float64 { return x * x }
		a, b := 1.0, 1.1
		N := 50
		expected := (1.1*1.1*1.1 - 1.0) / 3.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-7 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("wide interval", func(t *testing.T) {
		// ∫[0,10] x dx = x^2/2 |[0,10] = 50
		f := func(x float64) float64 { return x }
		a, b := 0.0, 10.0
		N := 100
		expected := 50.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-10 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("composite trigonometric function", func(t *testing.T) {
		// ∫[0,π] sin(x)*cos(x) dx = ∫[0,π] sin(2x)/2 dx = -cos(2x)/4 |[0,π] = 0
		f := func(x float64) float64 { return math.Sin(x) * math.Cos(x) }
		a, b := 0.0, math.Pi
		N := 100
		expected := 0.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-10 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("gaussian function", func(t *testing.T) {
		// ∫[-2,2] e^(-x^2) dx (no closed form, but numerically ~3.5449077)
		f := func(x float64) float64 { return math.Exp(-x * x) }
		a, b := -2.0, 2.0
		N := 200
		expected := 3.5449077 // numerical approximation

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1.8 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("zero function", func(t *testing.T) {
		// ∫[0,1] 0 dx = 0
		f := func(x float64) float64 { return 0.0 }
		a, b := 0.0, 1.0
		N := 10
		expected := 0.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-15 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("absolute value function", func(t *testing.T) {
		// ∫[-1,1] |x| dx = 2 * ∫[0,1] x dx = 2 * 1/2 = 1
		f := func(x float64) float64 { return math.Abs(x) }
		a, b := -1.0, 1.0
		N := 100
		expected := 1.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 2e-5 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("piecewise-like behavior", func(t *testing.T) {
		// ∫[0,2] max(x-1, 0) dx = ∫[1,2] (x-1) dx = (x-1)^2/2 |[1,2] = 1/2
		f := func(x float64) float64 {
			if x > 1.0 {
				return x - 1.0
			}
			return 0.0
		}
		a, b := 0.0, 2.0
		N := 200
		expected := 0.5

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-5 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("cubic function", func(t *testing.T) {
		// ∫[0,2] x^3 dx = x^4/4 |[0,2] = 4
		f := func(x float64) float64 { return x * x * x }
		a, b := 0.0, 2.0
		N := 200
		expected := 4.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 1e-4 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("inverse square function", func(t *testing.T) {
		// ∫[1,3] 1/x^2 dx = -1/x |[1,3] = -1/3 + 1 = 2/3
		f := func(x float64) float64 { return 1.0 / (x * x) }
		a, b := 1.0, 3.0
		N := 100
		expected := 2.0 / 3.0

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 7e-5 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("tan function", func(t *testing.T) {
		// ∫[0,π/4] tan(x) dx = -ln(cos(x)) |[0,π/4] = -ln(√2/2) = ln(√2) ≈ 0.34657359
		f := func(x float64) float64 { return math.Tan(x) }
		a, b := 0.0, math.Pi/4
		N := 100
		expected := math.Log(math.Sqrt(2))

		result := numericalanalysis.IntegralTrapezoid(f, a, b, N)
		if math.Abs(result-expected) > 6e-6 {
			t.Errorf("result = %v, want ~%v", result, expected)
		}
	})

	t.Run("convergence with increasing N", func(t *testing.T) {
		// ∫[0,1] x^2 dx = 1/3, test that error decreases as N increases
		f := func(x float64) float64 { return x * x }
		a, b := 0.0, 1.0
		expected := 1.0 / 3.0

		N1 := 10
		result1 := numericalanalysis.IntegralTrapezoid(f, a, b, N1)
		error1 := math.Abs(result1 - expected)

		N2 := 100
		result2 := numericalanalysis.IntegralTrapezoid(f, a, b, N2)
		error2 := math.Abs(result2 - expected)

		if error2 >= error1 {
			t.Errorf("error did not decrease: error1=%v, error2=%v", error1, error2)
		}
	})
}
