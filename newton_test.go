package numericalanalysis_test

import (
	"math"
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestNewton(t *testing.T) {
	t.Run("input validation - empty functions", func(t *testing.T) {
		f := []func(u []float64) float64{}
		u0 := []float64{1.0}
		deltaU := []float64{0.01}
		eps := 1e-6

		_, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("input validation - empty initial guess", func(t *testing.T) {
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] - 4 },
		}
		u0 := []float64{}
		deltaU := []float64{0.01}
		eps := 1e-6

		_, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("input validation - mismatched lengths", func(t *testing.T) {
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] - 4 },
			func(u []float64) float64 { return u[0] - 2 },
		}
		u0 := []float64{1.0}
		deltaU := []float64{0.01}
		eps := 1e-6

		_, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("input validation - non-positive tolerance", func(t *testing.T) {
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] - 4 },
		}
		u0 := []float64{1.0}
		deltaU := []float64{0.01}
		eps := 0.0

		_, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("single equation - quadratic", func(t *testing.T) {
		// f(x) = x^2 - 4 = 0, solutions: x = ±2
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] - 4 },
		}
		u0 := []float64{3.0} // start near positive root
		deltaU := []float64{0.01}
		eps := 1e-6

		result, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 1 {
			t.Errorf("len(result) = %v, want 1", len(result))
		}
		if math.Abs(result[0]-2.0) > 1e-5 {
			t.Errorf("result[0] = %v, want ~2.0", result[0])
		}
	})

	t.Run("single equation - negative root", func(t *testing.T) {
		// f(x) = x^2 - 4 = 0, solutions: x = ±2
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] - 4 },
		}
		u0 := []float64{-3.0} // start near negative root
		deltaU := []float64{0.01}
		eps := 1e-6

		result, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 1 {
			t.Errorf("len(result) = %v, want 1", len(result))
		}
		if math.Abs(result[0]+2.0) > 1e-5 {
			t.Errorf("result[0] = %v, want ~-2.0", result[0])
		}
	})

	t.Run("system of linear equations", func(t *testing.T) {
		// f1(x,y) = x + y - 3 = 0
		// f2(x,y) = x - y - 1 = 0
		// Solution: x = 2, y = 1
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0] + u[1] - 3 },
			func(u []float64) float64 { return u[0] - u[1] - 1 },
		}
		u0 := []float64{0.0, 0.0}
		deltaU := []float64{0.01, 0.01}
		eps := 1e-6

		result, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 2 {
			t.Errorf("len(result) = %v, want 2", len(result))
		}
		if math.Abs(result[0]-2.0) > 1e-5 {
			t.Errorf("result[0] = %v, want ~2.0", result[0])
		}
		if math.Abs(result[1]-1.0) > 1e-5 {
			t.Errorf("result[1] = %v, want ~1.0", result[1])
		}
	})

	t.Run("system of nonlinear equations", func(t *testing.T) {
		// f1(x,y) = x^2 + y^2 - 1 = 0  (unit circle)
		// f2(x,y) = x - y = 0           (line y = x)
		// Solution: x = y = ±1/√2
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] + u[1]*u[1] - 1 },
			func(u []float64) float64 { return u[0] - u[1] },
		}
		u0 := []float64{0.5, 0.5} // start near positive solution
		deltaU := []float64{0.01, 0.01}
		eps := 1e-6

		result, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 2 {
			t.Errorf("len(result) = %v, want 2", len(result))
		}
		expected := 1.0 / math.Sqrt(2)
		if math.Abs(result[0]-expected) > 1e-5 {
			t.Errorf("result[0] = %v, want ~%v", result[0], expected)
		}
		if math.Abs(result[1]-expected) > 1e-5 {
			t.Errorf("result[1] = %v, want ~%v", result[1], expected)
		}
	})

	t.Run("already at solution", func(t *testing.T) {
		// f(x) = x - 2 = 0, solution: x = 2
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0] - 2 },
		}
		u0 := []float64{2.0} // already at solution
		deltaU := []float64{0.01}
		eps := 1e-6

		result, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 1 {
			t.Errorf("len(result) = %v, want 1", len(result))
		}
		if math.Abs(result[0]-2.0) > 1e-5 {
			t.Errorf("result[0] = %v, want ~2.0", result[0])
		}
	})

	t.Run("different tolerance values", func(t *testing.T) {
		// f(x) = x^2 - 4 = 0
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] - 4 },
		}
		u0 := []float64{3.0}
		deltaU := []float64{0.01}

		// Test with looser tolerance
		eps1 := 1e-3
		result1, err := numericalanalysis.SENewton(f, u0, deltaU, eps1)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if math.Abs(result1[0]-2.0) > 1e-2 {
			t.Errorf("result1[0] = %v, want ~2.0", result1[0])
		}

		// Test with tighter tolerance
		eps2 := 1e-10
		result2, err := numericalanalysis.SENewton(f, u0, deltaU, eps2)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if math.Abs(result2[0]-2.0) > 1e-9 {
			t.Errorf("result2[0] = %v, want ~2.0", result2[0])
		}
	})

	t.Run("different step sizes", func(t *testing.T) {
		// f(x) = x^2 - 4 = 0
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] - 4 },
		}
		u0 := []float64{3.0}
		eps := 1e-6

		// Test with different step sizes
		testCases := []float64{0.001, 0.01, 0.1}
		for _, delta := range testCases {
			deltaU := []float64{delta}
			result, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
			if err != nil {
				t.Errorf("deltaU = %v: err = %v, want nil", delta, err)
				continue
			}
			if math.Abs(result[0]-2.0) > 1e-5 {
				t.Errorf("deltaU = %v: result[0] = %v, want ~2.0", delta, result[0])
			}
		}
	})

	t.Run("verify solution accuracy", func(t *testing.T) {
		// System: f1(x,y) = x^2 - y - 1 = 0, f2(x,y) = x + y^2 - 3 = 0
		f := []func(u []float64) float64{
			func(u []float64) float64 { return u[0]*u[0] - u[1] - 1 },
			func(u []float64) float64 { return u[0] + u[1]*u[1] - 3 },
		}
		u0 := []float64{1.0, 1.0}
		deltaU := []float64{0.01, 0.01}
		eps := 1e-8

		result, err := numericalanalysis.SENewton(f, u0, deltaU, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}

		// Verify that the result satisfies the original equations
		residual1 := f[0](result)
		residual2 := f[1](result)
		if math.Abs(residual1) > 1e-6 {
			t.Errorf("residual1 = %v, want ~0", residual1)
		}
		if math.Abs(residual2) > 1e-6 {
			t.Errorf("residual2 = %v, want ~0", residual2)
		}
	})
}

func TestDampedNewtonExtremum(t *testing.T) {
	t.Run("input validation - invalid C1 zero", func(t *testing.T) {
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{1.0}
		deltaX := []float64{0.01}
		alpha0 := 1.0
		C1 := 0.0
		eps := 1e-6

		_, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("input validation - invalid C1 negative", func(t *testing.T) {
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{1.0}
		deltaX := []float64{0.01}
		alpha0 := 100.0
		C1 := -0.5
		eps := 1e-6

		_, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("input validation - invalid C1 one", func(t *testing.T) {
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{1.0}
		deltaX := []float64{0.01}
		alpha0 := 100.0
		C1 := 1.0
		eps := 1e-6

		_, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("input validation - invalid C1 greater than one", func(t *testing.T) {
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{1.0}
		deltaX := []float64{0.01}
		alpha0 := 100.0
		C1 := 1.5
		eps := 1e-6

		_, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("input validation - non-positive tolerance", func(t *testing.T) {
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{1.0}
		deltaX := []float64{0.01}
		alpha0 := 1.0
		C1 := 0.5
		eps := 0.0

		_, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("single variable - quadratic minimum", func(t *testing.T) {
		// f(x) = [(x - 3)^2], minimum at x = 3
		f := func(x []float64) float64 {
			return (x[0] - 3) * (x[0] - 3)
		}
		x0 := []float64{1.0}
		deltaX := []float64{1e-4}
		alpha0 := 100.0
		C1 := 0.5
		eps := 1e-6

		result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 1 {
			t.Errorf("len(result) = %v, want 1", len(result))
		}
		if math.Abs(result[0]-3.0) > 1e-3 {
			t.Errorf("result[0] = %v, want ~3.0", result[0])
		}
	})

	t.Run("single variable - quadratic at origin", func(t *testing.T) {
		// f(x) = x^2, minimum at x = 0
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{5.0}
		deltaX := []float64{1e-4}
		alpha0 := 100.0
		C1 := 0.5
		eps := 1e-6

		result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 1 {
			t.Errorf("len(result) = %v, want 1", len(result))
		}
		if math.Abs(result[0]) > 1e-3 {
			t.Errorf("result[0] = %v, want ~0.0", result[0])
		}
	})

	t.Run("two variables - paraboloid minimum", func(t *testing.T) {
		// f(x,y) = x^2 + y^2, minimum at (0, 0)
		f := func(x []float64) float64 {
			return x[0]*x[0] + x[1]*x[1]
		}
		x0 := []float64{3.0, 4.0}
		deltaX := []float64{1e-4, 1e-4}
		alpha0 := 100.0
		C1 := 0.5
		eps := 1e-6

		result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 2 {
			t.Errorf("len(result) = %v, want 2", len(result))
		}
		if math.Abs(result[0]) > 1e-3 {
			t.Errorf("result[0] = %v, want ~0.0", result[0])
		}
		if math.Abs(result[1]) > 1e-3 {
			t.Errorf("result[1] = %v, want ~0.0", result[1])
		}
	})

	t.Run("two variables - shifted paraboloid", func(t *testing.T) {
		// f(x,y) = [(x-1)^2 + (y-2)^2, (x-1)^2 + (y-2)^2], minimum at (1, 2)
		f := func(x []float64) float64 {
			return (x[0]-1)*(x[0]-1) + (x[1]-2)*(x[1]-2)
		}
		x0 := []float64{5.0, 5.0}
		deltaX := []float64{1e-4, 1e-4}
		alpha0 := 100.0
		C1 := 0.5
		eps := 1e-6

		result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 2 {
			t.Errorf("len(result) = %v, want 2", len(result))
		}
		if math.Abs(result[0]-1.0) > 1e-3 {
			t.Errorf("result[0] = %v, want ~1.0", result[0])
		}
		if math.Abs(result[1]-2.0) > 1e-3 {
			t.Errorf("result[1] = %v, want ~2.0", result[1])
		}
	})

	t.Run("two variables - rosenbrock function", func(t *testing.T) {
		// f(x,y) = [(1-x)^2 + 100(y-x^2)^2, (1-x)^2 + 100(y-x^2)^2], minimum at (1, 1)
		f := func(x []float64) float64 {
			return (1-x[0])*(1-x[0]) + 100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])
		}
		x0 := []float64{0.0, 0.0}
		deltaX := []float64{1e-4, 1e-4}
		alpha0 := 100.0
		C1 := 0.5
		eps := 1e-4

		result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 2 {
			t.Errorf("len(result) = %v, want 2", len(result))
		}
		if math.Abs(result[0]-1.0) > 0.1 {
			t.Errorf("result[0] = %v, want ~1.0", result[0])
		}
		if math.Abs(result[1]-1.0) > 0.1 {
			t.Errorf("result[1] = %v, want ~1.0", result[1])
		}
	})

	t.Run("already at minimum", func(t *testing.T) {
		// f(x) = x^2, start at minimum
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{0.0}
		deltaX := []float64{1e-6}
		alpha0 := 100.0
		C1 := 0.5
		eps := 1e-4

		result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 1 {
			t.Errorf("len(result) = %v, want 1", len(result))
		}
		if math.Abs(result[0]) > 1e-5 {
			t.Errorf("result[0] = %v, want ~0.0", result[0])
		}
	})

	t.Run("different damping factors", func(t *testing.T) {
		// f(x) = x^2, minimum at x = 0
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{5.0}
		deltaX := []float64{1e-4}
		C1 := 0.5
		eps := 1e-6

		testCases := []float64{1, 10.0, 100.0}
		for _, alpha0 := range testCases {
			result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
			if err != nil {
				t.Errorf("alpha0 = %v: err = %v, want nil", alpha0, err)
				continue
			}
			if math.Abs(result[0]) > 1e-3 {
				t.Errorf("alpha0 = %v: result[0] = %v, want ~0.0", alpha0, result[0])
			}
		}
	})

	t.Run("different C1 values", func(t *testing.T) {
		// f(x) = x^2, minimum at x = 0
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{3.0}
		deltaX := []float64{1e-4}
		alpha0 := 1.0
		eps := 1e-6

		testCases := []float64{0.1, 0.3, 0.5, 0.7, 0.9}
		for _, C1 := range testCases {
			result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
			if err != nil {
				t.Errorf("C1 = %v: err = %v, want nil", C1, err)
				continue
			}
			if math.Abs(result[0]) > 1e-3 {
				t.Errorf("C1 = %v: result[0] = %v, want ~0.0", C1, result[0])
			}
		}
	})

	t.Run("different tolerance values", func(t *testing.T) {
		// f(x) = x^2, minimum at x = 0
		f := func(x []float64) float64 {
			return x[0] * x[0]
		}
		x0 := []float64{3.0}
		deltaX := []float64{1e-5}
		alpha0 := 100.0
		C1 := 0.5

		// Test with looser tolerance
		eps1 := 1e-3
		result1, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps1)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if math.Abs(result1[0]) > 1e-2 {
			t.Errorf("result1[0] = %v, want ~0.0", result1[0])
		}

		// Test with tighter tolerance
		eps2 := 1e-8
		result2, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps2)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if math.Abs(result2[0]) > 1e-5 {
			t.Errorf("result2[0] = %v, want ~0.0", result2[0])
		}
	})

	t.Run("three variables - sphere", func(t *testing.T) {
		// f(x,y,z) = [x^2 + y^2 + z^2, x^2 + y^2 + z^2, x^2 + y^2 + z^2], minimum at (0, 0, 0)
		f := func(x []float64) float64 {
			return x[0]*x[0] + x[1]*x[1] + x[2]*x[2]
		}
		x0 := []float64{2.0, 3.0, 4.0}
		deltaX := []float64{1e-3, 1e-3, 1e-3}
		alpha0 := 100.0
		C1 := 0.5
		eps := 1e-6

		result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}
		if len(result) != 3 {
			t.Errorf("len(result) = %v, want 3", len(result))
		}
		if math.Abs(result[0]) > 1e-3 {
			t.Errorf("result[0] = %v, want ~0.0", result[0])
		}
		if math.Abs(result[1]) > 1e-3 {
			t.Errorf("result[1] = %v, want ~0.0", result[1])
		}
		if math.Abs(result[2]) > 1e-3 {
			t.Errorf("result[2] = %v, want ~0.0", result[2])
		}
	})

	t.Run("verify minimum accuracy", func(t *testing.T) {
		// f(x,y) = [x^2 + y^2, x^2 + y^2], minimum at (0, 0)
		f := func(x []float64) float64 {
			return x[0]*x[0] + x[1]*x[1]
		}
		x0 := []float64{5.0, 5.0}
		deltaX := []float64{0.0001, 0.0001}
		alpha0 := 100.0
		C1 := 0.5
		eps := 1e-8

		result, err := numericalanalysis.DampedNewtonExtremum(f, x0, deltaX, alpha0, C1, eps)
		if err != nil {
			t.Errorf("err = %v, want nil", err)
		}

		// Verify function value is close to minimum
		fval := f(result)
		if math.Abs(fval) > 1e-4 {
			t.Errorf("f(result) = %v, want ~0.0", fval)
		}
	})
}
