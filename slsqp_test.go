package numericalanalysis_test

import (
	"math"
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestSLSQP(t *testing.T) {
	approx := func(a, b, tol float64) bool { return math.Abs(a-b) < tol }
	dx := []float64{1e-5, 1e-5}

	t.Run("input validation", func(t *testing.T) {
		f := func(x []float64) float64 { return x[0]*x[0] + x[1]*x[1] }

		_, err := numericalanalysis.SLSQP(f, []float64{1, 1}, []float64{-1e-5, 1e-5}, nil, nil, 1e-6, 100)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (negative deltaX)", err)
		}
		_, err = numericalanalysis.SLSQP(f, []float64{}, dx, nil, nil, 1e-6, 100)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (empty x0)", err)
		}
		_, err = numericalanalysis.SLSQP(f, []float64{1, 1}, dx, nil, nil, 0, 100)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (eps=0)", err)
		}
		_, err = numericalanalysis.SLSQP(f, []float64{1, 1}, dx, nil, nil, 1e-6, 0)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (maxIter=0)", err)
		}
	})

	t.Run("unconstrained — bowl", func(t *testing.T) {
		// min (x-1)² + (y-2)²  →  solution (1, 2)
		f := func(x []float64) float64 {
			return (x[0]-1)*(x[0]-1) + (x[1]-2)*(x[1]-2)
		}
		result, err := numericalanalysis.SLSQP(f, []float64{0, 0}, dx, nil, nil, 1e-6, 200)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(result[0], 1, 1e-4) || !approx(result[1], 2, 1e-4) {
			t.Errorf("result = %v, want [1 2]", result)
		}
	})

	t.Run("equality constraint — min x²+y² s.t. x+y=1", func(t *testing.T) {
		// Solution (0.5, 0.5)
		f := func(x []float64) float64 { return x[0]*x[0] + x[1]*x[1] }
		eq := func(x []float64) float64 { return x[0] + x[1] - 1 }
		result, err := numericalanalysis.SLSQP(f, []float64{0, 1}, dx, []func([]float64) float64{eq}, nil, 1e-6, 300)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(result[0], 0.5, 1e-4) || !approx(result[1], 0.5, 1e-4) {
			t.Errorf("result = %v, want [0.5 0.5]", result)
		}
	})

	t.Run("inequality constraint — inactive", func(t *testing.T) {
		// min (x-1)²+(y-2)² s.t. x+y >= 0  →  unconstrained min (1,2) is feasible
		f := func(x []float64) float64 { return (x[0]-1)*(x[0]-1) + (x[1]-2)*(x[1]-2) }
		ineq := func(x []float64) float64 { return x[0] + x[1] }
		result, err := numericalanalysis.SLSQP(f, []float64{0, 0}, dx, nil, []func([]float64) float64{ineq}, 1e-6, 200)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(result[0], 1, 1e-4) || !approx(result[1], 2, 1e-4) {
			t.Errorf("result = %v, want [1 2]", result)
		}
	})

	t.Run("inequality constraint — active", func(t *testing.T) {
		// min (x-3)²+(y-3)² s.t. x+y <= 4  (→  -(x+y) >= -4)
		// Unconstrained min (3,3) violates constraint (6 > 4); solution on x+y=4: (2,2)
		f := func(x []float64) float64 { return (x[0]-3)*(x[0]-3) + (x[1]-3)*(x[1]-3) }
		ineq := func(x []float64) float64 { return 4 - x[0] - x[1] }
		result, err := numericalanalysis.SLSQP(f, []float64{1, 1}, dx, nil, []func([]float64) float64{ineq}, 1e-6, 300)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(result[0], 2, 1e-3) || !approx(result[1], 2, 1e-3) {
			t.Errorf("result = %v, want [2 2]", result)
		}
	})

	t.Run("equality + inequality", func(t *testing.T) {
		// min x²+y² s.t. x+y=1 and x >= 0 and y >= 0
		// Solution: (0.5, 0.5) — both inequalities inactive
		f := func(x []float64) float64 { return x[0]*x[0] + x[1]*x[1] }
		eq := func(x []float64) float64 { return x[0] + x[1] - 1 }
		ineqX := func(x []float64) float64 { return x[0] }
		ineqY := func(x []float64) float64 { return x[1] }
		result, err := numericalanalysis.SLSQP(
			f, []float64{0.5, 0.5}, dx,
			[]func([]float64) float64{eq},
			[]func([]float64) float64{ineqX, ineqY},
			1e-6, 300,
		)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(result[0], 0.5, 1e-4) || !approx(result[1], 0.5, 1e-4) {
			t.Errorf("result = %v, want [0.5 0.5]", result)
		}
	})

	t.Run("ErrDidNotConverge - maxIter=1", func(t *testing.T) {
		f := func(x []float64) float64 { return (x[0]-5)*(x[0]-5) + (x[1]-5)*(x[1]-5) }
		_, err := numericalanalysis.SLSQP(f, []float64{0, 0}, dx, nil, nil, 1e-12, 1)
		if err != numericalanalysis.ErrDidNotConverge {
			t.Errorf("err = %v, want ErrDidNotConverge", err)
		}
	})

	t.Run("Rosenbrock unconstrained", func(t *testing.T) {
		// min (1-x)² + 100(y-x²)²  →  solution (1, 1)
		f := func(x []float64) float64 {
			return (1-x[0])*(1-x[0]) + 100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])
		}
		result, err := numericalanalysis.SLSQP(f, []float64{0, 0}, dx, nil, nil, 1e-5, 2000)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(result[0], 1, 1e-3) || !approx(result[1], 1, 1e-3) {
			t.Errorf("result = %v, want [1 1]", result)
		}
	})
}
