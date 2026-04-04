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

	t.Run("rocket staging — 3-stage payload maximisation", func(t *testing.T) {
		// Variant-22 constants (from real application).
		// Minimise f = ∏ (1-αᵢ)/(μᵢ-αᵢ)  s.t.  V_ch + ΣWᵢ·ln(μᵢ) = 0
		// and αᵢ < μᵢ < 1 (mass-ratio domain bounds).
		// This is a regression test for the BFGS y-vector bug: previously
		// y was computed as ∇L(x_new,λ_new) − ∇L(x_old,λ_old), mixing two
		// different multipliers; the cross-term (λ_new−λ_old)·∇g caused ‖y‖
		// to explode (up to 10¹⁵), corrupting B and producing a singular KKT.
		const vCh = 12.084
		W := []float64{3.0550, 3.1771, 3.3894}
		alpha := []float64{0.0990, 0.0934, 0.0979}
		n := 3

		f := func(x []float64) float64 {
			prod := 1.0
			for i := range n {
				prod *= (1 - alpha[i]) / (x[i] - alpha[i])
			}
			return prod
		}
		eq := func(x []float64) float64 {
			sum := vCh
			for i := range n {
				sum += W[i] * math.Log(x[i])
			}
			return sum
		}
		bounds := make([]func([]float64) float64, 2*n)
		for i := range n {
			i := i
			bounds[i] = func(x []float64) float64 { return x[i] - alpha[i] }
			bounds[n+i] = func(x []float64) float64 { return 1 - x[i] }
		}

		// Initial guess: equal delta-v split across stages.
		x0 := make([]float64, n)
		for i := range n {
			x0[i] = math.Exp(-vCh / (float64(n) * W[i]))
		}

		dxRocket := []float64{1e-5, 1e-5, 1e-5}
		result, err := numericalanalysis.SLSQP(f, x0, dxRocket,
			[]func([]float64) float64{eq}, bounds, 1e-6, 500)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		// Equality constraint satisfied.
		if math.Abs(eq(result)) > 1e-4 {
			t.Errorf("equality constraint violated: g=%v", eq(result))
		}
		// All mass ratios inside (αᵢ, 1).
		for i := range n {
			if result[i] <= alpha[i] {
				t.Errorf("result[%d]=%v <= alpha[%d]=%v", i, result[i], i, alpha[i])
			}
			if result[i] >= 1 {
				t.Errorf("result[%d]=%v >= 1", i, result[i])
			}
		}
		// Objective must not increase from the initial feasible point
		// (first-order optimality: we at least found a descent).
		if f(result) >= f(x0) {
			t.Errorf("objective did not decrease: f(result)=%v, f(x0)=%v", f(result), f(x0))
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
