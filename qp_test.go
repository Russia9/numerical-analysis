package numericalanalysis_test

import (
	"math"
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestQPEqualityConstrained(t *testing.T) {
	approx := func(a, b, tol float64) bool { return math.Abs(a-b) < tol }

	t.Run("unconstrained min x²+y²", func(t *testing.T) {
		// min x² + y²  →  x=y=0
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{0, 0}
		x, lam, err := numericalanalysis.QPEqualityConstrained(H, c, nil, nil)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(x[0], 0, 1e-10) || !approx(x[1], 0, 1e-10) {
			t.Errorf("x = %v, want [0 0]", x)
		}
		if len(lam) != 0 {
			t.Errorf("expected no multipliers")
		}
	})

	t.Run("min x²+y² s.t. x+y=1", func(t *testing.T) {
		// Solution: x=y=0.5, λ=-1
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{0, 0}
		A := numericalanalysis.Matrix{{1, 1}}
		b := []float64{1}
		x, lam, err := numericalanalysis.QPEqualityConstrained(H, c, A, b)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(x[0], 0.5, 1e-10) || !approx(x[1], 0.5, 1e-10) {
			t.Errorf("x = %v, want [0.5 0.5]", x)
		}
		if !approx(lam[0], -1, 1e-10) {
			t.Errorf("lambda = %v, want [-1]", lam)
		}
	})

	t.Run("min (x-1)²+(y-2)² s.t. x+y=1", func(t *testing.T) {
		// H=2I, c=[-2,-4]; KKT: 2x+λ=2, 2y+λ=4, x+y=1
		// Subtracting: x-y=-1; combined with x+y=1 → x=0, y=1
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{-2, -4}
		A := numericalanalysis.Matrix{{1, 1}}
		b := []float64{1}
		x, _, err := numericalanalysis.QPEqualityConstrained(H, c, A, b)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(x[0], 0, 1e-10) || !approx(x[1], 1, 1e-10) {
			t.Errorf("x = %v, want [0 1]", x)
		}
	})

	t.Run("wrong input — c length mismatch", func(t *testing.T) {
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{0, 0, 0} // wrong length
		_, _, err := numericalanalysis.QPEqualityConstrained(H, c, nil, nil)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})
}

func TestQPActiveSet(t *testing.T) {
	approx := func(a, b, tol float64) bool { return math.Abs(a-b) < tol }

	t.Run("no constraints — pure QP", func(t *testing.T) {
		// min (x-1)² + (y-2)²  →  x=1, y=2
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{-2, -4}
		x, _, _, err := numericalanalysis.QPActiveSet(H, c, nil, nil, nil, nil)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(x[0], 1, 1e-8) || !approx(x[1], 2, 1e-8) {
			t.Errorf("x = %v, want [1 2]", x)
		}
	})

	t.Run("equality only — min x²+y² s.t. x+y=1", func(t *testing.T) {
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{0, 0}
		A := numericalanalysis.Matrix{{1, 1}}
		b := []float64{1}
		x, _, _, err := numericalanalysis.QPActiveSet(H, c, A, b, nil, nil)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(x[0], 0.5, 1e-8) || !approx(x[1], 0.5, 1e-8) {
			t.Errorf("x = %v, want [0.5 0.5]", x)
		}
	})

	t.Run("inactive inequality — unconstrained min inside feasible region", func(t *testing.T) {
		// min (x-1)²+(y-2)²; unconstrained min (1,2) satisfies x+y >= 1 (3>=1)
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{-2, -4}
		AIneq := numericalanalysis.Matrix{{1, 1}}
		bIneq := []float64{1}
		x, _, mu, err := numericalanalysis.QPActiveSet(H, c, nil, nil, AIneq, bIneq)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(x[0], 1, 1e-8) || !approx(x[1], 2, 1e-8) {
			t.Errorf("x = %v, want [1 2]", x)
		}
		if mu[0] > 1e-10 { // inactive constraint → multiplier = 0
			t.Errorf("mu = %v, want 0 for inactive constraint", mu[0])
		}
	})

	t.Run("active inequality — min x²+y² s.t. x+y >= 1", func(t *testing.T) {
		// Unconstrained min is (0,0); constraint pushes to x=y=0.5
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{0, 0}
		AIneq := numericalanalysis.Matrix{{1, 1}}
		bIneq := []float64{1}
		x, _, mu, err := numericalanalysis.QPActiveSet(H, c, nil, nil, AIneq, bIneq)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if !approx(x[0], 0.5, 1e-8) || !approx(x[1], 0.5, 1e-8) {
			t.Errorf("x = %v, want [0.5 0.5]", x)
		}
		if mu[0] < 0 {
			t.Errorf("mu = %v, want >= 0 (active constraint)", mu[0])
		}
	})

	t.Run("wrong input — H size mismatch", func(t *testing.T) {
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}} // 2×2
		c := []float64{0, 0, 0}                       // len 3
		_, _, _, err := numericalanalysis.QPActiveSet(H, c, nil, nil, nil, nil)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput", err)
		}
	})

	t.Run("multiple inequality constraints", func(t *testing.T) {
		// min (x-3)²+(y-3)²  s.t. x<=4 (→ -x>=-4) and y<=4 (→ -y>=-4) and x+y<=5 (→ -x-y>=-5)
		H := numericalanalysis.Matrix{{2, 0}, {0, 2}}
		c := []float64{-6, -6}
		AIneq := numericalanalysis.Matrix{{-1, 0}, {0, -1}, {-1, -1}}
		bIneq := []float64{-4, -4, -5}
		x, _, _, err := numericalanalysis.QPActiveSet(H, c, nil, nil, AIneq, bIneq)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		// Unconstrained min (3,3): x+y=6 > 5, so constraint 3 is active → solution on x+y=5
		// On x+y=5: min (x-3)²+(5-x-3)² = (x-3)²+(2-x)² → x=2.5, y=2.5
		if !approx(x[0], 2.5, 1e-8) || !approx(x[1], 2.5, 1e-8) {
			t.Errorf("x = %v, want [2.5 2.5]", x)
		}
	})
}
