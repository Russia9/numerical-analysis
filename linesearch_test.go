package numericalanalysis_test

import (
	"math"
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

func TestLineSearchArmijo(t *testing.T) {
	t.Run("quadratic descent", func(t *testing.T) {
		// φ(α) = (1 - α)²  at x=1, direction d=-1
		// φ(0)=1, φ'(0)=-2
		phi := func(a float64) float64 { return (1 - a) * (1 - a) }
		alpha, err := numericalanalysis.LineSearchArmijo(phi, 1.0, -2.0, 1.0, 0.1, 50)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		// Armijo must be satisfied
		if phi(alpha) > 1.0+0.1*alpha*(-2.0) {
			t.Errorf("Armijo condition not met at alpha=%v", alpha)
		}
	})

	t.Run("already satisfied at alpha0", func(t *testing.T) {
		phi := func(a float64) float64 { return 1 - a } // linear descent
		alpha, err := numericalanalysis.LineSearchArmijo(phi, 1.0, -1.0, 1.0, 0.1, 10)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if math.Abs(alpha-1.0) > 1e-10 {
			t.Errorf("expected alpha=1.0, got %v", alpha)
		}
	})

	t.Run("wrong input - non-descent direction", func(t *testing.T) {
		phi := func(a float64) float64 { return a }
		_, err := numericalanalysis.LineSearchArmijo(phi, 0.0, 1.0, 1.0, 0.1, 10)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (dphi0 >= 0)", err)
		}
	})

	t.Run("wrong input - c1 out of range", func(t *testing.T) {
		phi := func(a float64) float64 { return (1 - a) * (1 - a) }
		_, err := numericalanalysis.LineSearchArmijo(phi, 1.0, -2.0, 1.0, 0.0, 10)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (c1=0)", err)
		}
		_, err = numericalanalysis.LineSearchArmijo(phi, 1.0, -2.0, 1.0, 1.0, 10)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (c1=1)", err)
		}
	})

	t.Run("wrong input - alpha0 non-positive", func(t *testing.T) {
		phi := func(a float64) float64 { return (1 - a) * (1 - a) }
		_, err := numericalanalysis.LineSearchArmijo(phi, 1.0, -2.0, 0.0, 0.1, 10)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (alpha0=0)", err)
		}
	})

	t.Run("wrong input - maxIter non-positive", func(t *testing.T) {
		phi := func(a float64) float64 { return (1 - a) * (1 - a) }
		_, err := numericalanalysis.LineSearchArmijo(phi, 1.0, -2.0, 1.0, 0.1, 0)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (maxIter=0)", err)
		}
	})

	t.Run("ErrDidNotConverge", func(t *testing.T) {
		// phi never satisfies Armijo with only 1 iteration: halve once, still no good
		// Use a function where the condition is never met in 1 step
		phi := func(a float64) float64 { return 1 - 0.001*a } // tiny decrease, c1=0.5 requires 0.5*a*(-1)
		_, err := numericalanalysis.LineSearchArmijo(phi, 1.0, -1.0, 1.0, 0.5, 1)
		if err != numericalanalysis.ErrDidNotConverge {
			t.Errorf("err = %v, want ErrDidNotConverge", err)
		}
	})
}

func TestLineSearchWolfe(t *testing.T) {
	t.Run("quadratic", func(t *testing.T) {
		// f(x) = x², at x=2, direction d=-1
		// φ(α) = (2-α)², φ'(α) = -2(2-α)
		// φ(0)=4, φ'(0)=-4
		phi := func(a float64) float64 { x := 2 - a; return x * x }
		dphi := func(a float64) float64 { return -2 * (2 - a) }
		alpha, err := numericalanalysis.LineSearchWolfe(phi, dphi, 4.0, -4.0, 1.0, 1e-4, 0.9, 30)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		// Verify strong Wolfe conditions
		phiA := phi(alpha)
		dphiA := dphi(alpha)
		if phiA > 4.0+1e-4*alpha*(-4.0) {
			t.Errorf("sufficient decrease not met at alpha=%v: phi=%v", alpha, phiA)
		}
		if math.Abs(dphiA) > 0.9*math.Abs(-4.0) {
			t.Errorf("curvature condition not met at alpha=%v: dphi=%v", alpha, dphiA)
		}
	})

	t.Run("shifted quadratic satisfies conditions", func(t *testing.T) {
		// φ(α) = (α-2)² + 1, minimum at α=2; φ(0)=5, φ'(0)=-4
		phi := func(a float64) float64 { return (a-2)*(a-2) + 1 }
		dphi := func(a float64) float64 { return 2 * (a - 2) }
		c1, c2 := 1e-4, 0.9
		alpha, err := numericalanalysis.LineSearchWolfe(phi, dphi, 5.0, -4.0, 1.0, c1, c2, 50)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if phi(alpha) > 5.0+c1*alpha*(-4.0) {
			t.Errorf("sufficient decrease not met at alpha=%v", alpha)
		}
		if math.Abs(dphi(alpha)) > c2*math.Abs(-4.0) {
			t.Errorf("curvature condition not met at alpha=%v", alpha)
		}
	})

	t.Run("wrong input - non-descent direction", func(t *testing.T) {
		phi := func(a float64) float64 { return a }
		dphi := func(a float64) float64 { return 1.0 }
		_, err := numericalanalysis.LineSearchWolfe(phi, dphi, 0.0, 1.0, 1.0, 1e-4, 0.9, 10)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (dphi0 >= 0)", err)
		}
	})

	t.Run("wrong input - c1 >= c2", func(t *testing.T) {
		phi := func(a float64) float64 { return (2 - a) * (2 - a) }
		dphi := func(a float64) float64 { return -2 * (2 - a) }
		_, err := numericalanalysis.LineSearchWolfe(phi, dphi, 4.0, -4.0, 1.0, 0.9, 0.5, 10)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("err = %v, want ErrWrongInput (c1 >= c2)", err)
		}
	})

	t.Run("ErrDidNotConverge", func(t *testing.T) {
		// phi(a) = -10a + sin(a): phi(0)=0, phi'(0)=-9 < 0 (descent).
		// |phi'(a)| = |cos(a)-10| >= 9 > 0.9*9 = 8.1, so curvature is never met.
		// With maxIter=2, the outer loop exhausts before finding a valid step.
		phi := func(a float64) float64 { return -10*a + math.Sin(a) }
		dphi := func(a float64) float64 { return -10 + math.Cos(a) }
		_, err := numericalanalysis.LineSearchWolfe(phi, dphi, 0, -9, 1.0, 1e-4, 0.9, 2)
		if err != numericalanalysis.ErrDidNotConverge {
			t.Errorf("err = %v, want ErrDidNotConverge", err)
		}
	})
}
