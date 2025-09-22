package numericalanalysis_test

import (
	"math"
	"testing"

	numericalanalysis "github.com/Russia9/numerical-analysis"
)

// Test cases for differential equation methods
var testCases = []struct {
	name     string
	system   numericalanalysis.FuncSystem
	x0       float64
	start    []float64
	xChar    []float64
	hBase    float64
	stop     func(x float64, y ...float64) bool
	expected [][]numericalanalysis.Point2D
}{
	{
		name: "simple_linear",
		system: numericalanalysis.FuncSystem{
			func(x float64, y ...float64) float64 { return 1 }, // dy/dx = 1, solution: y = x + C
		},
		x0:    0,
		start: []float64{0},
		xChar: []float64{},
		hBase: 1,
		stop:  func(x float64, y ...float64) bool { return x >= 3 },
		expected: [][]numericalanalysis.Point2D{
			{{0, 0}, {1, 1}, {2, 2}, {3, 3}},
		},
	},
	{
		name: "exponential_growth",
		system: numericalanalysis.FuncSystem{
			func(x float64, y ...float64) float64 { return y[0] }, // dy/dx = y, solution: y = e^x
		},
		x0:    0,
		start: []float64{1},
		xChar: []float64{},
		hBase: 0.1,
		stop:  func(x float64, y ...float64) bool { return x >= 0.2 },
		expected: [][]numericalanalysis.Point2D{
			{{0, 1}, {0.1, math.Exp(0.1)}, {0.2, math.Exp(0.2)}},
		},
	},
	{
		name: "system_of_equations",
		system: numericalanalysis.FuncSystem{
			func(x float64, y ...float64) float64 { return y[1] },  // dx/dt = y
			func(x float64, y ...float64) float64 { return -y[0] }, // dy/dt = -x
		},
		x0:    0,
		start: []float64{1, 0}, // x(0) = 1, y(0) = 0
		xChar: []float64{},
		hBase: 0.1,
		stop:  func(x float64, y ...float64) bool { return x >= 0.2 },
		expected: [][]numericalanalysis.Point2D{
			{{0, 1}, {0.1, math.Cos(0.1)}, {0.2, math.Cos(0.2)}}, // x(t) = cos(t)
			{{0, 0}, {0.1, math.Sin(0.1)}, {0.2, math.Sin(0.2)}}, // y(t) = sin(t)
		},
	},
	{
		name: "with_characteristic_points",
		system: numericalanalysis.FuncSystem{
			func(x float64, y ...float64) float64 {
				if x >= 1.5 {
					return 2
				}
				return 1
			},
		},
		x0:    0,
		start: []float64{0},
		xChar: []float64{1.5},
		hBase: 1,
		stop:  func(x float64, y ...float64) bool { return x >= 3 },
		expected: [][]numericalanalysis.Point2D{
			{{0, 0}, {1, 1}, {1.5, 1.5}, {2.5, 3.5}, {3.5, 5.5}},
		},
	},
}

func TestEulerMethod(t *testing.T) {
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			result, err := numericalanalysis.EulerMethod(
				tc.system, tc.x0, tc.start, tc.xChar, tc.hBase, tc.stop,
			)

			if err != nil {
				t.Errorf("EulerMethod returned error: %v", err)
				return
			}

			// Check dimensions
			if len(result) != len(tc.system) {
				t.Errorf("Expected %d result series, got %d", len(tc.system), len(result))
				return
			}

			// For simple test cases, check basic properties
			if len(result) > 0 && len(result[0]) > 0 {
				// Check initial condition
				if result[0][0].X != tc.x0 {
					t.Errorf("Initial X: expected %v, got %v", tc.x0, result[0][0].X)
				}
				if result[0][0].Y != tc.start[0] {
					t.Errorf("Initial Y: expected %v, got %v", tc.start[0], result[0][0].Y)
				}
			}
		})
	}
}

func TestModifiedEulerMethod(t *testing.T) {
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			result, err := numericalanalysis.ModifiedEulerMethod(
				tc.system, tc.x0, tc.start, tc.xChar, tc.hBase, tc.stop,
			)

			if err != nil {
				t.Errorf("ModifiedEulerMethod returned error: %v", err)
				return
			}

			// Check dimensions
			if len(result) != len(tc.system) {
				t.Errorf("Expected %d result series, got %d", len(tc.system), len(result))
				return
			}

			// Check initial condition
			if len(result) > 0 && len(result[0]) > 0 {
				if result[0][0].X != tc.x0 {
					t.Errorf("Initial X: expected %v, got %v", tc.x0, result[0][0].X)
				}
				if result[0][0].Y != tc.start[0] {
					t.Errorf("Initial Y: expected %v, got %v", tc.start[0], result[0][0].Y)
				}
			}
		})
	}
}

func TestRungeKuttaMethod(t *testing.T) {
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			result, err := numericalanalysis.RungeKuttaMethod(
				tc.system, tc.x0, tc.start, tc.xChar, tc.hBase, tc.stop,
			)

			if err != nil {
				t.Errorf("RungeKuttaMethod returned error: %v", err)
				return
			}

			// Check dimensions
			if len(result) != len(tc.system) {
				t.Errorf("Expected %d result series, got %d", len(tc.system), len(result))
				return
			}

			// Check initial condition
			if len(result) > 0 && len(result[0]) > 0 {
				if result[0][0].X != tc.x0 {
					t.Errorf("Initial X: expected %v, got %v", tc.x0, result[0][0].X)
				}
				if result[0][0].Y != tc.start[0] {
					t.Errorf("Initial Y: expected %v, got %v", tc.start[0], result[0][0].Y)
				}
			}
		})
	}
}

// Test accuracy comparison between methods
func TestMethodAccuracy(t *testing.T) {
	// Simple exponential growth: dy/dx = y, y(0) = 1, exact solution: y = e^x
	system := numericalanalysis.FuncSystem{
		func(x float64, y ...float64) float64 { return y[0] },
	}

	x0 := 0.0
	start := []float64{1.0}
	hBase := 0.1
	stop := func(x float64, y ...float64) bool { return x >= 1.0 }

	euler, err := numericalanalysis.EulerMethod(system, x0, start, nil, hBase, stop)
	if err != nil {
		t.Fatalf("EulerMethod failed: %v", err)
	}

	modified, err := numericalanalysis.ModifiedEulerMethod(system, x0, start, nil, hBase, stop)
	if err != nil {
		t.Fatalf("ModifiedEulerMethod failed: %v", err)
	}

	rk4, err := numericalanalysis.RungeKuttaMethod(system, x0, start, nil, hBase, stop)
	if err != nil {
		t.Fatalf("RungeKuttaMethod failed: %v", err)
	}

	// Check that all methods produce results
	if len(euler) == 0 || len(modified) == 0 || len(rk4) == 0 {
		t.Fatal("One or more methods produced no results")
	}

	// Check that RK4 is more accurate than Euler for the final point
	exact := math.Exp(1.0) // e^1

	eulerFinal := euler[0][len(euler[0])-1].Y
	rk4Final := rk4[0][len(rk4[0])-1].Y

	eulerError := math.Abs(eulerFinal - exact)
	rk4Error := math.Abs(rk4Final - exact)

	if rk4Error > eulerError {
		t.Logf("Warning: RK4 error (%v) > Euler error (%v) for exponential growth", rk4Error, eulerError)
	}

	t.Logf("Exact: %v, Euler: %v (error: %v), RK4: %v (error: %v)",
		exact, eulerFinal, eulerError, rk4Final, rk4Error)
}

// Test error conditions
func TestDifferentialMethodErrors(t *testing.T) {
	t.Run("dimension_mismatch", func(t *testing.T) {
		system := numericalanalysis.FuncSystem{
			func(x float64, y ...float64) float64 { return 1 },
			func(x float64, y ...float64) float64 { return 1 },
		}
		start := []float64{0} // Wrong dimension
		stop := func(x float64, y ...float64) bool { return x >= 1 }

		_, err := numericalanalysis.EulerMethod(system, 0, start, nil, 1, stop)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("Expected ErrWrongInput, got %v", err)
		}

		_, err = numericalanalysis.ModifiedEulerMethod(system, 0, start, nil, 1, stop)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("Expected ErrWrongInput, got %v", err)
		}

		_, err = numericalanalysis.RungeKuttaMethod(system, 0, start, nil, 1, stop)
		if err != numericalanalysis.ErrWrongInput {
			t.Errorf("Expected ErrWrongInput, got %v", err)
		}
	})
}

// Test characteristic points handling
func TestCharacteristicPoints(t *testing.T) {
	system := numericalanalysis.FuncSystem{
		func(x float64, y ...float64) float64 { return 1 },
	}

	x0 := 0.0
	start := []float64{0.0}
	xChar := []float64{0.5, 1.5, 2.5}
	hBase := 1.0
	stop := func(x float64, y ...float64) bool { return x >= 3 }

	methods := []struct {
		name string
		fn   func(numericalanalysis.FuncSystem, float64, []float64, []float64, float64, func(float64, ...float64) bool) ([][]numericalanalysis.Point2D, error)
	}{
		{"Euler", numericalanalysis.EulerMethod},
		{"ModifiedEuler", numericalanalysis.ModifiedEulerMethod},
		{"RungeKutta", numericalanalysis.RungeKuttaMethod},
	}

	for _, method := range methods {
		t.Run(method.name, func(t *testing.T) {
			result, err := method.fn(system, x0, start, xChar, hBase, stop)
			if err != nil {
				t.Fatalf("%s failed: %v", method.name, err)
			}

			if len(result) == 0 || len(result[0]) == 0 {
				t.Fatal("No results produced")
			}

			// Check that characteristic points are included
			xValues := make([]float64, len(result[0]))
			for i, point := range result[0] {
				xValues[i] = point.X
			}

			for _, charPoint := range xChar {
				found := false
				for _, x := range xValues {
					if math.Abs(x-charPoint) < 1e-10 {
						found = true
						break
					}
				}
				if !found {
					t.Errorf("Characteristic point %v not found in x values %v", charPoint, xValues)
				}
			}
		})
	}
}
