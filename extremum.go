package numericalanalysis

import (
	"math"
)

func BisectionExtremum(f Func1D, a, b float64, tol float64, max bool) (float64, error) {
	// Check input
	if a >= b || tol <= 0 {
		return 0, ErrWrongInput
	}

	// Loop until tolerance is met
	for math.Abs(b-a) > tol {
		// Calculate midpoint
		x := (a + b) / 2

		if f(x-tol) > f(x+tol) { // Check if the function value at x-tol is greater than at x+tol
			if max {
				b = x
			} else {
				a = x
			}
		} else {
			if max {
				a = x
			} else {
				b = x
			}
		}
	}

	return (a + b) / 2, nil
}

func BisectionValue(f Func1D, a, b, value float64, tol float64, backward bool) (float64, error) {
	// Check input
	if a >= b || tol <= 0 {
		return 0, ErrWrongInput
	}

	// Loop until tolerance is met
	for math.Abs(b-a) > tol {
		// Calculate midpoint
		x := (a + b) / 2

		if f(x) > value { // Check if the function value at x is greater than the target value
			if backward {
				b = x
			} else {
				a = x
			}
		} else {
			if backward {
				a = x
			} else {
				b = x
			}
		}
	}

	return (a + b) / 2, nil
}
