package numericalanalysis

import "math"

func BisectionExtremum(f Func1D, a, b float64, tol float64, max bool) float64 {
	for math.Abs(b-a) > tol {
		x := (a + b) / 2

		if f(x-tol) > f(x+tol) {
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

	return (a + b) / 2
}
