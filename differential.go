package numericalanalysis

import "math"

type FuncSystem = []FuncND

func EulerMethod(
	f FuncSystem,
	x0 float64,
	start []float64, // [y1(x0), y2(x0), ...]
	xChar []float64, // [x1, x2, ...] - Points are always included, regardless of the hBase. The step size is adjusted to ensure that.
	hBase float64,
	stop func(x float64, y ...float64) bool,
) ([][]Point2D, error) {
	// Check input dimensions
	if len(f) != len(start) {
		return nil, ErrWrongInput
	}

	// Initialize result
	x := make([]float64, 1)
	result := make([][]float64, 1)
	x[0] = x0
	result[0] = start
	i := 1

	// Iterate
	for {
		// Check if the last point satisfies the stop condition
		if stop(x[i-1], result[i-1]...) {
			break
		}

		// Initialize new point
		result = append(result, make([]float64, len(f)))
		h := hBase
		x = append(x, x[i-1]+h)

		// Check if we would step over any characteristic point
		for _, ch := range xChar {
			if x[i-1] < ch && ch < x[i] {
				h = ch - x[i-1]
				x[i] = ch
			} else if math.Abs(x[i-1]-ch) < 10e-6 { // or if the last point was characteristic
				h = x[i-2] + hBase - ch
				x[i] = x[i-2] + hBase
			}
		}

		// Calculate new point using Euler's method
		for j := range f {
			result[i][j] = result[i-1][j] + h*f[j](x[i-1], result[i-1]...)
		}

		i++
	}

	// Convert to Point2D
	points := make([][]Point2D, len(f))
	for i := range f {
		points[i] = make([]Point2D, len(result))
		for j := range result {
			points[i][j] = Point2D{X: x[j], Y: result[j][i]}
		}
	}

	return points, nil
}

func ModifiedEulerMethod(
	f FuncSystem,
	x0 float64,
	start []float64, // [y1(x0), y2(x0), ...]
	xChar []float64, // [x1, x2, ...] - Points are always included, regardless of the hBase. The step size is adjusted to ensure that.
	hBase float64,
	stop func(x float64, y ...float64) bool,
) ([][]Point2D, error) {
	// Check input dimensions
	if len(f) != len(start) {
		return nil, ErrWrongInput
	}

	// Initialize result
	x := make([]float64, 1)
	result := make([][]float64, 1)
	x[0] = x0
	result[0] = start
	i := 1

	// Iterate
	for {
		// Check if the last point satisfies the stop condition
		if stop(x[i-1], result[i-1]...) {
			break
		}

		// Initialize new point
		result = append(result, make([]float64, len(f)))
		h := hBase
		x = append(x, x[i-1]+h)

		// Check if we would step over any characteristic point
		for _, ch := range xChar {
			if x[i-1] < ch && ch < x[i] {
				h = ch - x[i-1]
				x[i] = ch
			} else if math.Abs(x[i-1]-ch) < 10e-6 { // or if the last point was characteristic
				h = x[i-2] + hBase - ch
				x[i] = x[i-2] + hBase
			}
		}

		// Calculate k1
		k1 := make([]float64, len(f))
		for j := range f {
			k1[j] = f[j](x[i-1], result[i-1]...)
		}

		// Calculate k2
		k2 := make([]float64, len(f))
		for j := range f {
			// Calculate input
			inp := make([]float64, len(f))
			for k := range f {
				inp[k] = result[i-1][k] + h*k1[k]
			}

			k2[j] = f[j](x[i-1], inp...)
		}

		// Calculate new point
		for j := range f {
			result[i][j] = result[i-1][j] + h*(k1[j]+k2[j])/2
		}

		i++
	}

	// Convert to Point2D
	points := make([][]Point2D, len(f))
	for i := range f {
		points[i] = make([]Point2D, len(result))
		for j := range result {
			points[i][j] = Point2D{X: x[j], Y: result[j][i]}
		}
	}

	return points, nil
}

func RungeKuttaMethod(
	f FuncSystem,
	x0 float64,
	start []float64, // [y1(x0), y2(x0), ...]
	xChar []float64, // [x1, x2, ...] - Points are always included, regardless of the hBase. The step size is adjusted to ensure that.
	hBase float64,
	stop func(x float64, y ...float64) bool,
) ([][]Point2D, error) {
	// Check input dimensions
	if len(f) != len(start) {
		return nil, ErrWrongInput
	}

	// Initialize result
	x := make([]float64, 1)
	result := make([][]float64, 1)
	x[0] = x0
	result[0] = start
	i := 1

	// Iterate
	for {
		// Check if the last point satisfies the stop condition
		if stop(x[i-1], result[i-1]...) {
			break
		}

		// Initialize new point
		result = append(result, make([]float64, len(f)))
		h := hBase
		x = append(x, x[i-1]+h)

		// Check if we would step over any characteristic point
		for _, ch := range xChar {
			if x[i-1] < ch && ch < x[i] {
				h = ch - x[i-1]
				x[i] = ch
			} else if math.Abs(x[i-1]-ch) < 10e-6 { // or if the last point was characteristic
				h = x[i-2] + hBase - ch
				x[i] = x[i-2] + hBase
			}
		}

		// Calculate k1
		k1 := make([]float64, len(f))
		for j := range f {
			k1[j] = h * f[j](x[i-1], result[i-1]...)
		}

		// Calculate k2
		k2 := make([]float64, len(f))
		for j := range f {
			// Calculate input
			inp := make([]float64, len(f))
			for k := range f {
				inp[k] = result[i-1][k] + k1[k]/2
			}

			k2[j] = h * f[j](x[i-1]+h/2, inp...)
		}

		// Calculate k3
		k3 := make([]float64, len(f))
		for j := range f {
			// Calculate input
			inp := make([]float64, len(f))
			for k := range f {
				inp[k] = result[i-1][k] + k2[k]/2
			}

			k3[j] = h * f[j](x[i-1]+h/2, inp...)
		}

		// Calculate k4
		k4 := make([]float64, len(f))
		for j := range f {
			// Calculate input
			inp := make([]float64, len(f))
			for k := range f {
				inp[k] = result[i-1][k] + k3[k]
			}

			k4[j] = h * f[j](x[i], inp...)
		}

		// Calculate new point
		for j := range f {
			result[i][j] = result[i-1][j] + (k1[j]+2*k2[j]+2*k3[j]+k4[j])/6
		}

		i++
	}

	// Convert to Point2D
	points := make([][]Point2D, len(f))
	for i := range f {
		points[i] = make([]Point2D, len(result))
		for j := range result {
			points[i][j] = Point2D{X: x[j], Y: result[j][i]}
		}
	}

	return points, nil
}
