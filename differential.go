package numericalanalysis

type FuncSystem = []func(x float64, y ...float64) float64

func EulerMethod(f FuncSystem, start []float64, x0, xN, h float64) ([][]Point2D, error) {
	// Check input dimensions
	if len(f) != len(start) {
		return nil, ErrWrongInput
	}

	// Calculate number of steps
	n := int((xN - x0) / h)

	// Initialize result
	result := make([][]float64, n)
	result[0] = start

	// Iterate
	for i := 1; i < n; i++ {
		// Initialize new point
		result[i] = make([]float64, len(f))

		// Calculate new point
		for j := 0; j < len(f); j++ {
			result[i][j] = result[i-1][j] + h*f[j](x0+float64(i-1)*h, result[i-1]...)
		}
	}

	// Convert to Point2D
	points := make([][]Point2D, len(f))
	for i := 0; i < len(f); i++ {
		points[i] = make([]Point2D, n)
		for j := 0; j < n; j++ {
			points[i][j] = Point2D{X: x0 + float64(j)*h, Y: result[j][i]}
		}
	}

	return points, nil
}

func ModifiedEulerMethod(f FuncSystem, start []float64, x0, xN, h float64) ([][]Point2D, error) {
	// Check input dimensions
	if len(f) != len(start) {
		return nil, ErrWrongInput
	}

	// Calculate number of steps
	n := int((xN - x0) / h)

	// Initialize result
	result := make([][]float64, n)
	result[0] = start

	// Iterate
	for i := 1; i < n; i++ {
		// Initialize new point
		result[i] = make([]float64, len(f))

		// Calculate k1
		k1 := make([]float64, len(f))
		for j := 0; j < len(f); j++ {
			k1[j] = f[j](x0+float64(i-1)*h, result[i-1]...)
		}

		// Calculate k2
		k2 := make([]float64, len(f))
		for j := 0; j < len(f); j++ {
			// Calculate input
			inp := make([]float64, len(f))
			for k := 0; k < len(f); k++ {
				inp[k] = result[i-1][k] + h*k1[k]
			}

			k2[j] = f[j](x0+float64(i)*h, inp...)
		}

		// Calculate new point
		for j := 0; j < len(f); j++ {
			result[i][j] = result[i-1][j] + h*(k1[j]+k2[j])/2
		}
	}

	// Convert to Point2D
	points := make([][]Point2D, len(f))
	for i := 0; i < len(f); i++ {
		points[i] = make([]Point2D, n)
		for j := 0; j < n; j++ {
			points[i][j] = Point2D{X: x0 + float64(j)*h, Y: result[j][i]}
		}
	}

	return points, nil
}

func RungeKuttaMethod(f FuncSystem, start []float64, x0, xN, h float64) ([][]Point2D, error) {
	// Check input dimensions
	if len(f) != len(start) {
		return nil, ErrWrongInput
	}

	// Calculate number of steps
	n := int((xN - x0) / h)

	// Initialize result
	result := make([][]float64, n)
	result[0] = start

	// Iterate
	for i := 1; i < n; i++ {
		// Initialize new point
		result[i] = make([]float64, len(f))

		// Calculate k1
		k1 := make([]float64, len(f))
		for j := 0; j < len(f); j++ {
			k1[j] = h * f[j](x0+float64(i-1)*h, result[i-1]...)
		}

		// Calculate k2
		k2 := make([]float64, len(f))
		for j := 0; j < len(f); j++ {
			// Calculate input
			inp := make([]float64, len(f))
			for k := 0; k < len(f); k++ {
				inp[k] = result[i-1][k] + k1[k]/2
			}

			k2[j] = h * f[j](x0+float64(i-1)*h+h/2, inp...)
		}

		// Calculate k3
		k3 := make([]float64, len(f))
		for j := 0; j < len(f); j++ {
			// Calculate input
			inp := make([]float64, len(f))
			for k := 0; k < len(f); k++ {
				inp[k] = result[i-1][k] + k2[k]/2
			}

			k3[j] = h * f[j](x0+float64(i-1)*h+h/2, inp...)
		}

		// Calculate k4
		k4 := make([]float64, len(f))
		for j := 0; j < len(f); j++ {
			// Calculate input
			inp := make([]float64, len(f))
			for k := 0; k < len(f); k++ {
				inp[k] = result[i-1][k] + k3[k]
			}

			k4[j] = h * f[j](x0+float64(i)*h, inp...)
		}

		// Calculate new point
		for j := 0; j < len(f); j++ {
			result[i][j] = result[i-1][j] + (k1[j]+2*k2[j]+2*k3[j]+k4[j])/6
		}
	}

	// Convert to Point2D
	points := make([][]Point2D, len(f))
	for i := 0; i < len(f); i++ {
		points[i] = make([]Point2D, n)
		for j := 0; j < n; j++ {
			points[i][j] = Point2D{X: x0 + float64(j)*h, Y: result[j][i]}
		}
	}

	return points, nil
}
