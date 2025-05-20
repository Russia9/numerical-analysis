package numericalanalysis

func LagrangeInterpolation1D(points []Point2D) Func1D {
	n := len(points)

	// Lagrange basis polynomials
	l := make([]Func1D, n)

	for i := range n {
		l[i] = func(x float64) float64 {
			result := 1.0
			for j := range n {
				if j != i {
					result *= (x - points[j].X) / (points[i].X - points[j].X)
				}
			}
			return result * points[i].Y
		}
	}

	// Lagrange interpolation
	return func(x float64) float64 {
		result := 0.0
		for i := range n {
			result += l[i](x)
		}
		return result
	}
}

func LinearInterpolation1D(points []Point2D) Func1D {
	n := len(points)

	return func(x float64) float64 {
		i := 0
		if x < points[0].X {
			i = 0
		} else if x >= points[n-1].X {
			i = n - 2
		} else {
			for i = range n - 1 {
				if points[i].X <= x && x < points[i+1].X {
					break
				}
			}
		}

		return points[i].Y + (points[i+1].Y-points[i].Y)/(points[i+1].X-points[i].X)*(x-points[i].X)
	}
}

func QuadraticInterpolation1D(points []Point2D) Func1D {
	n := len(points)
	if n < 3 {
		return nil
	}

	return func(x float64) float64 {
		i := 1
		if x < points[1].X {
			i = 1
		} else if x >= points[n-2].X {
			i = n - 2
		} else {
			// Find center point of the point region
			for j := 1; j < n-1; j += 2 {
				if points[j-1].X <= x && x < points[j+1].X {
					i = j
					break
				}
			}
		}

		return points[i-1].Y*(x-points[i].X)*(x-points[i+1].X)/((points[i-1].X-points[i].X)*(points[i-1].X-points[i+1].X)) +
			points[i].Y*(x-points[i-1].X)*(x-points[i+1].X)/((points[i].X-points[i-1].X)*(points[i].X-points[i+1].X)) +
			points[i+1].Y*(x-points[i-1].X)*(x-points[i].X)/((points[i+1].X-points[i-1].X)*(points[i+1].X-points[i].X))
	}
}
