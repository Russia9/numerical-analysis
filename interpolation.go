package numericalanalysis

func LagrangeInterpolation(points []Point2D) Func1D {
	n := len(points)

	// Lagrange basis polynomials
	l := make([]Func1D, n)

	for i := 0; i < n; i++ {
		l[i] = func(x float64) float64 {
			result := 1.0
			for j := 0; j < n; j++ {
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
		for i := 0; i < n; i++ {
			result += l[i](x)
		}
		return result
	}
}
