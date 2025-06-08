package numericalanalysis

import "sort"

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

	// Sort points by X
	sort.Slice(points, func(i, j int) bool {
		return points[i].X < points[j].X
	})

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

func BilinearInterpolation2D(points []Point3D) Func2D {
	// Convert points to map[y]Point2D(x,z)
	rows := map[float64][]Point2D{}
	var yValues []float64
	for _, pt := range points {
		rows[pt.Y] = append(rows[pt.Y], Point2D{X: pt.X, Y: pt.Z})
	}

	// Sort Y values
	for y := range rows {
		yValues = append(yValues, y)
	}
	sort.Float64s(yValues)

	// Sort each row by X
	for _, row := range rows {
		sort.Slice(row, func(i, j int) bool {
			return row[i].X < row[j].X
		})
	}

	return func(x, y float64) float64 {
		// Find bounding Y values
		var y0, y1 float64
		switch {
		case y <= yValues[0]:
			y0, y1 = yValues[0], yValues[1]
		case y >= yValues[len(yValues)-1]:
			y0, y1 = yValues[len(yValues)-2], yValues[len(yValues)-1]
		default:
			for i := range len(yValues) - 1 {
				if yValues[i] <= y && y < yValues[i+1] {
					y0, y1 = yValues[i], yValues[i+1]
					break
				}
			}
		}

		// Interpolate along x for each y-row
		f0 := LinearInterpolation1D(rows[y0])
		f1 := LinearInterpolation1D(rows[y1])

		z0 := f0(x)
		z1 := f1(x)

		// Interpolate along y between z0 and z1
		return z0 + (z1-z0)*(y-y0)/(y1-y0)
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
