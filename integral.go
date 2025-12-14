package numericalanalysis

import "math"

// IntegralTrapezoid calculates the integral of a function f(x) using the trapezoid rule.
func IntegralTrapezoid(f func(float64) float64, a, b float64, N int) float64 {
	h := (b - a) / float64(N)
	res := 0.5 * (f(a) + f(b))
	for i := 1; i < N; i++ {
		x := a + float64(i)*h
		res += f(x)
	}
	res *= h
	return res
}

// IntegralChebychev calculates the integral of a function f(x) using the Chebychev quadrature method.
func IntegralChebychev(f func(float64) float64, a, b float64, N int) float64 {
	res := 0.0
	C_i := math.Pi / float64(N)
	for i := 1; i <= N; i++ {
		t_i := math.Cos((2*float64(i) - 1) / (2 * float64(N)) * math.Pi)
		res += f((a+b)/2+(b-a)/2*t_i) * math.Sqrt(1-t_i*t_i)
	}
	res *= C_i * (b - a) / 2
	return res
}
