package main

import (
	"image/color"

	numericalanalysis "github.com/Russia9/numerical-analysis"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
)

func main() {
	p := plot.New()

	p.Title.Text = "Interpolation test"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"
	p.Add(plotter.NewGrid())

	// --- Cxa(M)
	p.X.Min = 0
	p.X.Max = 11
	p.Y.Min = 0
	p.Y.Max = 1
	points := []numericalanalysis.Point2D{
		{X: 0.01, Y: 0.30},
		{X: 0.55, Y: 0.30},
		{X: 0.8, Y: 0.55},
		{X: 0.9, Y: 0.70},
		{X: 1.0, Y: 0.84},
		{X: 1.06, Y: 0.86},
		{X: 1.1, Y: 0.87},
		{X: 1.2, Y: 0.83},
		{X: 1.3, Y: 0.80},
		{X: 1.4, Y: 0.79},
		{X: 2.0, Y: 0.65},
		{X: 2.6, Y: 0.55},
		{X: 3.4, Y: 0.50},
		{X: 6.2, Y: 0.45},
		{X: 10.2, Y: 0.40},
	}

	// --- C^a_ya(M)
	// p.X.Min = 0
	// p.X.Max = 11
	// p.Y.Min = 0
	// p.Y.Max = 1
	// points := []numericalanalysis.Point2D{
	// 	{X: 0.01, Y: 0.25},
	// 	{X: 0.55, Y: 0.25},
	// 	{X: 0.8, Y: 0.25},
	// 	{X: 0.9, Y: 0.20},
	// 	{X: 1.0, Y: 0.30},
	// 	{X: 1.06, Y: 0.31},
	// 	{X: 1.1, Y: 0.25},
	// 	{X: 1.2, Y: 0.25},
	// 	{X: 1.3, Y: 0.25},
	// 	{X: 1.4, Y: 0.25},
	// 	{X: 2.0, Y: 0.25},
	// 	{X: 2.6, Y: 0.25},
	// 	{X: 3.4, Y: 0.25},
	// 	{X: 6.2, Y: 0.25},
	// 	{X: 10.2, Y: 0.25},
	// }

	// p.X.Min = 0
	// p.X.Max = 10
	// p.Y.Min = -4
	// p.Y.Max = 4
	// points := []numericalanalysis.Point2D{
	// 	{X: 0.0, Y: 0.0},
	// 	{X: 1.0, Y: 2.0},
	// 	{X: 2.0, Y: -1.0},
	// 	{X: 3.0, Y: 1.0},
	// 	{X: 2.0, Y: -1.0},
	// }

	lagrange := numericalanalysis.LagrangeInterpolation1D(points)
	linear := numericalanalysis.LinearInterpolation1D(points)
	quadratic := numericalanalysis.QuadraticInterpolation1D(points)

	lagrangePlotter := plotter.NewFunction(lagrange)
	lagrangePlotter.Color = color.RGBA{B: 255, A: 255}
	lagrangePlotter.Samples = 1000

	linearPlotter := plotter.NewFunction(linear)
	linearPlotter.Color = color.RGBA{R: 255, A: 255}
	linearPlotter.Samples = 1000

	quadraticPlotter := plotter.NewFunction(quadratic)
	quadraticPlotter.Color = color.RGBA{G: 255, A: 255}
	quadraticPlotter.Samples = 1000

	// p.Add(lagrangePlotter)
	p.Add(linearPlotter)
	p.Add(quadraticPlotter)

	// Create test table
	// for i := 0.0; i < 10; i += 0.1 {
	// 	x := float64(i)
	// 	println(x, lagrange(x), linear(x))
	// }

	xys := make(plotter.XYs, len(points))
	for i, point := range points {
		xys[i] = struct{ X, Y float64 }{point.X, point.Y}
	}
	scatter, err := plotter.NewScatter(xys)
	if err != nil {
		panic(err)
	}
	p.Add(scatter)

	p.Save(400, 400, "lagrange.png")
}
