package main

import (
	"image/color"

	numericalanalysis "github.com/Russia9/numerical-analysis"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
)

func main() {
	p := plot.New()

	p.Title.Text = "Lagrangian Interpolation test"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"
	p.Add(plotter.NewGrid())
	p.X.Min = 0
	p.X.Max = 10
	p.Y.Min = -10
	p.Y.Max = 10

	points := []numericalanalysis.Point2D{
		{X: 1, Y: 2},
		{X: 4, Y: 3},
		{X: 5, Y: -2},
		{X: 7, Y: 2},
	}

	interpFunc := numericalanalysis.LagrangeInterpolation(points)

	// Create test table
	for i := 0; i < 10; i++ {
		x := float64(i)
		y := interpFunc(x)
		println(x, y)
	}

	interp := plotter.NewFunction(numericalanalysis.LagrangeInterpolation(points))
	interp.Color = color.RGBA{B: 255, A: 255}

	p.Add(interp)

	xys := make(plotter.XYs, len(points))
	for i, point := range points {
		xys[i] = struct{ X, Y float64 }{point.X, point.Y}
	}
	plotutil.AddLinePoints(p, xys)

	p.Save(400, 400, "lagrange.png")
}
