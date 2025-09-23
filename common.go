package numericalanalysis

import "errors"

type Point2D struct {
	X float64
	Y float64
}

type Func1D func(x float64) float64

type Func2D func(x, y float64) float64

type FuncND func(x float64, y ...float64) float64

type Point3D struct {
	X float64
	Y float64
	Z float64
}

var ErrNoSolution = errors.New("no solution")

var ErrWrongInput = errors.New("wrong input")
