# numerical-analysis

Go library (`package numericalanalysis`, module `github.com/Russia9/numerical-analysis`) implementing numerical methods for coursework/study.

## Files

| File | Contents |
|------|----------|
| `common.go` | Shared types: `Point2D`, `Point3D`, `Func1D`, `Func2D`, `FuncND`; sentinel errors: `ErrNoSolution`, `ErrWrongInput`, `ErrSingularMatrix`, `ErrDidNotConverge` |
| `matrix.go` | `Matrix` type (`[][]float64`); `Norm`, `IdentityMatrix`, `Column`, `Equal`, `Add`, `Transpose`, `Det` (recursive cofactor), `Mul`, `MulNumber`, `Inverse` (adjugate/det) |
| `sle.go` | `Cramer` — linear system solver via Cramer's rule |
| `bisection.go` | `BisectionExtremum` (finds min/max of 1D function), `BisectionValue` (finds x where f(x)=value) |
| `newton.go` | `DampedNewtonExtremum` (multivariable minimisation, Hessian + Levenberg–Marquardt-style backtracking); `SENewton` (Newton for nonlinear systems) |
| `sde.go` | ODE solvers: `EulerMethod`, `ModifiedEulerMethod`, `RungeKuttaMethod` (RK4) |
| `integral.go` | `IntegralTrapezoid` (composite trapezoid, N intervals), `IntegralChebychev` (Chebyshev-Gauss quadrature) |
| `interpolation.go` | `LagrangeInterpolation1D`, `LinearInterpolation1D`, `QuadraticInterpolation1D` (piecewise), `BilinearInterpolation2D` |

## Design conventions

- All solvers return `(result, error)` and validate inputs, returning sentinel errors from `common.go`.
- ODE solvers support **characteristic points** (`xChar []float64`): the solver lands exactly on these x-values and signals discontinuities via `fromRight bool` in `FuncND`.
- ODE solvers accept an adaptive **stop/half callback**: `(half=true, stop=false)` halves the step and retries; `stop=true` terminates.
- `DampedNewtonExtremum` regularises with `(H + αI)^{-1} ∇f`; α decreases on acceptance, increases on backtrack failure.
