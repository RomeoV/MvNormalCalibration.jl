# MvNormalCalibration.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://romeov.github.io/MvNormalCalibration.jl/dev/)
[![codecov](https://codecov.io/gh/RomeoV/ProbabilisticParameterEstimators.jl/graph/badge.svg?token=5J82UXPL8I)](https://codecov.io/gh/RomeoV/ProbabilisticParameterEstimators.jl)
[![Aqua](https://img.shields.io/badge/Tested_with-Aqua-turquoise)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/Tested_with-JET-violet)](https://github.com/aviatesk/JET.jl)


This package implements a simple way to measure calibration and sharpness for multivariate normal distributions.
In particular, it exports two functions with (roughly) signatures
- `computecalibration( preds::Vector{<:MvNormal}, truevals::Vector{Vector{Float64}} )`
- `sharpness( pred::MvNormal )`

`computecalibration(...)` iterates over `(prediction, trueval)` tuples and measures
whether `trueval` falls into a (central) prediction set for a series of coverage levels.
In particular, if the predictions are well calibrated, the following should hold:
```julia
(; pvals, calibrationvals) = computecalibration(preds, truevals)
@test pvals ≈ calibrationvals
```

`sharpness(...)` computes the (hyper-)volume of one standard deviation of a multivariate normal.
