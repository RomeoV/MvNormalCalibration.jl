# MvNormalCalibration.jl
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
