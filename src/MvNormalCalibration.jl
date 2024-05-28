module MvNormalCalibration
using Distributions, LinearAlgebra
using SpecialFunctions

export CentralPredictionSet
export computecalibration, sharpness

abstract type PredictionSetType end

struct CentralPredictionSet <: PredictionSetType end

"""
    computecalibration(preds::AbstractVector{<:Normal}, truevals::AbstractVector{<:Real}; pvals)
    computecalibration(preds::AbstractVector{<:MvNormal}, truevals::AbstractVector{<:AbstractVector{<:Real}}; pvals)

Compute calibration for a series of predicted (uni- or multivariate) normal distributions
given a series of true observations using central prediction sets.

Returns a named tuple `(; pvals, calibrationvals)`.
If the predictions are well calibrated, then `pvals ≈ calibrationvals` for all indices.
Plotting `plot(pvals, calibrationvals)` should then give a straight line from (0, 0), to (1, 1).

# kwargs
- `pvals`: The probabilities to evaluate the coverage at. Defaults to `0:0.05:1`.
"""
function computecalibration(preds::AbstractVector{<:Normal},
        truevals::AbstractVector{<:Real};
        kwargs...)
    computecalibration(CentralPredictionSet, preds, truevals; kwargs...)
end
function computecalibration(::Type{CentralPredictionSet},
        preds::AbstractVector{<:Normal},
        truevals::AbstractVector{<:Real};
        pvals = 0.0:0.05:1.0)
    @assert length(preds) .== length(truevals)
    d = length(first(truevals))
    @assert all(0 .<= pvals .<= 1)

    dist_thresholds = quantile.([Normal()], 1 / 2 .+ pvals ./ 2)
    containments = map(zip(preds, truevals)) do (pred, trueval)
        normalized_trueval = (trueval - mean(pred)) / std(pred)
        abs(normalized_trueval) .<= dist_thresholds
    end

    # Average over samples.
    # The result is a vector with an average containment rate for each `pval`,
    # indicating what percentage of samples actually fell in the pvals[0], pvals[1], 
    # pvals[2], ... percentile prediction set.
    # If the predictions are well calibrated, then `pvals ≈ calibrationvals` for all indices.
    # Plotting `plot(pvals, calibrationvals)` should then give a straight line 
    # from (0, 0), to (1, 1).
    calibrationvals = mean(containments)
    @assert length(calibrationvals) == length(pvals)
    return (; pvals, calibrationvals)
end

function computecalibration(preds::AbstractVector{<:AbstractMvNormal},
        truevals::AbstractVector{<:AbstractVector{<:Real}};
        kwargs...)
    computecalibration(CentralPredictionSet, preds, truevals; kwargs...)
end
function computecalibration(::Type{CentralPredictionSet},
        preds::AbstractVector{<:AbstractMvNormal},
        truevals::AbstractVector{<:AbstractVector{<:Real}};
        pvals = 0.0:0.05:1.0)
    @assert allequal(length.(truevals))
    @assert length.(preds) == length.(truevals)
    d = length(first(truevals))
    @assert all(0 .<= pvals .<= 1)

    dist_thresholds = quantile.([Normal()], 1 / 2 .+ pvals .^ (1 / d) / 2)
    containments = map(zip(preds, truevals)) do (pred, trueval)
        Λvec, Q = eigen(cov(pred))
        normalized_trueval = Q * (trueval - mean(pred)) ./ sqrt.(Λvec)
        # all normalized values have to be in a given distance threshold
        all(abs.(normalized_trueval) .< dist_thresholds'; dims = 1)[:]
    end

    calibrationvals = mean(containments)
    (; pvals, calibrationvals)
end

"""
    sharpness(pred::Normal)
    sharpness(pred::AbstractMvNormal)

Compute sharpness for (uni- or multivariate) normal distribution, which we define as the
(hyper-)volume of the ellipsoid that contains one standard deviation of the distribution.

# Examples
```julia-repl
julia> sharpness(Normal())
2
julia> sharpness(MvNormal(zeros(2), I(2)))  # recall area of circle=πr^2
π
```
"""
function sharpness(pred::Normal)
    2 * std(pred)
end

function sharpness(pred::AbstractMvNormal)
    Λvec, _Q = eigen(cov(pred))
    d = length(pred)
    # https://en.wikipedia.org/wiki/Ellipsoid#Standard_equation_2
    π^(d / 2) / gamma(d / 2 + 1) * prod(sqrt.(Λvec))
end

end # module MvNormalCalibration
