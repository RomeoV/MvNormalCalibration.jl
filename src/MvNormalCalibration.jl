module MvNormalCalibration
using Distributions, LinearAlgebra
using SpecialFunctions

export CentralPredictionSet
export computecalibration, sharpness

greet() = print("Hello World!")

abstract type PredictionSetType end

struct CentralPredictionSet <: PredictionSetType end


computecalibration( preds::AbstractVector{<:AbstractMvNormal},
                    truevals::AbstractVector{<:AbstractVector{<:Real}};
                    kwargs...) = computecalibration(CentralPredictionSet, preds, truevals; kwargs...)
function computecalibration(::Type{CentralPredictionSet},
                            preds::AbstractVector{<:AbstractMvNormal},
                            truevals::AbstractVector{<:AbstractVector{<:Real}};
                            pvals=0.05:0.05:0.95)
  @assert allequal(length.(truevals))
  d = length(first(truevals))
  @assert all(0 .< pvals .< 1)

  distthresholds = quantile.([Normal()], 1/2 .+ pvals.^(1/d)/2)
  containments = map(zip(preds, truevals)) do (pred, trueval)
    Λvec, Q = eigen(pred.Σ)
    normalizedtrueval = Q*(trueval - pred.μ) ./ sqrt.(Λvec)
    [all(abs.(normalizedtrueval) .< thresh) for thresh in distthresholds]
  end

  calibrationvals = mean(stack(containments; dims=1); dims=1)[:]
  (; pvals, calibrationvals)
end

# Measures volume of 1std ellipsoid.
function sharpness(pred::AbstractMvNormal)
    Λvec, _Q = eigen(pred.Σ)
    d = length(pred)
    # https://en.wikipedia.org/wiki/Ellipsoid#Standard_equation_2
    π^(d/2)/gamma(d/2+1) * prod(sqrt.(Λvec))
end

end # module MvNormalCalibration
