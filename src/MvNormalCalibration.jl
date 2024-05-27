module MvNormalCalibration
using Distributions, LinearAlgebra
using SpecialFunctions

export CentralPredictionSet
export computecalibration, sharpness

greet() = print("Hello World!")

abstract type PredictionSetType end

struct CentralPredictionSet <: PredictionSetType end

computecalibration( preds::AbstractVector{<:Normal},
                    truevals::AbstractVector{<:Real};
                    kwargs...) = computecalibration(CentralPredictionSet, preds, truevals; kwargs...)
function computecalibration(::Type{CentralPredictionSet},
                            preds::AbstractVector{<:Normal},
                            truevals::AbstractVector{<:Real};
                            pvals=0.0:0.05:1.0)
  @assert length(preds) .== length(truevals)
  d = length(first(truevals))
  @assert all(0 .<= pvals .<= 1)

  dist_thresholds = quantile.([Normal()], 1/2 .+ pvals./2)
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



computecalibration( preds::AbstractVector{<:AbstractMvNormal},
                    truevals::AbstractVector{<:AbstractVector{<:Real}};
                    kwargs...) = computecalibration(CentralPredictionSet, preds, truevals; kwargs...)
function computecalibration(::Type{CentralPredictionSet},
                            preds::AbstractVector{<:AbstractMvNormal},
                            truevals::AbstractVector{<:AbstractVector{<:Real}};
                            pvals=0.0:0.05:1.0)
    @assert allequal(length.(truevals))
    @assert length.(preds) == length.(truevals)
    d = length(first(truevals))
    @assert all(0 .<= pvals .<= 1)
  
    dist_thresholds = quantile.([Normal()], 1/2 .+ pvals.^(1/d)/2)
    containments = map(zip(preds, truevals)) do (pred, trueval)
        Λvec, Q = eigen(cov(pred))
        normalized_trueval = Q*(trueval - mean(pred)) ./ sqrt.(Λvec)
        # all normalized values have to be in a given distance threshold
        all(abs.(normalized_trueval) .< dist_thresholds'; dims=1)[:]
    end
  
    calibrationvals = mean(containments)
    (; pvals, calibrationvals)
end

function sharpness(pred::Normal)
    2*std(pred)
end

# Measures volume of 1std ellipsoid.
function sharpness(pred::AbstractMvNormal)
    Λvec, _Q = eigen(cov(pred))
    d = length(pred)
    # https://en.wikipedia.org/wiki/Ellipsoid#Standard_equation_2
    π^(d/2)/gamma(d/2+1) * prod(sqrt.(Λvec))
end

end # module MvNormalCalibration
