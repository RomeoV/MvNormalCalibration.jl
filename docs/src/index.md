# MvNormalCalibration.jl

This package implements computing *calibration* and a measure of *sharpness*,
specifically for uni- and multivariate normal distributions, i.e. "Gaussians".
This page will present a brief usage example and rationale for the 
[`computecalibration`](@ref) and [`sharpness`](@ref) functions.
We also provide some [Mathematical Background](@ref) on calibration and an [API Reference](@ref).

## A Brief Usage Example
Before we go into the mathematical background, here's a quick usage example for measuring
calibration using [`computecalibration`](@ref).
In this example, we make a constant 2D prediction; however, the true distribution slightly 
changes at every time step.

```@example
using Distributions
using MvNormalCalibration: computecalibration, sharpness
using LinearAlgebra: I, hermitianpart

nobs = 500

# We predict that every observation is a zero mean, iso-normal...
predictions = [MvNormal(zeros(2), I(2)) for _ in 1:nobs]
# ... but actually, each observation has nonzero mean and small covariate terms.
# Note: Usually we don't know the true distributions.
true_distributions = [MvNormal((rand(2).-0.5)./10, I(2) + hermitianpart(rand(2,2)))
                      for _ in 1:nobs]
observations = rand.(true_distributions)

# Let's compute calibration for our predictions.
(; pvals, calibrationvals) = computecalibration(predictions, observations)

# As expected, there is a significant gap between theoretical and actual calibration values.
@info """
Maximum calibration discrepancy for (somewhat inaccurate) predictions
$(maximum(abs, pvals - calibrationvals) |> x->round(x; sigdigits=2))
"""

# However, if we had managed to predict each distribution accurately, things would look better.
(; pvals, calibrationvals) = computecalibration(true_distributions, observations)
@info """
and for true distributions
$(maximum(abs, pvals - calibrationvals) |> x->round(x; sigdigits=2)).
"""
```

### Plotting calibration results
It can also be useful to plot the resulting calibration values again the theoretical values.
Indeed, for perfect calibration, we expect the theoretical coverage probabilities (`pvals`)
to exactly match the observed coverage probabilities (`calibrationvals`) computed by [`computecalibration`](@ref).
Let's compute and plot calibration for overly optimistic predictions (variance too small,
assumes zero-mean), for "perfect" predictions, and for overly pessimistic predictions.

```@example
using Distributions, MvNormalCalibration, LinearAlgebra, Plots
plot_kwargs = (; label="Predictions", aspect_ratio=:equal, xlabel="Theoretical coverage",
                 ylabel="Observed coverage", xlims=(0,1), ylims=(0,1), size=(750, 280))
nobs = 500
optimistic_predictions  = [MvNormal(zeros(2), I(2)) for _ in 1:nobs]
true_distributions      = [MvNormal(rand(2)./10, I(2) + hermitianpart(rand(2,2)))
                           for _ in 1:nobs]
pessimistic_predictions = [1.3^2 * I(2)] .* optimistic_predictions
observations = rand.(true_distributions)

p1 = plot(computecalibration(optimistic_predictions, observations)...;
          title="Overconfident", plot_kwargs...)
p2 = plot(computecalibration(true_distributions, observations)...;
          title="Well calibrated", plot_kwargs...)
p3 = plot(computecalibration(pessimistic_predictions, observations)...;
          title="Pessimistic", plot_kwargs...)
for p in [p1, p2, p3]; plot!(p, [0,1], [0,1]; label="Ground Truth"); end
plot(p1, p2, p3; layout=(1,3))
```

## Measuring Sharpness
Calibration paints only half the picture. We can have well calibrated, but still imprecise
predictions if they have (unnecessarily) large uncertainty, i.e. they are *not sharp*.
Here is an example of such a setup, and how to compute and interpret [`sharpness`](@ref).

```@example
using Distributions, MvNormalCalibration, LinearAlgebra
preds_poor  = [MvNormal(zeros(2), I(2)) for _ in 1:1_000]
preds_sharp = [MvNormal( rand(D), (0.1^2).*I(2)) for D in preds_poor]
observations = rand.(preds_sharp)
@info """Average sharpness of poor predictions:
$(mean(sharpness.(preds_poor)) |> x->round(x; sigdigits=2))
"""
@info """Average sharpness of sharp predictions:
$(mean(sharpness.(preds_sharp)) |> x->round(x; sigdigits=2))
"""
```

Indeed we see that the one set of predictions has one hundred times lower uncertainty than
the other. But do they both predict the observations equally well? Let's compare both
predictions' calibration and sharpness. We will see that both predictions have almost 
identical calibration curves, despite their difference in uncertainty.

### Comparing predictions via calibration and sharpness

```@example
using Distributions, MvNormalCalibration, Plots, LinearAlgebra
import Random: seed!
preds_poor  = [MvNormal(zeros(2), I(2)) for _ in 1:1_000]
preds_sharp = [MvNormal( rand(D), (0.1^2).*I(2)) for D in preds_poor]

observations = rand.(preds_sharp)
pvals = 0:0.01:1
calibrationvals_poor  = computecalibration(preds_poor, observations; pvals)[:calibrationvals]
calibrationvals_sharp = computecalibration(preds_sharp, observations; pvals)[:calibrationvals]

plt1 = plot(pvals, calibrationvals_poor;
            label="Poor (constant) predictions", title="Calibration curve",
            aspect_ratio=:equal, xlabel="Theoretical coverage", ylabel="Observed coverage",
            xlims=(0,1), ylims=(0,1), size=(700, 360))
plot!(pvals, calibrationvals_sharp; label="Sharp predictions")

seed!(1)
xs = LinRange(-3, 3, 100)
plt2 = plot(xs, x->pdf(Normal(), x);
            label="Poor (constant) predictions", title="Predictions", xlims=(-3, 3),ylims=(0, 6), aspect_ratio=:equal, xlabel=raw"$x$", ylabel="Probability "*raw"$pdf(x)$")
for i in 1:6
   D = Normal(randn(), 0.1)
   xs_=LinRange(mean(D)-5*std(D), mean(D)+5*std(D), 100)
   plot!(plt2, xs_, x->pdf(D, x); color=2, label=(i==1 ? "Sharp predictions (i=1,2,..)" : nothing))
end
plot(plt1, plt2)
```

We can see that both sets of predictions have practically identical calibration curves, but 
extremely different sharpness -- i.e. the sharp predictions are much more precise for
each sample.
