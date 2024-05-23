using Test
using MvNormalCalibration
using Distributions
using LinearAlgebra: Diagonal, I

@testset "sampling from ground truths" begin
  # When the ground truths are sampled exactly from the predictions,
  # we will always be calibrated, no matter the predictions.
  samples = map(1:10_000) do _
    pred = MvNormal(rand(2), Diagonal(rand(2)))
    measurement = rand(pred)
    (; pred, measurement)
  end
  preds = getfield.(samples, :pred)
  measurements = getfield.(samples, :measurement)
  (; pvals, calibrationvals) = computecalibration(preds, measurements; pvals=0.01:0.01:0.99)
  @test pvals ≈ calibrationvals rtol=1e-2
end

@testset "sharpness values" begin
  # we know the volume of a circle in 2D / 3D
  r = 1.5
  @test sharpness(MvNormal(ones(2), r^2*I(2))) ≈ π*r^2
  @test sharpness(MvNormal(ones(3), r^2*I(3))) ≈ (4/3*π*r^3)
end
