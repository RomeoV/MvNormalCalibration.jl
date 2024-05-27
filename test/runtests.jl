using Test
using MvNormalCalibration
using Distributions
using LinearAlgebra: Diagonal, I, isposdef, hermitianpart!
using StructArrays
using Aqua
using JET

@testset "MvNormalCalibration.jl" begin
@testset "Code quality (Aqua.jl)" begin
    # gotta split this: see https://github.com/JuliaTesting/Aqua.jl/issues/77
    Aqua.test_all(MvNormalCalibration, ambiguities = false)
    Aqua.test_ambiguities(MvNormalCalibration)
end

@testset "Code linting (JET.jl)" begin
    JET.test_package(MvNormalCalibration;
        target_defined_modules = true)
end

@testset "sampling from ground truths" begin
    @testset "univariate" begin
        # When the ground truths are sampled exactly from the predictions,
        # we will always be calibrated, no matter the predictions.
        (; preds, measurements) = map(1:500_000) do _
            pred = Normal(rand(), rand()+1.0)
            measurement = rand(pred)
            (; preds=pred, measurements=measurement)
        end |> StructArray
        (; pvals, calibrationvals) = computecalibration(preds, measurements;
                                                        pvals=0.00:0.05:1.0)
        @test pvals ≈ calibrationvals rtol=1e-2
    end

    @testset "compare univariate and multivariate " begin
        (; preds, measurements) = map(1:100) do _
            pred = Normal(rand(), rand()+1)
            measurement = rand(pred)
            (; preds=pred, measurements=measurement)
        end |> StructArray
        preds_mv = map(preds) do pred
            MvNormal([mean(pred)], Diagonal([std(pred)^2]))
        end
        measurements_mv = map(measurements) do measurement
            [measurement]
        end
        (; pvals, calibrationvals) = computecalibration(preds, measurements)
        (; pvals_mv, calibrationvals_mv) = let pvals, calibrationvals
            (; pvals, calibrationvals) = computecalibration(preds_mv, measurements_mv)
            (; pvals_mv=pvals, calibrationvals_mv = calibrationvals)
        end
        @test pvals == pvals_mv
        @test calibrationvals ≈ calibrationvals_mv
    end

    @testset "multivariate" begin
        @testset for d in 1:5
            # When the ground truths are sampled exactly from the predictions,
            # we will always be calibrated, no matter the predictions.
            (; preds, measurements) = map(1:100_000) do _
                Σ = zeros(d,d)
                while !isposdef(Σ); Σ .= 3*I(d) + hermitianpart!(rand(d,d)); end
                pred = MvNormal(rand(d), Σ)
                measurement = rand(pred)
                (; preds=pred, measurements=measurement)
            end |> StructArray
            (; pvals, calibrationvals) = computecalibration(preds, measurements;
                                                            pvals=0.0:0.02:1.0)
            @test pvals ≈ calibrationvals rtol=1e-1
        end
    end
end

@testset "sharpness values" begin
    @testset "univariate" begin
        for _ in 1:10
            D = Normal(rand(), rand() + 0.5)
            @test sharpness(D) ≈ 2*std(D)
        end
    end
    @testset "multivariate" begin
        # we know the volume of a circle in 2D / 3D
        r = 1.5
        @test sharpness(MvNormal(rand(1), r^2*I(1))) ≈ 2*r^1
        @test sharpness(MvNormal(rand(2), r^2*I(2))) ≈ π*r^2
        @test sharpness(MvNormal(rand(2), Diagonal([3^2, 2^2]))) ≈ π*3*2
        @test sharpness(MvNormal(rand(3), r^2*I(3))) ≈ (4/3*π*r^3)
    end
end

end
