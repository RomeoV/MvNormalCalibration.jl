@testset "sampling from ground truths" begin
    @testset "univariate" begin
        @testset for scheme in [CentralPredictionSet, LeftPredictionSet]
            # When the ground truths are sampled exactly from the predictions,
            # we will always be calibrated, no matter the predictions.
            (; preds, measurements) = map(1:500_000) do _
                pred = Normal(rand(), rand() + 1.0)
                measurement = rand(pred)
                (; preds = pred, measurements = measurement)
            end |> StructArray
            (; pvals, calibrationvals) = computecalibration(scheme,
                preds, measurements;
                pvals = 0.00:0.05:1.0)
            @test pvals≈calibrationvals rtol=1e-2
        end
    end

    @testset "compare univariate and multivariate " begin
        (; preds, measurements) = map(1:100) do _
            pred = Normal(rand(), rand() + 1)
            measurement = rand(pred)
            (; preds = pred, measurements = measurement)
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
            (; pvals_mv = pvals, calibrationvals_mv = calibrationvals)
        end
        @test pvals == pvals_mv
        @test calibrationvals ≈ calibrationvals_mv
    end

    @testset "multivariate" begin
        @testset for d in 1:5
            # When the ground truths are sampled exactly from the predictions,
            # we will always be calibrated, no matter the predictions.
            (; preds, measurements) = map(1:100_000) do _
                Σ = zeros(d, d)
                while !isposdef(Σ)
                    Σ .= 3 * I(d) + hermitianpart!(rand(d, d))
                end
                pred = MvNormal(rand(d), Σ)
                measurement = rand(pred)
                (; preds = pred, measurements = measurement)
            end |> StructArray
            (; pvals, calibrationvals) = computecalibration(preds, measurements;
                pvals = 0.0:0.02:1.0)
            @test pvals≈calibrationvals rtol=1e-1
        end
    end
end
