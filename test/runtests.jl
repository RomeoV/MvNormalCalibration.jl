using Test
using MvNormalCalibration
using Distributions
using LinearAlgebra: Diagonal, I, isposdef, hermitianpart!
using StructArrays
using Aqua
using JET

@testset "MvNormalCalibration.jl" begin
    include("test_prelude.jl")
    include("calibration.jl")
    include("sharpness.jl")
end
