@testset "Code quality (Aqua.jl)" begin
    # gotta split this: see https://github.com/JuliaTesting/Aqua.jl/issues/77
    Aqua.test_all(MvNormalCalibration, ambiguities = false)
    Aqua.test_ambiguities(MvNormalCalibration)
end

@testset "Code linting (JET.jl)" begin
    JET.test_package(MvNormalCalibration;
        target_defined_modules = true)
end
