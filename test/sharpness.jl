@testset "sharpness values" begin
    @testset "univariate" begin
        for _ in 1:10
            D = Normal(rand(), rand() + 0.5)
            @test sharpness(D) ≈ 2 * std(D)
        end
    end
    @testset "multivariate" begin
        # we know the volume of a circle in 2D / 3D
        r = 1.5
        @test sharpness(MvNormal(rand(1), r^2 * I(1))) ≈ 2 * r^1
        @test sharpness(MvNormal(rand(2), r^2 * I(2))) ≈ π * r^2
        @test sharpness(MvNormal(rand(2), Diagonal([3^2, 2^2]))) ≈ π * 3 * 2
        @test sharpness(MvNormal(rand(3), r^2 * I(3))) ≈ (4 / 3 * π * r^3)
    end
end
