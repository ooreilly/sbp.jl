import sbp
using sbp.Grid:grid_xp, grid_xm
using Test

@testset "Grid properties" begin
        n = 10
        x, h = grid_xp(n)
        @test isapprox(x[2] - x[1], x[end] - x[end-1])
        @test isapprox(x[1], 0)
        @test isapprox(x[end], 1.0)

        x, h = grid_xm(n)
        @test isapprox(x[2] - x[1], x[end] - x[end-1])
        @test isapprox(x[1], 0)
        @test isapprox(x[end], 1.0)

end
