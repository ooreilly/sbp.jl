import sbp
using sbp.Grid:grid_xp, grid_xm, grid_2d_x, grid_2d_y
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

@testset "2D grid vector" begin
        nx = 10
        ny = 12
        x1, h = grid_xp(nx)
        y1, h = grid_xp(ny)
        x = grid_2d_x(x1, ny)
        @test size(x) == (nx * ny,)
        @test isapprox(x[1:nx], x1[1] * ones(nx))
        
        y = grid_2d_y(y1, nx)
        @test size(y) == (nx * ny,)
        @test isapprox(y[1:ny], y1)

end
