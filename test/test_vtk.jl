using Test
import sbp
using sbp.Grid:grid_xp, grid_xm, grid_2d_x, grid_2d_y
using sbp.VTK:vtk_write, vtk_read


@testset "Structured Grid 2D data" begin
        nx = 20
        ny = 10
        x1, h = grid_xm(nx)
        y1, h = grid_xp(ny)
        x = grid_2d_x(x1, ny)
        y = grid_2d_y(y1, nx)
        z = x .* y

        vtk_write("test.vtk", x, y, z, nx, ny)
        xr, yr, zr = vtk_read("test.vtk")
        @test size(x, 1) == size(xr, 1)
        @test size(y, 1) == size(yr, 1)
        @test size(z, 1) == size(zr, 1)
        @test isapprox(x, xr, rtol=1e-6)
        @test isapprox(y, yr, rtol=1e-6)
        @test isapprox(z, zr, rtol=1e-6)
        rm("test.vtk")
end




