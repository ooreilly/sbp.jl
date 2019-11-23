using Test
import sbp
using sbp.Metrics: build_jacobian, build_contravariant_basis
using sbp.Staggered: build_all_operators_2d, build_covariant_basis
using sbp.Grid: grid_xp, grid_xm, grid_2d_x, grid_2d_y

include("../operators/oreilly_petersson_2019/operators.jl")
using . Operators: build_operators
nx=10
ny=10


pp, mm, pm, mp = build_all_operators_2d(build_operators, nx, ny)

x1, h = grid_xp(nx)
y1, h = grid_xp(ny)

x = grid_2d_x(x1, ny)
y = grid_2d_y(y1, ny)

fx = x
fy = y

a = build_covariant_basis(fx, fy, mm, pm)


@testset "Jacobian" begin
        Jans = ones(mm.nx * mm.ny)
        J = build_jacobian(a)
        @test isapprox(J, Jans)
end

@testset "Contravariant basis" begin
        J = build_jacobian(a)
        b = build_contravariant_basis(J, a)

        one = ones(mm.nx * mm.ny)
        zero = zeros(mm.nx * mm.ny)

        @test isapprox(a.x_r1 .* b.r1_x + a.y_r1 .* b.r1_y, one)
        @test isapprox(a.x_r1 .* b.r2_x + a.y_r1 .* b.r2_y, zero, atol=1.0)
        @test isapprox(a.x_r2 .* b.r2_x + a.y_r2 .* b.r2_y, one)
        @test isapprox(a.x_r2 .* b.r1_x + a.y_r2 .* b.r1_y, zero, atol=1.0)

end

