using Test
using LinearAlgebra
import sbp
using sbp.Config: Tv
using sbp.Metrics: build_jacobian, build_contravariant_basis, 
                   build_covariant_metric_tensor, 
                   build_contravariant_metric_tensor
using sbp.Staggered: build_all_operators_2d, build_covariant_basis_mm
using sbp.Grid: grid_xp, grid_xm, grid_2d_x, grid_2d_y

sop = sbp.OP2019
nx=10
ny=20


pp, mm, pm, mp = build_all_operators_2d(sop.build_operators, nx, ny)

x1, h = grid_xp(nx)
y1, h = grid_xp(ny)

x = grid_2d_x(x1, ny)
y = grid_2d_y(y1, nx)

fx = x
fy = y

a = build_covariant_basis_mm(fx, fy, mm, pm)
Ga = build_covariant_metric_tensor(a)


@testset "Jacobian" begin
        Jans = ones(mm.nx * mm.ny)
        J = build_jacobian(a)
        @test isapprox(J, Jans)
end

@testset "Covariant Metric tensor" begin
         Ga = build_covariant_metric_tensor(a)
         # Check symmetry
         @test isapprox(Ga.a12, Ga.a21)
        
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

        @testset "Contravariant Metric tensor" begin
                Gb = build_contravariant_metric_tensor(b)
                # Check symmetry
                @test isapprox(Gb.b12, Gb.b21)
                # Check orthogonality
                @test isapprox(Ga.a11 .* Gb.b11 + Ga.a12 .* Gb.b12, one)
                @test isapprox(Ga.a11 .* Gb.b12 + Ga.a12 .* Gb.b22, zero, 
                               atol=1.0)
                @test isapprox(Ga.a12 .* Gb.b11 + Ga.a22 .* Gb.b12, zero, 
                               atol=1.0)
                @test isapprox(Ga.a21 .* Gb.b12 + Ga.a22 .* Gb.b22, one)
        
        end

end



