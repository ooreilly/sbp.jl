using Test
import sbp
using sbp.Collocated: boundary_matrix, build_operators_2d, build_covariant_basis


st=sbp.Strand1994

order=2
nx = 10
ny = 20
desc_D, desc_H = st.load_description(order=order)
x, h, Dx, Hx = st.build_operators(desc_D, desc_H, nx)
y, h, Dy, Hy = st.build_operators(desc_D, desc_H, ny)

@testset "Boundary matrix"  begin
        Bp = boundary_matrix(nx)
        @test size(Bp, 1) == nx
        @test size(Bp, 2) == nx
        @test isapprox(sum(Bp), 0)

        Bp = boundary_matrix(nx, false)
        @test size(Bp, 1) == nx
        @test size(Bp, 2) == nx
        @test isapprox(sum(Bp), 0)
end

@testset "2D operators" begin
        Bx = boundary_matrix(nx)
        By = boundary_matrix(ny)
        op = build_operators_2d(Dx, Dy, Hx, Hy, Bx, By)
        @test op.nx == nx
        @test op.ny == ny
        @test size(op.Dx, 1) == nx * ny
        @test size(op.Dx, 2) == nx * ny
        @test size(op.Dy, 1) == nx * ny
        @test size(op.Dy, 2) == nx * ny
end

@testset "Covariant basis" begin
        Bx = boundary_matrix(nx)
        By = boundary_matrix(ny)
        op = build_operators_2d(Dx, Dy, Hx, Hy, Bx, By)
        fx = ones(op.nx*op.ny)
        fy = ones(op.nx*op.ny)
        a = build_covariant_basis(fx, fy, op)

        x_r1 = zeros(op.nx * op.ny) 
        x_r2 = zeros(op.nx * op.ny) 
        y_r1 = zeros(op.nx * op.ny) 
        y_r2 = zeros(op.nx * op.ny) 

        @test size(a.x_r1) == size(a.x_r2)
        @test size(a.x_r1) == size(a.y_r1)
        @test size(a.x_r1) == size(a.y_r2)

        @test isapprox(a.x_r1, x_r1, atol=1.0)
        @test isapprox(a.x_r2, x_r2, atol=1.0)
        @test isapprox(a.y_r1, y_r1, atol=1.0)
        @test isapprox(a.y_r2, y_r2, atol=1.0)
end

