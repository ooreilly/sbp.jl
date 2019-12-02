using Test
import sbp
using sbp.Staggered: boundary_matrix_p, boundary_matrix_m, build_operators_2d,
                     build_all_operators_2d, build_covariant_basis_mm, 
                     build_covariant_basis_pm


include("../operators/oreilly_petersson_2019/operators.jl")
using . Operators: build_operators

nx = 10
ny = 20
xp, xm, h, Dxp, Dxm, Hxp, Hxm, Pxp, Pxm = build_operators(nx)
yp, ym, h, Dyp, Dym, Hyp, Hym, Pyp, Pym = build_operators(ny)

@testset "Boundary matrix"  begin
        Bp = boundary_matrix_p(nx)
        @test size(Bp, 1) == nx
        @test size(Bp, 2) == nx + 1
        @test isapprox(sum(Bp), 0)

        Bp = boundary_matrix_p(nx, false)
        @test size(Bp, 1) == nx
        @test size(Bp, 2) == nx + 1
        @test isapprox(sum(Bp), 0)

        Bm = boundary_matrix_m(nx+1)
        @test size(Bm, 1) == nx + 1
        @test size(Bm, 2) == nx
        @test isapprox(sum(Bm), 0)
end

@testset "2D operators" begin
        Bxp = boundary_matrix_p(nx)
        Bxm = boundary_matrix_m(nx+1)
        Byp = boundary_matrix_p(ny)
        Bym = boundary_matrix_m(ny+1)
        pp = build_operators_2d(Dxp, Dyp, Pxp, Pyp, Hxp, Hyp, Bxp, Byp)
        @test pp.nx == nx
        @test pp.ny == ny
        @test size(pp.Dx, 1) == nx * ny
        @test size(pp.Dx, 2) == (nx + 1) * ny
        @test size(pp.Dy, 1) == nx * ny
        @test size(pp.Dy, 2) == nx * (ny + 1)

        pp, mm, pm, mp = build_all_operators_2d(build_operators, nx, ny)
        @test pp.nx == nx
        @test pp.ny == ny

        @test mm.nx == nx + 1
        @test mm.ny == ny + 1

        @test pm.nx == nx
        @test pm.ny == ny + 1

        @test mp.nx == nx + 1
        @test mp.ny == ny 

        @test size(pp.Dx) == size(pp.Px)
        @test size(pp.Dy) == size(pp.Py)
        @test size(pp.Bx) == size(pp.Bx)
        @test size(pp.By) == size(pp.By)

        @test size(mm.Dx) == size(mm.Px)
        @test size(mm.Dy) == size(mm.Py)
        @test size(mm.Bx) == size(mm.Bx)
        @test size(mm.By) == size(mm.By)

        @test size(pm.Dx) == size(pm.Px)
        @test size(pm.Dy) == size(pm.Py)
        @test size(pm.Bx) == size(pm.Bx)
        @test size(pm.By) == size(pm.By)
        
        @test size(mp.Dx) == size(mp.Px)
        @test size(mp.Dy) == size(mp.Py)
        @test size(mp.Bx) == size(mp.Bx)
        @test size(mp.By) == size(mp.By)

end

@testset "Covariant basis (cell-centers)" begin

        pp, mm, pm, mp = build_all_operators_2d(build_operators, nx, ny)
        fx = ones(pp.nx*pp.ny)
        fy = ones(pp.nx*pp.ny)
        a = build_covariant_basis_mm(fx, fy, mm, pm)

        x_r1 = zeros(mm.nx * mm.ny) 
        x_r2 = zeros(mm.nx * mm.ny) 
        y_r1 = zeros(mm.nx * mm.ny) 
        y_r2 = zeros(mm.nx * mm.ny) 

        @test size(a.x_r1) == size(a.x_r2)
        @test size(a.x_r1) == size(a.y_r1)
        @test size(a.x_r1) == size(a.y_r2)

        @test isapprox(a.x_r1, x_r1, atol=1.0)
        @test isapprox(a.x_r2, x_r2, atol=1.0)
        @test isapprox(a.y_r1, y_r1, atol=1.0)
        @test isapprox(a.y_r2, y_r2, atol=1.0)
end

@testset "Covariant basis (v1)" begin

        pp, mm, pm, mp = build_all_operators_2d(build_operators, nx, ny)
        fx = ones(pp.nx * pp.ny)
        fy = ones(pp.nx * pp.ny)
        a = build_covariant_basis_pm(fx, fy, mm, pm, mp)

        x_r1 = zeros(pm.nx * pm.ny) 
        x_r2 = zeros(pm.nx * pm.ny) 
        y_r1 = zeros(pm.nx * pm.ny) 
        y_r2 = zeros(pm.nx * pm.ny) 

        @test size(a.x_r1) == size(a.x_r2)
        @test size(a.x_r1) == size(a.y_r1)
        @test size(a.x_r1) == size(a.y_r2)

        @test isapprox(a.x_r1, x_r1, atol=1.0)
        @test isapprox(a.x_r2, x_r2, atol=1.0)
        @test isapprox(a.y_r1, y_r1, atol=1.0)
        @test isapprox(a.y_r2, y_r2, atol=1.0)
end

