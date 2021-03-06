using Test
import sbp


st = sbp.Staggered
sop = sbp.OP2019

nx = 10
ny = 20

desc = sop.load_description()
xp, xm, h, Dxp, Dxm, Hxp, Hxm, Pxp, Pxm = sop.build_operators(desc, nx)
yp, ym, h, Dyp, Dym, Hyp, Hym, Pyp, Pym = sop.build_operators(desc, ny)

@testset "Boundary matrix"  begin
        Bp = st.boundary_matrix_p(nx)
        @test size(Bp, 1) == nx
        @test size(Bp, 2) == nx + 1
        @test isapprox(sum(Bp), 0)

        Bp = st.boundary_matrix_p(nx, false)
        @test size(Bp, 1) == nx
        @test size(Bp, 2) == nx + 1
        @test isapprox(sum(Bp), 0)

        Bm = st.boundary_matrix_m(nx+1)
        @test size(Bm, 1) == nx + 1
        @test size(Bm, 2) == nx
        @test isapprox(sum(Bm), 0)
end

@testset "2D operators" begin
        Bxp = st.boundary_matrix_p(nx)
        Bxm = st.boundary_matrix_m(nx+1)
        Byp = st.boundary_matrix_p(ny)
        Bym = st.boundary_matrix_m(ny+1)
        pp = st.build_operators_2d(Dxp, Dyp, Pxp, Pyp, Hxp, Hyp, Bxp, Byp)
        @test pp.nx == nx
        @test pp.ny == ny
        @test size(pp.Dx, 1) == nx * ny
        @test size(pp.Dx, 2) == (nx + 1) * ny
        @test size(pp.Dy, 1) == nx * ny
        @test size(pp.Dy, 2) == nx * (ny + 1)

        builder = n -> sop.build_operators(desc, n)
        pp, mm, pm, mp = st.build_all_operators_2d(builder, nx, ny)
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

        builder = n -> sop.build_operators(desc, n)
        pp, mm, pm, mp = st.build_all_operators_2d(builder, nx, ny)
        fx = ones(pp.nx*pp.ny)
        fy = ones(pp.nx*pp.ny)
        a = st.build_covariant_basis_mm(fx, fy, mm, pm)

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
                                                          
        builder = n -> sop.build_operators(desc, n)
        pp, mm, pm, mp = st.build_all_operators_2d(builder, nx, ny)
        fx = ones(pp.nx * pp.ny)
        fy = ones(pp.nx * pp.ny)
        a = st.build_covariant_basis_pm(fx, fy, mm, pm, mp)

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

@testset "Covariant basis (v2)" begin

        builder = n -> sop.build_operators(desc, n)
        pp, mm, pm, mp = st.build_all_operators_2d(builder, nx, ny)
        fx = ones(pp.nx * pp.ny)
        fy = ones(pp.nx * pp.ny)
        a = st.build_covariant_basis_mp(fx, fy, mm, pm, mp)

        x_r1 = zeros(mp.nx * mp.ny) 
        x_r2 = zeros(mp.nx * mp.ny) 
        y_r1 = zeros(mp.nx * mp.ny) 
        y_r2 = zeros(mp.nx * mp.ny) 

        @test size(a.x_r1) == size(a.x_r2)
        @test size(a.x_r1) == size(a.y_r1)
        @test size(a.x_r1) == size(a.y_r2)

        @test isapprox(a.x_r1, x_r1, atol=1.0)
        @test isapprox(a.x_r2, x_r2, atol=1.0)
        @test isapprox(a.y_r1, y_r1, atol=1.0)
        @test isapprox(a.y_r2, y_r2, atol=1.0)
end

