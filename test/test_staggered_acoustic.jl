using Test
using SparseArrays
using LinearAlgebra
import sbp
using sbp.Sparse: diag_block_matrix_2x2
using sbp.Grid: grid_xp, grid_xm, grid_2d_x, grid_2d_y

st = sbp.StaggeredAcoustic
sop = sbp.OP2019


nx=10
ny=20

function mapping(nx, ny)
        x1, h = grid_xp(nx)
        y1, h = grid_xp(ny)
        
        p = grid_2d_x(x1, ny)
        q = grid_2d_y(y1, nx)
        x = p
        y = (1 .+ 0.1 * sin.(p)).*q
        return x, y
end

fx, fy = mapping(nx, ny)
desc = sop.load_description()
builder = n -> sop.build_operators(desc, n)
ops = st.init_operators(nx, ny, builder, fx, fy)

HJp = st.pressure_norm(ops)
HJ = st.velocity_norm(ops)
Gp = st.contravariant_metric_tensor(ops)
Ek = HJ*Gp
Ap = st.divergence(ops)
Av = st.gradient(ops, Gp)
S = st.pressure_sat(ops, Gp)
H = st.energy_norm(HJp, HJ, Gp)
A = st.spatial_discretization(Ap, Av - S)
d = st.pressure_source(ops, 0.5, 0.5, 4, 4) 


@testset "Grids" begin
        xp, yp = st.grids("p", nx, ny)
        @test size(xp,1) == (nx + 1) * (ny + 1)
        @test size(yp,1) == (nx + 1) * (ny + 1)
        x1, y1 = st.grids("v1", nx, ny)
        @test size(x1,1) == (nx + 0) * (ny + 1)
        @test size(y1,1) == (nx + 0) * (ny + 1)
        x2, y2 = st.grids("v2", nx, ny)
        @test size(x2,1) == (nx + 1) * (ny + 0)
        @test size(y2,1) == (nx + 1) * (ny + 0)
        xn, yn = st.grids("node", nx, ny)
        @test size(xn,1) == (nx + 0) * (ny + 0)
        @test size(yn,1) == (nx + 0) * (ny + 0)
end

@testset "Kinetic energy tensor" begin
        # Check symmetry
        @test isapprox(Ek, (Ek)')
        # Check positivity
        lam = real(eigvals(Matrix(Ek)))
        @test minimum(lam) > 0

end

@testset "Stability" begin
        # Check that the discretization is energy conservative to rounding error
        lam = real(eigvals(Matrix(A)))
        @test minimum(lam) > -1e-13
        @test maximum(lam) < 1e-13

        E = H*A + (H*A)'
        lam = real(eigvals(Matrix(E)))
        @test minimum(lam) > -1e-13
        @test maximum(lam) < 1e-13

end


