using Test
using SparseArrays
using LinearAlgebra
import sbp
using sbp.CollocatedAcoustic: init_operators, pressure_norm, velocity_norm,
                             energy_norm, contravariant_metric_tensor, 
                             divergence, gradient,
                             spatial_discretization, pressure_sat, grid,
                             pressure_source   
using sbp.Sparse: diag_block_matrix_2x2
using sbp.Grid: grid_xp, grid_xm, grid_2d_x, grid_2d_y

st=sbp.Strand1994

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
desc_D, desc_H = st.load_description(order=2)
ops = init_operators(nx, ny, st, desc_D, desc_H, fx, fy)

HJp = pressure_norm(ops)
HJ = velocity_norm(ops)
Gp = contravariant_metric_tensor(ops)
Ek = HJ*Gp
Ap = divergence(ops)
Av = gradient(ops, Gp)
S = pressure_sat(ops, Gp)
H = energy_norm(HJp, HJ, Gp)
A = spatial_discretization(Ap, Av - S)
d = pressure_source(ops, 0.5, 0.5, 4, 4) 


@testset "Grids" begin
        x, y = grid(nx, ny)
        @test size(x,1) == nx * ny
        @test size(y,1) == nx * ny
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
        @test minimum(lam) > -1e13
        @test maximum(lam) < 1e13

        E = H*A + (H*A)'
        lam = real(eigvals(Matrix(E)))
        @test minimum(lam) > -1e13
        @test maximum(lam) < 1e13

end


