using Test
using SparseArrays
using LinearAlgebra
import sbp
using sbp.StaggeredAcoustic: init_operators, pressure_norm, velocity_norm,
                             energy_norm, contravariant_metric_tensor, 
                             divergence, gradient,
                             spatial_discretization, pressure_sat, grids, 
                             pressure_source   
using sbp.Sparse: diag_block_matrix_2x2
using sbp.Grid: grid_xp, grid_xm, grid_2d_x, grid_2d_y

include("../operators/oreilly_petersson_2019/operators.jl")
using . Operators: build_operators

nx=10
ny=10

function mapping(nx, ny)
        x1, h = grid_xp(nx)
        y1, h = grid_xp(ny)
        
        x = grid_2d_x(x1, ny)
        y = grid_2d_y(y1, ny)
        return x, y
end


fx, fy = mapping(nx, ny)
ops = init_operators(nx, ny, build_operators, fx, fy)

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
        xp, yp = grids("p", nx, ny)
        @test size(xp,1) == (nx + 1) * (ny + 1)
        @test size(yp,1) == (nx + 1) * (ny + 1)
        x1, y1 = grids("v1", nx, ny)
        @test size(x1,1) == (nx + 0) * (ny + 1)
        @test size(y1,1) == (nx + 0) * (ny + 1)
        x2, y2 = grids("v2", nx, ny)
        @test size(x2,1) == (nx + 1) * (ny + 0)
        @test size(y2,1) == (nx + 1) * (ny + 0)
        xn, yn = grids("node", nx, ny)
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
        @test minimum(lam) > -1e13
        @test maximum(lam) < 1e13

        E = H*A + (H*A)'
        lam = real(eigvals(Matrix(E)))
        @test minimum(lam) > -1e13
        @test maximum(lam) < 1e13

end


