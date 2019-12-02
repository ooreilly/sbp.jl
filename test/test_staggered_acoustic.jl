using Test
using SparseArrays
using LinearAlgebra
import sbp
using sbp.Staggered: build_all_operators_2d, 
                   build_covariant_basis_mm, 
                   build_covariant_basis_pm,
                   build_covariant_basis_mp
using sbp.StaggeredAcoustic: init_operators, pressure_norm, velocity_norm,
        contravariant_metric_tensor
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


ops = init_operators(nx, ny, build_operators, mapping)

HJp = pressure_norm(ops)
HJ = velocity_norm(ops)
Gp = contravariant_metric_tensor(ops)
Ek = HJ*Gp

@testset "Kinetic energy tensor" begin
        # Check symmetry
        @test isapprox(Ek, (Ek)')
        # Check positivity
        lam = real(eigvals(Matrix(Ek)))
        @test minimum(lam) > 0

end
