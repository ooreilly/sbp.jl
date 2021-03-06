using Test

import sbp
using sbp.Operator: read_operator, to_matrix
using sbp.Grid:grid_xp, grid_xm
using sbp.Test:test_first_derivative, test_quadrature, test_interpolation,
        test_first_derivative_sbp, test_interpolation_sbp
include("operators.jl")

using . OP2019: build_operators, load_description

m = 10 

desc = load_description()
xp, xm, h, Dp, Dm, Hp, Hm, Pp, Pm = build_operators(desc, m)

@time begin

@testset "Dp_42" begin
        ml = size(desc.Dp.left, 1)
        mr = size(desc.Dp.right, 1)
        test_first_derivative(Dp, xp, xm, 4, 2, ml, mr)
end

@testset "Dm_42" begin
        ml = size(desc.Dm.left, 1)
        mr = size(desc.Dm.right, 1)
        test_first_derivative(Dm, xm, xp, 4, 2, ml, mr)
end

@testset "Pp_42" begin
        ml = size(desc.Pp.left, 1)
        mr = size(desc.Pp.right, 1)
        test_interpolation(Pp, xp, xm, 3, 1, ml, mr)
end

@testset "Pm_42" begin
        ml = size(desc.Pm.left, 1)
        mr = size(desc.Pm.right, 1)
        test_interpolation(Pm, xm, xp, 3, 1, ml, mr)
end

@testset "Hp_42" begin
        test_quadrature(Hp, xp, 3)

end

@testset "Hm_42" begin
        test_quadrature(Hm, xm, 3)

end

@testset "D_SBP" begin
        B = 0*Dp
        B[1,1] = -1
        B[end,end] = 1
        test_first_derivative_sbp(Hp, Dp, Hm, Dm, B)
end

@testset "P_SBP" begin
        test_interpolation_sbp(Hp, Pp, Hm, Pm)
end

end
