using Test

import sbp
using sbp.Operator: read_operator, to_matrix
using sbp.Grid:grid_xp, grid_xm
using sbp.Test:test_first_derivative, test_quadrature, test_interpolation,
        test_first_derivative_sbp, test_interpolation_sbp
include("operators.jl")

using . Operators: build_operators, ODp, ODm, OPp, OPm

m = 10 

xp, xm, h, Dp, Dm, Hp, Hm, Pp, Pm = build_operators(m)

@time begin

@testset "Dp_42" begin
        ml = size(ODp.left, 1)
        mr = size(ODp.right, 1)
        test_first_derivative(Dp, xp, xm, 4, 2, ml, mr)
end

@testset "Dm_42" begin
        ml = size(ODm.left, 1)
        mr = size(ODm.right, 1)
        test_first_derivative(Dm, xm, xp, 4, 2, ml, mr)
end

@testset "Pp_42" begin
        ml = size(OPp.left, 1)
        mr = size(OPp.right, 1)
        test_interpolation(Pp, xp, xm, 3, 1, ml, mr)
end

@testset "Pm_42" begin
        ml = size(OPm.left, 1)
        mr = size(OPm.right, 1)
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


