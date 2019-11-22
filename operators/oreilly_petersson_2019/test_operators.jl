using Test

import sbp
using sbp.Operator: read_operator, to_matrix
using sbp.Grid:grid_xp, grid_xm
using sbp.Test:test_first_derivative, test_quadrature

path = string(@__DIR__, "/")

m = 10 
n = 10 
xm, h = grid_xm(m+1)
xp, h = grid_xp(m)

@time begin

@testset "Dp_42" begin
        op = read_operator(string(path,  "Dp_42.txt"))
        Dp = to_matrix(op, m, m + 1) / h
        ml = size(op.left, 1)
        mr = size(op.right, 1)
        test_first_derivative(Dp, xp, xm, 4, 2, ml, mr)
end

@testset "Dm_42" begin
        op = read_operator(string(path,  "Dm_42.txt"))
        Dm = to_matrix(op, m + 1, m) / h
        ml = size(op.left, 1)
        mr = size(op.right, 1)
        test_first_derivative(Dm, xm, xp, 4, 2, ml, mr)
end

@testset "Hp_42" begin
        op = read_operator(string(path,  "Hp_42.txt"))
        Hp = to_matrix(op, m, m) * h
        test_quadrature(Hp, xp, 3)

end

end


