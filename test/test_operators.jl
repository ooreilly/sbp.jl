import sbp
using sbp.Operator: read_operator, to_matrix
using sbp.Grid:grid_xp, grid_xm
using sbp.Test:test_first_derivative
using Test

op_path = ["../operators/ooreilly_petersson_2019/"]

m = 10 
n = 10 
xm, h = grid_xm(m+1)
x, h = grid_xp(m)
xp, h = grid_xp(m)

@time begin

for path in op_path
@testset "D1_12" begin
        op = read_operator(string(path,  "D1_21.txt"))
        D = to_matrix(op, m, n) / h
        ml = size(op.left, 1)
        mr = size(op.right, 1)
        test_first_derivative(D, xp, xp, 2, 1, ml, mr)
end

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

@testset "Hm_42" begin
        op = read_operator(string(path,  "Hm_42.txt"))
        Hm = to_matrix(op, m + 1, m + 1) * h
        test_quadrature(Hm, xm, 3)

        #flipped = reverse(reverse(op.left, dims=1), dims=2)
        #show(IOContext(stdout, :compact=>false), "text/plain", flipped)
        #display(op.right)
        #println()
        

end

end

end # time

