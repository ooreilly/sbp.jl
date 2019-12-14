using Test
using SparseArrays

import sbp
using sbp.Operator: read_operator, to_matrix
using sbp.Grid:grid_xp
using sbp.Test:test_first_derivative, test_quadrature, test_first_derivative_sbp
include("operators.jl")

st = sbp.Strand1994

m = 10 


B = spzeros(m, m)
B[1,1] = -1
B[end,end] = 1

@test_throws UndefVarError st.load_descrption(order=10)
desc_D21, desc_H21 = st.load_description(order=2)

@testset "D_21, H_21" begin
        x, h, D, H = st.build_operators(desc_D21, desc_H21, m)
        ml = size(desc_D21.left, 1)
        mr = size(desc_D21.right, 1)
        test_first_derivative(D, x, x, 2, 1, ml, mr)
        test_quadrature(H, x, 1)
        test_first_derivative_sbp(H, D, B)
end



