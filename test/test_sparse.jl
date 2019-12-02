using Test
import sbp
using sbp.Sparse: block_matrix, block_matrix_insert

m = 4
rows = [m, m]
cols = [m, m]
A = ones(m,m)

@testset "insert" begin
        B = block_matrix(rows, cols)
        B = block_matrix_insert(B, rows, cols, 1, 1, A)
        @test isapprox(A, B[1:m,1:m])
end

