using Test
import sbp
using sbp.Sparse: block_matrix, block_matrix_insert, block_matrix_2x2,
                  diag_block_matrix_2x2

m = 4
rows = [m, m]
cols = [m, m]
A = ones(m,m)

@testset "block_matrix_insert" begin
        B = block_matrix(rows, cols)
        B = block_matrix_insert(B, rows, cols, 1, 1, A)
        @test isapprox(A, B[1:m,1:m])
end

@testset "2x2_block_matrix" begin
        B = block_matrix_2x2([1],[2],[3],[4])
        @test B[1,1] == 1
        @test B[1,2] == 2
        @test B[2,1] == 3
        @test B[2,2] == 4
end

@testset "2x2_block_matrix_diag" begin
        B = diag_block_matrix_2x2([1],[2])
        @test B[1,1] == 1
        @test B[1,2] == 0
        @test B[2,1] == 0
        @test B[2,2] == 2
end
        

