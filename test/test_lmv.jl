using Test
using SparseArrays
import sbp

m = 2 << 5
n = 2 << 4


@testset "read/write vector" begin
        file = "x.bin"

        for type in [Int32, Float32, Float64]
                x = ones(type, n)
                sbp.lmv.write_vector(file, x, verbose=true)
                sbp.lmv.write_vector(file, x, verbose=false)
                y = sbp.lmv.read_vector(type, file)
                @test typeof(y) == Array{type, 1}
                @test isapprox(x, y)
        end
        
        rm(file)
end

@testset "read/write CSR matrix" begin
        file = "A.bin"
        A = spzeros(m, n)
        A[1,1] = 1.0
        A[1,2] = 2.0
        A[1,3] = 3.0
        A[2,1] = 4.0
        A[2,2] = 5.0
        A[2,3] = 6.0
        sbp.lmv.write_csr_matrix(file, A, verbose=true)
        sbp.lmv.write_csr_matrix(file, A, verbose=false)
        B = sbp.lmv.read_csr_matrix(Float64, Int64, file, verbose=true)
        B = sbp.lmv.read_csr_matrix(Float64, Int64, file, verbose=false)
        @test typeof(B) == SparseMatrixCSC{Float64, Int64}
        @test isapprox(A, B)
        rm(file)
end

@testset "read/write dense matrix" begin
        file = "A.bin"
        A = zeros(n, n)
        A[1,1] = 1.0
        A[1,2] = 2.0
        A[1,3] = 3.0
        A[2,1] = 4.0
        A[2,2] = 5.0
        A[2,3] = 6.0
        sbp.lmv.write_dense_matrix(file, A, verbose=true)
        sbp.lmv.write_dense_matrix(file, A, verbose=false)
        B = sbp.lmv.read_dense_matrix(file, verbose=true)
        B = sbp.lmv.read_dense_matrix(file, verbose=false)
        @test typeof(B) == Array{Float64, 2}
        @test isapprox(A, B)

        rm(file)
end

@testset "read/write config" begin
        file = "config.txt"
        cfg = Dict("nt" => Int64(1), "dt" => Float64(0.1))
        sbp.lmv.write_config(file, cfg, verbose=false, skip_types=true)
        sbp.lmv.write_config(file, cfg, verbose=false)
        sbp.lmv.write_config(file, cfg, verbose=true)
        cfg2 = sbp.lmv.read_config(file, verbose=false)
        cfg2 = sbp.lmv.read_config(file, verbose=true)
        @test cfg == cfg2

        rm(file)


end
