import sbp
err = sbp.MMS
using Test
using LinearAlgebra
using DataStructures

@testset "extract_vars" begin
        n = [4, 8, 12]
        data = Array(1:sum(n))
        vars = OrderedDict("u" => n[1], "v" => n[2], "w" => n[3])
        ans_u = Array(1:n[1])
        ans_v = Array(n[1]+1:n[1] + n[2])
        ans_w = Array(n[1]+n[2]+1:sum(n))
        expected = Dict("u" => ans_u,
                         "v" => ans_v,
                         "w" => ans_w)
        extracted_vars = err.extract_vars(data, vars)
        @test extracted_vars == expected
end

@testset "l2_errors" begin
        u = ones(4)
        v = ones(4)
        H = Matrix(I, 4, 4)
        solution = Dict("u" => u, "v" => v)
        norms = Dict("u" => H, "v" => H)
        mms = Dict("u" => u, "v" => v)
        errors = err.l2_errors(solution, mms, norms)
        expected = Dict("u" => 0.0, "v" => 0.0)
        @test errors == expected
end

@testset "functional_errors" begin
        u = ones(4)
        v = ones(4)
        H = Matrix(I, 4, 4)
        solution = Dict("u" => u, "v" => v)
        norms = Dict("u" => H, "v" => H)
        mms = Dict("u" => 4.0, "v" => 4.0)
        errors = err.functional_errors(solution, mms, norms)
        expected = Dict("u" => 0.0, "v" => 0.0)
        @test errors == expected
end

@testset "log2_convergence_rates" begin
        u = [16, 8, 4, 2]
        q = err.log2_convergence_rates(u)
        ans = [Inf, 1, 1, 1]
        @test q == ans
end
