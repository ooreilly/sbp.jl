import sbp
using sbp.Operator: read_operator, write_operator, OperatorData, equals
using sbp.Operator: read_array, write_array, write_stencil, read_stencil
using sbp.Operator: version, parse_version, Version, check_version, to_matrix,
                    to_vector
using Test

@testset "Version" begin
        u = Version(1,2,3)
        v = parse_version("v1.2.3")
        @test u == v
        u = Version(version.major, version.minor, version.patch)
        check_version(u)
end

@testset "Array" begin
        m = 4
        n = 3
        a = zeros(m, n)
        f = open("array.txt", "w")
        write_array(f, "array", a)
        close(f)
        
        
        f = open("array.txt")
        label, b = read_array(f) 
        @test label == "array"
        @test a == b
        close(f)
        rm("array.txt")

end

@testset "Stencil" begin
        f = open("stencil.txt", "w")
        stencil = zeros(5)
        offset = 2
        write_stencil(f, "stencil", stencil, offset) 
        close(f)
        

        f = open("stencil.txt")
        label, st, off = read_stencil(f)
        @test label == "stencil"
        @test off == offset
        @test st == stencil
        close(f)
        rm("stencil.txt")

end

@testset "Operator" begin

        left = zeros(1,2)
        left[1,1] = -1.0
        left[end,end] = 1.0
        right = zeros(1,2)
        right[1,1] = -1.0
        right[end,end] = 1.0
        interior = zeros(3)
        interior[1] = -0.5
        interior[end] = 0.5
        offset = -2
        
        Op = OperatorData(version, left, right, interior, offset)
        write_operator("test.txt", Op)
        Op2 = read_operator("test.txt")
        @test equals(Op, Op2)
        rm("test.txt")

        @testset "to_matrix" begin
                @test_throws ErrorException("Too few rows") D = to_matrix(Op, 1, 1)
                m = 3
                n = 3
                D = to_matrix(Op, m, n, false)
                ans = zeros(m, n)
                ans[1,1] = -1.0
                ans[1,2] = 1.0
                ans[2,1] = -0.5
                ans[2,3] = 0.5
                ans[end,end-1] = -1.0
                ans[end,end] = 1.0
                @test isapprox(D, ans)
                Dsparse = to_matrix(Op, m, n, true)
                @test isapprox(Dsparse, ans)
        end

        @testset "to_vector" begin

                
                left = zeros(1, 2)
                interior = zeros(1)
                left[1] = -1.0
                left[2] = 1.0
                right = zeros(1, 2)
                right[1] = -1.0
                right[2] = 1.0
                interior[1] = 2.0
                offset = 0
                Op = OperatorData(version, left, right, interior, offset)
                @test_throws ErrorException("Too few elements") D = to_vector(Op, 3)
                m = 5
                ans = zeros(m)
                ans[1] = -1.0
                ans[2] = 1.0
                ans[3] = 2.0
                ans[end-1] = -1.0
                ans[end] = 1.0
                u = to_vector(Op, m)
                @test isapprox(u, ans)

        end

end
