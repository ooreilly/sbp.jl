using Test
import sbp
using sbp.StaggeredRK4: staggered_rk4

y = zeros(1)

@testset "Integrate cubic monomial" begin

        e = ones(1)

        function f(y::AbstractArray, t::Float64)
                return e * t ^ 3
        end

        function g(y::AbstractArray, t::Float64)
                return e * t ^ 3
        end
        
        nt = 10
        dt = 0.1
        u = zeros(1)
        v = zeros(1)
        U = zeros(nt+1)
        V = zeros(nt+1)
        
        ans_u = zeros(nt+1)
        ans_v = zeros(nt+1)
        for k in 1:(nt+1)
                ans_u[k] = 0.25 * ((k - 1) * dt) ^ 4 
                ans_v[k] = 0.25 * ((k - 0.5) * dt) ^ 4 
        end

        v[1] = ans_v[1]
        
        tk=0.0
        V[1] = v[1]
        for k in 1:nt
                tk = (k - 1) * dt
                @time u, v = staggered_rk4(f, g, u, v, tk, dt)
                U[k+1] = u[1]
                V[k+1] = v[1]
        end
        @test isapprox(ans_u, U)
        @test isapprox(ans_v, V)
end

