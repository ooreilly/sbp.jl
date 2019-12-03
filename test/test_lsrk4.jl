using Test
import sbp
using sbp.LSRK4: lsrk4

y = zeros(1)

@testset "Integrate cubic monomial" begin

        e = ones(1)

        function g(y::AbstractArray, t::Float64)
                return e * t ^ 3
        end
        
        nt = 10
        dt = 0.1
        Y = zeros(nt+1)
        u = zeros(1)
        dy = zeros(1)
        
        ans = zeros(nt+1)
        for k in 1:nt
                ans[k+1] = 0.25 * (k * dt) ^ 4 
        end
        
        tk=0.0
        for k in 1:nt
                tk = (k - 1) * dt
                @time global y = lsrk4(g, y, dy, tk, dt)
                Y[k+1] = y[1]
        end
        @test isapprox(ans, Y)
end

