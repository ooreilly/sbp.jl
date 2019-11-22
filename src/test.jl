module Test

using sbp.Config: Tv
using Test

function taylor_poly(x::AbstractArray, p::Int64)
        y = x.^p / factorial(p)
        return y
end

function test_first_derivative(D::AbstractArray, xl::AbstractArray, xr::AbstractArray, p_int::Int64, p_bnd::Int64, ml::Int64, mr::Int64)

        @test isapprox(interior(D*xr.^0, ml, mr), interior(eps(Tv) .* xl, ml, mr), atol=1.0)
        for p in 1:p_int
                yr = taylor_poly(xr, p) 
                yl = taylor_poly(xl, p - 1) 
                @test isapprox(interior(D*yr, ml, mr), interior(yl, ml, mr))
        end
        
        @test isapprox(boundary(D*xr.^0, ml, mr), boundary(eps(Tv) .* xl, ml, mr), atol=1.0)
        for p in 1:p_bnd
                yr = taylor_poly(xr, p) 
                yl = taylor_poly(xl, p - 1) 
                @test isapprox(boundary(D*yr, ml, mr), boundary(yl, ml, mr))
        end

end

function test_interpolation(P::AbstractArray, xl::AbstractArray, xr::AbstractArray, p_int::Int64, p_bnd::Int64, ml::Int64, mr::Int64)

        for p in 0:p_int
                yr = taylor_poly(xr, p) 
                yl = taylor_poly(xl, p) 
                @test isapprox(interior(P*yr, ml, mr), interior(yl, ml, mr))
        end
        
        for p in 0:p_bnd
                yr = taylor_poly(xr, p) 
                yl = taylor_poly(xl, p) 
                @test isapprox(boundary(P*yr, ml, mr), boundary(yl, ml, mr))
        end

end

function test_quadrature(H::AbstractArray, x::AbstractArray, pv::Int64)

        for p in 0:pv
                xr = taylor_poly(x, p) 
                xl = taylor_poly(x, p + 1) 
                @test isapprox( sum(H*xr), xl[end] - xl[1])
        end
end

function interior(x::AbstractArray, ml::Int64, mr::Int64)
        return x[ml+1:size(x,1) - mr]
end


function boundary(x::AbstractArray, ml::Int64, mr::Int64)
        y = zeros(ml + mr)
        for i in 1:ml
                @inbounds y[i] = x[i]
        end
        m = size(x, 1)
        for i in 1:mr
                @inbounds y[i] = x[m - mr + i]
        end
        return y
end


end
