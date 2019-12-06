import sbp
using sbp.Source: source_smoothness, source_discretize_1d, source_discretize_2d
using Test

tol=1e-12;
n = 8;
x = collect(1:n);
N = zeros(n)
for i=1:n
        N[i] = (-1) ^ i
end
        
moments = (coef, stencil, p, xs) -> abs(sum(coef.*x[stencil].^p) - xs^p);
smoothness = (coef, stencil, s, xs) -> abs(sum(coef.*N[stencil].*x[stencil].^s));

"""
 Source stencil spanning 3 cells (nearest stencil point to the left)

             |     *     |
 o---o---o---o---o---o---o---o
 1   2   3   4   5   6   7   8
 * : source location
"""

@testset "Source stencil spanning 3 cells, left" begin
p = 3
s = 0
xs = 5.5
stencil, coef, alpha = source_smoothness(x, xs, p, s)
for i=0:p
        @test (moments(coef, stencil, i, xs) < tol) 
end
@test (alpha == -0.5)
@test (length(stencil) == 4)
@test (stencil[1] == 4)
end

"""
 Source stencil spanning 3 cells (nearest stencil point to the right)
             |     *     |
 o---o---o---o---o---o---o---o
 1   2   3   4   5   6   7   8
 * : source location
"""

@testset "Source stencil spanning 3 cells, right" begin
p = 3
s = 0
xs = 5.6
stencil, coef, alpha = source_smoothness(x, xs, p, s)
for i=0:p
        @test (moments(coef, stencil, i, xs) < tol) 
end
@test (abs(alpha - 0.4) < tol)
@test (length(stencil) == 4)
@test (stencil[1] == 4)
end

"""
 Source stencil spanning 3 cells (nearest stencil point close to the right
 boundary)

                 |          *|
 o---o---o---o---o---o---o---o
 1   2   3   4   5   6   7   8
 * : source location
"""

@testset "Source stencil spanning 3 cells, right boundary" begin
p = 3
s = 0
xs = 7.7
stencil, coef, alpha = source_smoothness(x, xs, p, s)
for i=0:p
        @test (moments(coef, stencil, i, xs) < tol)
end
@test (abs(alpha - 0.3) < tol)
@test (length(stencil) == 4)
@test (stencil[1] == 5)
end

"""
 Test of smoothness conditions
 Source stencil spanning 4 cells (nearest stencil point to the right)
     |         *         |
 o---o---o---o---o---o---o---o
 1   2   3   4   5   6   7   8
 * : source location
"""

@testset "Source stencil spanning 3 cells, smoothness conditions" begin
p = 3
s = 2;
xs = 4.5;
stencil, coef, alpha = source_smoothness(x, xs, p, s);
for i=0:p
        @test (moments(coef, stencil, i, xs) < tol) 
end
for i=0:(s-1)
        @test (smoothness(coef, stencil, i, xs) < tol) 
end
end

@testset "1D source discretization" begin
nx = 10
xp, h = sbp.Grid.grid_xp(nx)
xs = 0.5
p = 3
s = 2
d = source_discretize_1d(xs, p, s, xp, h) 
@test isapprox(sum(d) * h, 1)

end

@testset "2D source discretization" begin
nx = 10
ny = 10
xp, h1 = sbp.Grid.grid_xp(nx)
yp, h2 = sbp.Grid.grid_xp(ny)
xm, h1 = sbp.Grid.grid_xm(nx + 1)
ym, h2 = sbp.Grid.grid_xm(ny + 1)
xs = 0.5
ys = 0.5
p = 3
s = 2
d = source_discretize_2d(xs, ys, p, s, xp, yp, h1, h2) 
@test isapprox(sum(d) * h1 * h2, 1)

end

