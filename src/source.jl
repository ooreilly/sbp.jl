module Source

using SparseArrays


"""
Determine the indices and values of the source discretization stencil that
spans `p` cells.

Input arguments:
x : List of coordinates
xs : Source location
p : Polynomial degree. 
s : Number of smoothness conditions. 

Return values:
stencil : An array of the stencil indices that maps to the `x` coordinates.
coef : An array of stencil coefficients (values).
alpha : The signed distance from the source location to the nearest grid point
"""
function source_smoothness(x, xs, p, s)

        idx = argmin(abs.(x .- xs))
        
        # Determine signed distance from source location to nearest point.
        alpha = x[idx] - xs
        
        shift::Int64 = 0
        if alpha < 0
                shift = 0
        else
                shift = -1
        end
        
        if p == 0
                coef = 1
                stencil = idx
                return
        end
        
        if mod(p,1) == 1
                error("Polynomial degree must be odd")
        end
        
        left = round((p-1)/2) + round(s/2);
        right = p + s - left - 1;
        
        idx_l::Int64 = idx - left + shift;
        idx_r::Int64 = idx + right + shift + 1; 
        
        adjust::Int64 = 0;
        if idx_l < 1
                adjust = 1 - idx_l;
        end
        if idx_r > length(x)
                adjust = length(x) - idx_r;
        end
        
        stencil = collect((adjust + idx_l):(adjust + idx_r));

        # Moment conditions
        e = ones(size(stencil))
        xl = x[stencil] .- x[idx];
        M = zeros(p+1, size(xl, 1))
        for i in 1:p+1
        for j in 1:size(xl, 1)
                M[i,j] = xl[j] ^ (i - 1)
        end
        end
        b1 = (xs - x[idx]).^(0:p); 

        if s > 0
                # Smoothness conditions
                S = zeros(s,size(xl, 1))
                for i in 1:s
                for j in 1:size(xl, 1)
                        S[i,j] = (-1)^j * xl[j] ^ (i - 1)
                end
                end
                A = [M; S]
                b = [b1; zeros(s,1)]
        else
                A = M
                b = b1
        end
        coef = A\b;
return stencil, coef, alpha
end


function source_discretize_1d(r1_s, p, s, r1, h1)
        d = spzeros(length(r1))
        stencil1, coef1, alpha1 = source_smoothness(r1, r1_s, p, s)
        for i=1:length(stencil1)
                d[stencil1[i]] = coef1[i] / h1
        end

        return d
end
"""
d = source_disc(xs,ys, h, gr, J)
Discretize point source using nearest grid point

Arguments:
xs, ys : Source location
h : Grid spacing in parameter space
gr : grid object that contains X,Y meshgrids.
J : Jacobian

Returns
A sparse vector that is zero everywhere except for one index, which is set to
a weight d = 1/(cell-area).
"""
function source_discretize_2d(r1_s, r2_s, p, s, r1, r2, h1, h2)
        d = spzeros(length(r1)*length(r2))
        stencil1, coef1, alpha1 = source_smoothness(r1, r1_s, p, s)
        stencil2, coef2, alpha2 = source_smoothness(r2, r2_s, p, s)

        n2 = length(r2)
        reci_area = 1 / (h1 * h2)
        for i=1:length(stencil1)
        for j=1:length(stencil2)
          d[stencil2[j] + n2 * stencil1[i]] = coef1[i] * coef2[j] * reci_area
        end
        end

        return d
end

end
