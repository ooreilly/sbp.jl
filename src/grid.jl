module Grid

function grid_xp(n::Int64)
        x = zeros(n)
        h = 1.0 / (n - 1)
        for i in 1:n
                x[i] = (i - 1) * h
        end

        return x, h
end

function grid_xm(n::Int64)
        x = zeros(n)
        h = 1.0 / (n - 2)
        for i in 1:n
                x[i] = (i - 1 - 0.5) * h
        end
        x[1] = 0.0
        x[end] = 1.0

        return x, h
end

"""
Construct a 2D grid vector from a 1D grid vector `x1`. The resulting vector
repeats all x-values: x1[1], x1[1], ... , x1[1], x1[2], x1[2] ... x1[2] `ny`
times.
"""
function grid_2d_x(x1::AbstractArray, ny::Int64)
        nx = size(x1, 1) 
        x = zeros(nx * ny) 
        for i=1:nx
                for j=1:ny
                        @inbounds x[j + (i - 1) * ny] = x1[i]
                end
        end
        return x
end

"""
Construct a 2D grid vector from a 1D grid vector `y1`. The resulting vector
repeats all y-values: y1[1], y1[1], ... , y1[ny], y1[1], y1[2] ... y1[ny] `nx`
times.
"""
function grid_2d_y(y1::AbstractArray, nx::Int64)
        ny = size(y1, 1) 
        y = zeros(nx * ny) 
        for i=1:nx
                for j=1:ny
                        @inbounds y[j + (i - 1) * ny] = y1[j]
                end
        end
        return y
end

end
