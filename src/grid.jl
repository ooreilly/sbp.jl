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

end
