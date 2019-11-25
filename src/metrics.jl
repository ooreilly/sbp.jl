module Metrics


struct CovariantBasis
        # a1
        x_r1::AbstractArray
        y_r1::AbstractArray
        # a2
        x_r2::AbstractArray
        y_r2::AbstractArray
end

struct ContravariantBasis
        # b1 
        r1_x::AbstractArray 
        r1_y::AbstractArray 
        # b2 
        r2_x::AbstractArray 
        r2_y::AbstractArray 
end


struct CovariantMetricTensor
        a11::AbstractArray
        a12::AbstractArray
        a21::AbstractArray
        a22::AbstractArray
end

struct ContravariantMetricTensor
        b11::AbstractArray
        b12::AbstractArray
        b21::AbstractArray
        b22::AbstractArray
end

function build_jacobian(a::CovariantBasis)
        return a.x_r1 .* a.y_r2 - a.x_r2 .* a.y_r1
end

function build_contravariant_basis(J::AbstractArray, a::CovariantBasis)
        return ContravariantBasis( a.y_r2 ./ J, 
                                 - a.y_r1 ./ J, 
                                 - a.x_r2 ./ J, 
                                   a.x_r1 ./ J)

end

function build_covariant_metric_tensor(a::CovariantBasis)
        return CovariantMetricTensor(a.x_r1 .* a.x_r1 + a.y_r1 .* a.y_r1, 
                                     a.x_r1 .* a.x_r2 + a.y_r1 .* a.y_r2,
                                     a.x_r2 .* a.x_r1 + a.y_r2 .* a.y_r1,
                                     a.x_r2 .* a.x_r2 + a.y_r2 .* a.y_r2)
end

function build_contravariant_metric_tensor(b::ContravariantBasis)
        return ContravariantMetricTensor(b.r1_x .* b.r1_x + b.r1_y .* b.r1_y, 
                                         b.r1_x .* b.r2_x + b.r1_y .* b.r2_y,
                                         b.r2_x .* b.r1_x + b.r2_y .* b.r1_y,
                                         b.r2_x .* b.r2_x + b.r2_y .* b.r2_y)
end

end
