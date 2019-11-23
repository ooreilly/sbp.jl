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

function build_jacobian(a::CovariantBasis)
        return a.x_r1 .* a.y_r2 - a.x_r2 .* a.y_r1
end

function build_contravariant_basis(J::AbstractArray, a::CovariantBasis)
        return ContravariantBasis( a.y_r2 ./ J, 
                                 - a.y_r1 ./ J, 
                                 - a.x_r2 ./ J, 
                                   a.x_r1 ./ J)

end


function covariant_metric_tensor()


end

end
