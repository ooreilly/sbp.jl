module BoundaryOperator
"""

Module for constructing the primal boundary operator 
L = X^- - X^+ R_1^T,
the dual operator La, and their complementary operators. 

In addition, this module can construct the projection operator P that gives 
\$ \\delta u = P u \$.

"""

using LinearAlgebra

function virtual_solution_matrix(dLp::Matrix, L::Matrix)
        return dLp * inv(L' * dLp) * L'
end

function dual_consistency_condition(dS::Matrix, A::Matrix, dSa::Matrix)
        return dS' * A * dSa
end

function boundary_operator(Xp::Matrix, Xm::Matrix, R1)
        L = Xm - Xp * transpose(R1)
        return L
end

function complementary_boundary_operator(Xp::Matrix, Xm::Matrix, R1::Matrix)
        Lp = Xp + Xm * R1
        return Lp
end

function dual_boundary_operator(Xp::Matrix, Xm::Matrix, R1a::Matrix)
        La = Xp - Xm * R1a
        return La
end

function complementary_dual_boundary_operator(Xp::Matrix, Xm::Matrix, R1a::Matrix)
        Lap = Xm + Xp * transpose(R1a)
        return Lap
end

function virtual_boundary_operator(Xp::Matrix, Xm::Matrix, dR1::Matrix)
        dL = dual_boundary_operator(Xp, Xm, dR1)
        return dL
end

function virtual_complementary_boundary_operator(Xp::Matrix, Xm::Matrix, dR1::Matrix)
        dLp = complementary_dual_boundary_operator(Xp, Xm, dR1)
        return dLp
end

function virtual_dual_boundary_operator(Xp::Matrix, Xm::Matrix, dR1a::Matrix)
        dLa = boundary_operator(Xp, Xm, dR1a)
        return dLa
end

function virtual_complementary_dual_boundary_operator(Xp::Matrix, Xm::Matrix, dR1a::Matrix)
        dLap = complementary_boundary_operator(Xp, Xm, dR1a)
        return dLap
end

function dual_reflection_matrix(Lp::Matrix, Lm::Matrix, R1::Matrix)
        R = scale_reflection_matrix(Lp, Lm, R1)
        Ra = R
        return inverse_scale_dual_reflection_matrix(Lp, Lm, Ra)
end

function scale_reflection_matrix(Lp::Matrix, Lm::Matrix, R1::Matrix)
        R =  abs.(Lm) .^ (1/2) * R1 * Lp .^ (-1/2)
        return R
end

function inverse_scale_reflection_matrix(Lp::Matrix, Lm::Matrix, R::Matrix)
        R1 =  abs.(Lm) .^ (-1/2) * R * Lp .^ (1/2)
        return R1
end

function scale_dual_reflection_matrix(Lp::Matrix, Lm::Matrix, R1a)
        Ra = abs.(Lm) .^ (-1/2) * R1a * Lp .^ (1/2)
        return Ra
end

function inverse_scale_dual_reflection_matrix(Lp::Matrix, Lm::Matrix, Ra)
        R1a = abs.(Lm) .^ (1/2) * Ra * Lp .^ (-1/2)
        return R1a
end
function virtual_scale_reflection_matrix(Lp::Matrix, Lm::Matrix, dR1::Matrix)
        dR = scale_dual_reflection_matrix(Lp, Lm, dR1)
        return dR
end

function virtual_inverse_scale_reflection_matrix(Lp::Matrix, Lm::Matrix, dR::Matrix)
        dR1 = inverse_scale_dual_reflection_matrix(Lp, Lm, dR)
        return dR1
end

function virtual_scale_dual_reflection_matrix(Lp::Matrix, Lm::Matrix, dR1a::Matrix)
        dRa = inverse_scale_reflection_matrix(Lp, Lm, dR1a)
        return dRa
end

function virtual_inverse_scale_dual_reflection_matrix(Lp::Matrix, Lm::Matrix, dRa::Matrix)
        dR1a = scale_reflection_matrix(Lp, Lm, dRa)
        return dR1a
end

function virtual_dual_reflection_matrix(Lp::Matrix, Lm::Matrix, dR1::Matrix)
        dR = virtual_scale_reflection_matrix(Lp, Lm, dR1)
        dRa = dR
        dR1a = virtual_inverse_scale_dual_reflection_matrix(Lp, Lm, dRa)
        return dR1a
end

function is_bounded_operator(R::Matrix)
        return norm(R) <= 1
end

function eigen_positive(A::Array)
        return eigen_split(A, x -> x > 1e-12)
end

function eigen_negative(A::Array)
        return eigen_split(A, x -> x < -1e-12)
end

function eigen_split(A::Array, comp::Function)

        F = eigen(A)
        X = F.vectors
        eig = F.values

        num_eig_s = 0

        for i in 1:length(eig)
                if comp(eig[i])
                        num_eig_s += 1
                end
        end

        if num_eig_s == 0
                num_eig_s = 1
        end

        Xs = zeros(size(X, 1), num_eig_s)
        Ls = zeros(num_eig_s, num_eig_s)
        
        k = 1
        for i in 1:length(eig)
                if comp(eig[i])
                        Xs[:,k] = X[:,i]
                        Ls[k,k] = eig[i]
                        k += 1
                end
        end

        return Xs, Ls

end

end
