module StaggeredAcoustic
using SparseArrays
import sbp
using sbp.Sparse: diag_block_matrix_2x2, block_matrix_2x2      
using sbp.Staggered: build_all_operators_2d, 
                   build_covariant_basis_mm, 
                   build_covariant_basis_pm,
                   build_covariant_basis_mp
using sbp.Metrics: build_jacobian, 
                   build_contravariant_basis, 
                   contravariant_metric_tensor_to_matrix,
                   ContravariantMetricTensor     


struct AcousticOperators
        Jp::AbstractArray
        J1::AbstractArray
        J2::AbstractArray
        Ji1::AbstractArray
        Ji2::AbstractArray
        Pc1::AbstractArray
        P1c::AbstractArray
        Pc2::AbstractArray
        P2c::AbstractArray
        D1p::AbstractArray
        D2p::AbstractArray
        D1::AbstractArray
        D2::AbstractArray
        Hp::AbstractArray
        H1::AbstractArray
        H2::AbstractArray
        Gp::ContravariantMetricTensor
end

function init_operators(nx::Int64, ny::Int64, build_operators::Function,
                        mapping::Function)
        fx, fy = mapping(nx, ny)
        pp, mm, pm, mp = build_all_operators_2d(build_operators, nx, ny)
        ap = build_covariant_basis_mm(fx, fy, mm, pm)
        a1 = build_covariant_basis_pm(fx, fy, mm, pm, mp)
        a2 = build_covariant_basis_mp(fx, fy, mm, pm, mp)
        Jp = build_jacobian(ap)
        J1 = build_jacobian(a1)
        J2 = build_jacobian(a2)
        bp = build_contravariant_basis(Jp, ap)
        Jp = spdiagm(0 => Jp)
        Ji1 = spdiagm(0 => 1 ./ J1)
        Ji2 = spdiagm(0 => 1 ./ J2)
        J1 = spdiagm(0 => J1)
        J2 = spdiagm(0 => J2)
        Gp = sbp.Metrics.build_contravariant_metric_tensor(bp)
        Gp = contravariant_metric_tensor_to_matrix(Gp)
        P1c = pm.Px
        P2c = mp.Py
        Pc1 = mm.Px
        Pc2 = mm.Py
        D1p = mm.Dx
        D2p = mm.Dy
        D1 = pm.Dx
        D2 = mp.Dy
        Hp = mm.Hx * mm.Hy
        H1 = pm.Hx * pm.Hy
        H2 = mp.Hx * mp.Hy
        return AcousticOperators(Jp, J1, J2, Ji1, Ji2, Pc1, P1c, Pc2, P2c, D1p,
                                     D2p, D1, D2, Hp, H1, H2, Gp)

end

function pressure_norm(op::AcousticOperators)
        return op.Hp * op.Jp
end

function velocity_norm(op::AcousticOperators)
        H = diag_block_matrix_2x2(op.H1, op.H2)
        J = diag_block_matrix_2x2(op.J1, op.J2)
        return H * J
end

function contravariant_metric_tensor(op::AcousticOperators)
        G11 = op.Ji1 * op.P1c * op.Jp * op.Gp.b11 * op.Pc1
        G12 = op.Ji1 * op.P1c * op.Jp * op.Gp.b12 * op.Pc2
        G21 = op.Ji2 * op.P2c * op.Jp * op.Gp.b21 * op.Pc1
        G22 = op.Ji2 * op.P2c * op.Jp * op.Gp.b22 * op.Pc2
        return block_matrix_2x2(G11, G12, G21, G22)
end

end
