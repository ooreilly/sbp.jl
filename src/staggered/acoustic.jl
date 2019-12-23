module StaggeredAcoustic
using LinearAlgebra
using SparseArrays
import sbp
using sbp.Sparse: diag_block_matrix_2x2, block_matrix_2x2, block_matrix,
                  block_matrix_insert
using sbp.Staggered: build_all_operators_2d, 
                   build_covariant_basis_mm, 
                   build_covariant_basis_pm,
                   build_covariant_basis_mp
using sbp.Metrics: build_jacobian, 
                   build_contravariant_basis, 
                   contravariant_metric_tensor_to_matrix,
                   ContravariantMetricTensor     


struct AcousticOperators
        nx::Int64
        ny::Int64
        np::Int64
        n1::Int64
        n2::Int64
        Jp::AbstractArray
        J1::AbstractArray
        J2::AbstractArray
        Jip::AbstractArray
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
        G1::ContravariantMetricTensor
        G2::ContravariantMetricTensor
        B1::AbstractArray
        B2::AbstractArray
end

function init_operators(nx::Int64, ny::Int64, build_operators::Function,
                        fx::AbstractArray, fy::AbstractArray)
        pp, mm, pm, mp = build_all_operators_2d(build_operators, nx, ny)
        np = (nx + 1) * (ny + 1)
        n1 = nx * (ny + 1)
        n2 = (nx + 1) * ny
        ap = build_covariant_basis_mm(fx, fy, mm, pm)
        a1 = build_covariant_basis_pm(fx, fy, mm, pm, mp)
        a2 = build_covariant_basis_mp(fx, fy, mm, pm, mp)
        Jp = build_jacobian(ap)
        J1 = build_jacobian(a1)
        J2 = build_jacobian(a2)
        bp = build_contravariant_basis(Jp, ap)
        b1 = build_contravariant_basis(J1, a1)
        b2 = build_contravariant_basis(J2, a2)
        Jip = spdiagm(0 => 1 ./Jp)
        Ji1 = spdiagm(0 => 1 ./ J1)
        Ji2 = spdiagm(0 => 1 ./ J2)
        Jp = spdiagm(0 => Jp)
        J1 = spdiagm(0 => J1)
        J2 = spdiagm(0 => J2)
        Gp = sbp.Metrics.build_contravariant_metric_tensor(bp)
        Gp = contravariant_metric_tensor_to_matrix(Gp)
        G1 = sbp.Metrics.build_contravariant_metric_tensor(b1)
        G1 = contravariant_metric_tensor_to_matrix(G1)
        G2 = sbp.Metrics.build_contravariant_metric_tensor(b2)
        G2 = contravariant_metric_tensor_to_matrix(G2)
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
        Hxi = spdiagm(0 => 1 ./ diag(pm.Hx))
        Hyi = spdiagm(0 => 1 ./ diag(mp.Hy))
        B1 = Hxi * pm.Bx
        B2 = Hyi * mp.By
        return AcousticOperators(nx, ny, np, n1, n2, Jp, J1, J2, Jip, Ji1, Ji2,
                                     Pc1, P1c, Pc2, P2c, D1p, D2p, D1, D2, Hp,
                                     H1, H2, Gp, G1, G2, B1, B2)

end

function grid_points(field::String, nx::Int64, ny::Int64)
        if field == "p"
                return nx+1, ny+1
        elseif field == "v1"
                return nx+1, ny
        elseif field == "v2"
                return nx, ny+1
        elseif field == "node"
                return nx, ny
        end
end


function grids(field::String, nx::Int64, ny::Int64; pad=false:Bool)
        xp1, h = sbp.Grid.grid_xp(nx, pad=pad)
        yp1, h = sbp.Grid.grid_xp(ny, pad=pad)

        xm1, h = sbp.Grid.grid_xm(nx+1)
        ym1, h = sbp.Grid.grid_xm(ny+1)

        if field == "p"
                xm = sbp.Grid.grid_2d_x(xm1, length(ym1))
                ym = sbp.Grid.grid_2d_y(ym1, length(xm1))
                return xm, ym
        elseif field == "v1"
                xp = sbp.Grid.grid_2d_x(xp1, length(ym1))
                ym = sbp.Grid.grid_2d_y(ym1, length(xp1))
                return xp, ym
        elseif field == "v2"
                xm = sbp.Grid.grid_2d_x(xm1, length(yp1))
                yp = sbp.Grid.grid_2d_y(yp1, length(xm1))
                return xm, yp
        elseif field == "node"
                xp = sbp.Grid.grid_2d_x(xp1, length(yp1))
                yp = sbp.Grid.grid_2d_y(yp1, length(xp1))
                return xp, yp
        end
end

        

function pressure_norm(op::AcousticOperators)
        return op.Hp * op.Jp
end

function velocity_norm(op::AcousticOperators)
        H = diag_block_matrix_2x2(op.H1, op.H2)
        J = diag_block_matrix_2x2(op.J1, op.J2)
        return H * J
end

function energy_norm(HJp::AbstractArray, HJ::AbstractArray, Gp::AbstractArray)
        rows = [size(HJp,1), size(HJ,1)]
        E = block_matrix(rows, rows)
        E = block_matrix_insert(E, rows, rows, 1, 1, HJp)
        Gip = inv(Matrix(Gp))
        E = block_matrix_insert(E, rows, rows, 2, 2, HJ * Gip)
        return E
end

function contravariant_metric_tensor(op::AcousticOperators)
        G11 = op.Ji1 * op.P1c * op.Jp * op.Gp.b11 * op.Pc1
        G12 = op.Ji1 * op.P1c * op.Jp * op.Gp.b12 * op.Pc2
        G21 = op.Ji2 * op.P2c * op.Jp * op.Gp.b21 * op.Pc1
        G22 = op.Ji2 * op.P2c * op.Jp * op.Gp.b22 * op.Pc2
        return block_matrix_2x2(G11, G12, G21, G22)
end

function modified_contravariant_metric_tensor(op::AcousticOperators)
        G11 = op.G1.b11
        G12 = op.Ji1 * op.P1c * op.Jp * op.Gp.b12 * op.Pc2
        G21 = op.Ji2 * op.P2c * op.Jp * op.Gp.b21 * op.Pc1
        G22 = op.G2.b22
        return block_matrix_2x2(G11, G12, G21, G22)
end

function divergence(op::AcousticOperators)
        rows = [op.np]
        cols = [op.n1, op.n2]
        Ap = block_matrix(rows, cols)
        Ap1 = op.Jip * op.D1p * op.J1
        Ap2 = op.Jip * op.D2p * op.J2
        Ap = block_matrix_insert(Ap, rows, cols, 1, 1, Ap1)
        Ap = block_matrix_insert(Ap, rows, cols, 1, 2, Ap2)
        return Ap
end

function gradient(op::AcousticOperators, G::AbstractArray)
        rows = [op.n1, op.n2]
        cols = [op.np]
        Av = block_matrix(rows, cols)
        A1 = op.D1
        A2 = op.D2
        Av = block_matrix_insert(Av, rows, cols, 1, 1, A1)
        Av = block_matrix_insert(Av, rows, cols, 2, 1, A2)
        return G * Av

end

function spatial_discretization(Ap::AbstractArray, Av::AbstractArray)
        rows = [size(Ap,1), size(Av,1)]
        cols = [size(Av,2), size(Ap,2)]
        A = block_matrix(rows, cols)
        A = block_matrix_insert(A, rows, cols, 1, 2, Ap)
        A = block_matrix_insert(A, rows, cols, 2, 1, Av)
        return A
end

function pressure_sat(op::AcousticOperators, G::AbstractArray)
        rows = [op.n1, op.n2]
        cols = [op.np]
        S = block_matrix(rows, cols)
        S = block_matrix_insert(S, rows, cols, 1, 1, op.B1)
        S = block_matrix_insert(S, rows, cols, 2, 1, op.B2)
        return G * S
end

function pressure_source(op::AcousticOperators, 
                         r1_s::Float64,
                         r2_s::Float64,
                         num_moment::Int64,
                         num_smoothness::Int64)
        rows = [op.np, op.n1, op.n2]
        cols = [1]
        r1, h1 = sbp.Grid.grid_xm(op.nx + 1)
        r2, h2 = sbp.Grid.grid_xm(op.ny + 1)
        d = sbp.Source.source_discretize_2d(r1_s, r2_s, num_moment,
                                            num_smoothness, r1, r2, h1,
                                            h2) 
        S = block_matrix(rows, cols)
        S = block_matrix_insert(S, rows, cols, 1, 1, d)
        return S
end

function pressure_receiver(op::AcousticOperators, 
                         r1_s::Float64,
                         r2_s::Float64,
                         num_moment::Int64,
                         num_smoothness::Int64)
        rows = [op.np, op.n1, op.n2]
        cols = [1]
        r1, h1 = sbp.Grid.grid_xm(op.nx + 1)
        r2, h2 = sbp.Grid.grid_xm(op.ny + 1)
        r = sbp.Source.source_discretize_2d(r1_s, r2_s, num_moment,
                                            num_smoothness, r1, r2, 1.0,
                                            1.0) 
        S = block_matrix(rows, cols)
        S = block_matrix_insert(S, rows, cols, 1, 1, r)
        return S
end

end
