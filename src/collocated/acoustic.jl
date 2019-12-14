module CollocatedAcoustic
using LinearAlgebra
using SparseArrays
import sbp
using sbp.Operator: OperatorData
using sbp.Sparse: diag_block_matrix_2x2, block_matrix_2x2, block_matrix,
                  block_matrix_insert
using sbp.Collocated: build_operators_2d, boundary_matrix,
                   build_covariant_basis
using sbp.Metrics: build_jacobian, 
                   build_contravariant_basis, 
                   contravariant_metric_tensor_to_matrix,
                   ContravariantMetricTensor     


struct AcousticOperators
        n::Int64
        nx::Int64
        ny::Int64
        J::AbstractArray
        Ji::AbstractArray
        D1::AbstractArray
        D2::AbstractArray
        H::AbstractArray
        G::ContravariantMetricTensor
        B1::AbstractArray
        B2::AbstractArray
end

function init_operators(nx::Int64, ny::Int64, builder::Module, desc_D::OperatorData,
                        desc_H::OperatorData, fx::AbstractArray,
                        fy::AbstractArray)

        x, h, Dx, Hx = builder.build_operators(desc_D, desc_H, nx)
        y, h, Dy, Hy = builder.build_operators(desc_D, desc_H, ny)
        Bx = boundary_matrix(nx)
        By = boundary_matrix(ny)
        op = build_operators_2d(Dx, Dy, Hx, Hy, Bx, By)
        n = nx * ny
        n1 = nx
        n2 = ny
        a = build_covariant_basis(fx, fy, op)
        J = build_jacobian(a)
        b = build_contravariant_basis(J, a)
        Ji = spdiagm(0 => 1 ./J)
        J = spdiagm(0 => J)
        G = sbp.Metrics.build_contravariant_metric_tensor(b)
        G = contravariant_metric_tensor_to_matrix(G)
        D1 = op.Dx
        D2 = op.Dy
        H = op.Hx * op.Hy
        Hxi = spdiagm(0 => 1 ./ diag(op.Hx))
        Hyi = spdiagm(0 => 1 ./ diag(op.Hy))
        B1 = Hxi * op.Bx
        B2 = Hyi * op.By
        return AcousticOperators(n, nx, ny, J, Ji, D1, D2, H, G, B1, B2)

end

function grid(nx::Int64, ny::Int64)
        x1, h = sbp.Grid.grid_xp(nx)
        y1, h = sbp.Grid.grid_xp(ny)
        x = sbp.Grid.grid_2d_x(x1, ny)
        y = sbp.Grid.grid_2d_y(y1, nx)
        return x, y
end

function pressure_norm(op::AcousticOperators)
        return op.H * op.J
end

function velocity_norm(op::AcousticOperators)
        H = diag_block_matrix_2x2(op.H, op.H)
        J = diag_block_matrix_2x2(op.J, op.J)
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
        G11 = op.G.b11
        G12 = op.G.b12
        G21 = op.G.b21
        G22 = op.G.b22
        return block_matrix_2x2(G11, G12, G21, G22)
end

function divergence(op::AcousticOperators)
        rows = [op.n]
        cols = [op.n, op.n]
        Ap = block_matrix(rows, cols)
        Ap1 = op.Ji * op.D1 * op.J
        Ap2 = op.Ji * op.D2 * op.J
        Ap = block_matrix_insert(Ap, rows, cols, 1, 1, Ap1)
        Ap = block_matrix_insert(Ap, rows, cols, 1, 2, Ap2)
        return Ap
end

function gradient(op::AcousticOperators, G::AbstractArray)
        rows = [op.n, op.n]
        cols = [op.n]
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
        rows = [op.n, op.n]
        cols = [op.n]
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
        rows = [op.n, op.n, op.n]
        cols = [1]
        r1, h1 = sbp.Grid.grid_xp(op.nx)
        r2, h2 = sbp.Grid.grid_xp(op.ny)
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
        rows = [op.n, op.n, op.n]
        cols = [1]
        r1, h1 = sbp.Grid.grid_xm(op.nx)
        r2, h2 = sbp.Grid.grid_xm(op.ny)
        r = sbp.Source.source_discretize_2d(r1_s, r2_s, num_moment,
                                            num_smoothness, r1, r2, 1.0,
                                            1.0) 
        S = block_matrix(rows, cols)
        S = block_matrix_insert(S, rows, cols, 1, 1, r)
        return S
end

end

