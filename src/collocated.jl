module Collocated
using SparseArrays, LinearAlgebra
import sbp
using sbp.Metrics: CovariantBasis
using sbp.Sparse: block_matrix, block_matrix_insert

struct Operators2D
        nx::Int64
        ny::Int64
        Dx::AbstractArray
        Dy::AbstractArray
        Hx::AbstractArray
        Hy::AbstractArray
        Bx::AbstractArray
        By::AbstractArray
end

function boundary_matrix(n, is_sparse=true)
        if is_sparse
                B = spzeros(n, n)
        else
                B = zeros(n, n)
        end

        B[1,1] = -1
        B[end,end] = 1

        return B
end

function build_operators_2d(Dx::AbstractArray, 
                            Dy::AbstractArray, 
                            Hx::AbstractArray,
                            Hy::AbstractArray,
                            Bx::AbstractArray,
                            By::AbstractArray
                            )
        nx = size(Dx, 1)
        ny = size(Dy, 1)

        Ix = sparse(1.0I, nx, nx)
        Iy = sparse(1.0I, ny, ny)

        Dx2 = kron(Dx, Iy)
        Dy2 = kron(Ix, Dy)

        Hx2 = kron(Hx, Iy)
        Hy2 = kron(Ix, Hy)

        Bx2 = kron(Bx, Iy)
        By2 = kron(Ix, By)

        return Operators2D(nx, ny, Dx2, Dy2, Hx2, Hy2, Bx2, By2)
                
end


"""
Construct the covariant basis at the cell centers (mm) using SBP operators.
The mapping function f(x,y) = (fx, fy) must be defined at the nodes (pp). 
"""
function build_covariant_basis(fx::AbstractArray, fy::AbstractArray,
                               op::Operators2D)
        x_r1 = op.Dx * fx
        x_r2 = op.Dy * fx
        y_r1 = op.Dx * fy
        y_r2 = op.Dy * fy

        return CovariantBasis(x_r1, y_r1, x_r2, y_r2)
end


end
