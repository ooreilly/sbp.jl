module Staggered
using SparseArrays, LinearAlgebra
import sbp
using sbp.Operator: Operators2D
using sbp.Metrics: CovariantBasis
using sbp.Sparse: block_matrix, block_matrix_insert

function boundary_matrix_p(n, is_sparse=true)
        if is_sparse
                B = spzeros(n, n + 1)
        else
                B = zeros(n, n + 1)
        end

        B[1,1] = -1
        B[end,end] = 1

        return B
end

function boundary_matrix_m(n, is_sparse=true)
        return boundary_matrix_p(n - 1, is_sparse)'
end

function build_operators_2d(Dx::AbstractArray, 
                            Dy::AbstractArray, 
                            Px::AbstractArray, 
                            Py::AbstractArray,
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

        Px2 = kron(Px, Iy)
        Py2 = kron(Ix, Py)

        Hx2 = kron(Hx, Iy)
        Hy2 = kron(Ix, Hy)

        Bx2 = kron(Bx, Iy)
        By2 = kron(Ix, By)

        return Operators2D(nx, ny, Dx2, Dy2, Px2, Py2, Hx2, Hy2, Bx2, By2)
                
end

function build_all_operators_2d(builder::Function, nx::Int64, ny::Int64)
        xp, xm, h, Dxp, Dxm, Hxp, Hxm, Pxp, Pxm = builder(nx)
        yp, ym, h, Dyp, Dym, Hyp, Hym, Pyp, Pym = builder(ny)

        Bxp = boundary_matrix_p(nx)
        Bxm = boundary_matrix_m(nx+1)
        Byp = boundary_matrix_p(ny)
        Bym = boundary_matrix_m(ny+1)

        pp = build_operators_2d(Dxp, Dyp, Pxp, Pyp, Hxp, Hyp, Bxp, Byp)
        mm = build_operators_2d(Dxm, Dym, Pxm, Pym, Hxm, Hym, Bxm, Bym)
        pm = build_operators_2d(Dxp, Dym, Pxp, Pym, Hxp, Hym, Bxp, Bym)
        mp = build_operators_2d(Dxm, Dyp, Pxm, Pyp, Hxm, Hyp, Bxm, Byp)

        return pp, mm, pm, mp
end


"""
Construct the covariant basis at the cell centers (mm) using SBP operators.
The mapping function f(x,y) = (fx, fy) must be defined at the nodes (pp). 
"""
function build_covariant_basis_mm(fx::AbstractArray, fy::AbstractArray,
                               mm::Operators2D, pm::Operators2D)
        x_r1 = mm.Dx * pm.Py * fx
        x_r2 = mm.Px * pm.Dy * fx
        y_r1 = mm.Dx * pm.Py * fy
        y_r2 = mm.Px * pm.Dy * fy

        return CovariantBasis(x_r1, y_r1, x_r2, y_r2)
end

"""
Construct the covariant basis at the v^1 components (pm) using SBP operators.
The mapping function f(x,y) = (fx, fy) must be defined at the nodes (pp). 
"""
function build_covariant_basis_pm(fx::AbstractArray, fy::AbstractArray,
                               mm::Operators2D, pm::Operators2D,
                               mp::Operators2D)
        x_r1 = pm.Px * mm.Py * mp.Dx * fx
        x_r2 = pm.Dy * fx
        y_r1 = pm.Px * mm.Py * mp.Dx * fy
        y_r2 = pm.Dy * fy

        return CovariantBasis(x_r1, y_r1, x_r2, y_r2)
end

"""
Construct the covariant basis at the v^2 components (mp) using SBP operators.
The mapping function f(x,y) = (fx, fy) must be defined at the nodes (pp). 
"""
function build_covariant_basis_mp(fx::AbstractArray, fy::AbstractArray,
                               mm::Operators2D, pm::Operators2D,
                               mp::Operators2D)
        x_r1 = mp.Dx * fx
        x_r2 = mp.Py * mm.Px * pm.Dy * fx
        y_r1 = mp.Dx * fy
        y_r2 = mp.Py * mm.Px * pm.Dy * fy

        return CovariantBasis(x_r1, y_r1, x_r2, y_r2)
end


end
