"""
Discretize the acoustic wave equation in 1D staggered SBP-SAT
"""

using SparseArrays
using LinearAlgebra
using Test

import sbp

# Number of grid points
nx = 100

"""
Load staggered grid operators and grids
xp denotes the regular grid with equidistantly spaced grid points, and xm
denotes the cell-centered grid, including boundary points.
"""
ops = sbp.OP2019
desc = ops.load_description()
xp, xm, h, Dp, Dm, Hp, Hm, Pp, Pm = ops.build_operators(desc, nx)

"""

Write the discretization as dq/dt = A * q, q = [v, p], where 

A = [ 0, Dp, 
      Dm, 0], 

v is stored on the xp grid, and p is stored on the xm grid.

"""
rows = [length(xp), length(xm)]
cols = [length(xp), length(xm)]

A = sbp.Sparse.block_matrix(rows, cols)
A12 = Dp
A21 = Dm
A = sbp.Sparse.block_matrix_insert(A, rows, rows, 1, 2, A12)
A = sbp.Sparse.block_matrix_insert(A, rows, rows, 2, 1, A21)

"""
Impose the boundary condition p = 0 on each boundary weakly in an energy
conserving manner.
"""
Bp = sbp.Staggered.boundary_matrix_p(nx)
Hpi = spdiagm(0 => 1 ./diag(Hp))
SAT = -Hpi*Bp

A = sbp.Sparse.block_matrix_insert(A, rows, rows, 1, 2, A12 + SAT)

"""
Check that the discretization is energy conservative to machine precision by
eigenvalue analysis

"""
lam = real(eigvals(Matrix(A)))
@test minimum(lam) > -1e-13
@test maximum(lam) < 1e-13
