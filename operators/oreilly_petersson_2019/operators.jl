module OP2019

import sbp

using sbp.Operator: read_operator, to_matrix, OperatorData
using sbp.Grid:grid_xp, grid_xm

default_path = string(@__DIR__, "/")


# First derivative operators
ODp = read_operator(string(default_path,  "Dp_42.txt"))
ODm = read_operator(string(default_path,  "Dm_42.txt"))

# Quadrature rules
OHp = read_operator(string(default_path,  "Hp_42.txt"))
OHm = read_operator(string(default_path,  "Hm_42.txt"))

# Interpolation operators
OPp = read_operator(string(default_path,  "Pp_42.txt"))
OPm = read_operator(string(default_path,  "Pm_42.txt"))


function build_operators(m::Int64, is_sparse::Bool=true)

        xm, h = grid_xm(m+1)
        xp, h = grid_xp(m)

        Dp = to_matrix(ODp, m, m + 1, is_sparse) / h
        Dm = to_matrix(ODm, m + 1, m, is_sparse) / h

        Hp = to_matrix(OHp, m, m, is_sparse) * h
        Hm = to_matrix(OHm, m + 1, m + 1, is_sparse) * h

        Pp = to_matrix(OPp, m, m + 1, is_sparse)
        Pm = to_matrix(OPm, m + 1, m, is_sparse)

        return xp, xm, h, Dp, Dm, Hp, Hm, Pp, Pm

end




end
