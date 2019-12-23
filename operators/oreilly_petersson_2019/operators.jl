module OP2019

import sbp

using sbp.Operator: read_operator, to_matrix, OperatorData
using sbp.Grid:grid_xp, grid_xm

default_path = string(@__DIR__, "/")

struct Description
        Dp::OperatorData
        Hp::OperatorData
        Pp::OperatorData
        Dm::OperatorData
        Hm::OperatorData
        Pm::OperatorData
end

function load_description(;path=""::String, order=4::Int64, verbose=false::Bool)
        if path == ""
                path = default_path
        end

        if order==4
                desc_Dp = sbp.Operator.read_operator(string(path, "Dp_42.txt"))
                desc_Hp = sbp.Operator.read_operator(string(path, "Hp_42.txt"))
                desc_Pp = sbp.Operator.read_operator(string(path, "Pp_42.txt"))
                desc_Dm = sbp.Operator.read_operator(string(path, "Dm_42.txt"))
                desc_Hm = sbp.Operator.read_operator(string(path, "Hm_42.txt"))
                desc_Pm = sbp.Operator.read_operator(string(path, "Pm_42.txt"))
        else
           error("No operators implemented for order = ", order)
        end

        return Description(desc_Dp, desc_Hp, desc_Pp, desc_Dm, desc_Hm, desc_Pm)
end

function build_operators(desc, m::Int64, is_sparse::Bool=true)

        xm, h = grid_xm(m+1)
        xp, h = grid_xp(m)

        Dp = to_matrix(desc.Dp, m, m + 1, is_sparse) / h
        Dm = to_matrix(desc.Dm, m + 1, m, is_sparse) / h

        Hp = to_matrix(desc.Hp, m, m, is_sparse) * h
        Hm = to_matrix(desc.Hm, m + 1, m + 1, is_sparse) * h

        Pp = to_matrix(desc.Pp, m, m + 1, is_sparse)
        Pm = to_matrix(desc.Pm, m + 1, m, is_sparse)

        return xp, xm, h, Dp, Dm, Hp, Hm, Pp, Pm

end

end
