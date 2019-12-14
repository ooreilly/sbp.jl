module Strand1994

import sbp
using sbp.Operator: OperatorData, to_matrix

default_path = string(@__DIR__, "/")

function load_description(;path=""::String, order=2::Int32, verbose=false::Bool)
        if path == ""
                path = default_path
        end

        if order==2
                desc_D = sbp.Operator.read_operator(string(path, "D_21.txt"))
                desc_H = sbp.Operator.read_operator(string(path, "H_21.txt"))
        elseif order==4
                desc_D = sbp.Operator.read_operator(string(path, "D_42.txt"))
                desc_H = sbp.Operator.read_operator(string(path, "H_42.txt"))
        else
           error("No operators implemented for order = ", order)
        end

        return desc_D, desc_H
end

function build_operators(desc_D::OperatorData, desc_H::OperatorData, m::Int64;
                                               is_sparse=true::Bool)
        x, h = sbp.Grid.grid_xp(m)

        D = to_matrix(desc_D, m, m, is_sparse) / h
        H = to_matrix(desc_H, m, m, is_sparse) * h

        return x, h, D, H

end


end
