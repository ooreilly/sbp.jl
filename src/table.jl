# Module for creating tables
module Table

"""
Generate table typeset in LaTeX 

A table that looks like

        a | b |
        --|---|
        0 | - |
        0 |   |

is created by passing header = ["a", "b"], a = [0, 0], b = [Inf, 0], 
variables = [a, b] to the function. The generated table is returned as a string.
"""
function latex(header, variables; format="c")
        @assert length(header) == length(variables)
        n = length(header)
        fmt = repeat(format, n)
        table = []
        txt = "\\begin{tabular}{$fmt}\n"
        append!(table, txt)
        m = length(variables[1]) 
        for j in 1:n
                @assert length(variables[j]) == m
                hi = header[j]
                if j != n
                        txt = "$hi & "
                else
                        txt = "$hi \\\\ \n"
                end
                append!(table, txt)
        end
        txt = "\\hline\n"
        append!(table, txt)

        for i in 1:m
                for j in 1:n
                        var = variables[j][i]
                        if var != Inf
                                val = "\$$var\$"
                        else
                                val = "-"
                        end
                        if j != n
                                txt = "$val & "
                        else
                                txt = "$val"
                        end
                        append!(table, txt)
                end

                if i != m
                        txt = " \\\\ \n"
                else
                        txt = "\n"
                end
                append!(table, txt)
        end
        
        txt = "\\end{tabular}\n"
        append!(table, txt)
        return table

end


end
