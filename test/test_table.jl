using Test
import sbp
#using sbp.Table: latex

header = ["a", "b"]
a = zeros(4)
b = zeros(4)
b[1] = Inf
vars = [a, b]

@testset "Latex table" begin
table = sbp.Table.latex(header, vars, format="c")
ans = """
\\begin{tabular}{cc}
a & b \\\\ 
\\hline
\$0.0\$ & - \\\\ 
\$0.0\$ & \$0.0\$ \\\\ 
\$0.0\$ & \$0.0\$ \\\\ 
\$0.0\$ & \$0.0\$
\\end{tabular}
"""
@test join(table) == ans

end
