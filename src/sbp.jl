module sbp

include("config.jl")
include("grid.jl")
include("operator.jl")
include("metrics.jl")
include("sparse.jl")
include("vtk.jl")
include("test.jl")
include("source.jl")
include("table.jl")
include("lmv.jl")
include("mms.jl")
include("boundary_operator.jl")

# Time integrators
include("time_integrators/lsrk4.jl")
include("time_integrators/staggered_rk4.jl")

# Operators
include("../operators/oreilly_petersson_2019/operators.jl")
include("../operators/strand_1994/operators.jl")

# Staggered operators and discretizations
include("staggered.jl")
include("staggered/acoustic.jl")

# Collocated operators and discretizations
include("collocated.jl")
include("collocated/acoustic.jl")

end
