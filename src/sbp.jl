module sbp

include("config.jl")
include("grid.jl")
include("operator.jl")
include("metrics.jl")
include("sparse.jl")
include("vtk.jl")
include("test.jl")
include("source.jl")
include("lmv.jl")

# Time integrators
include("time_integrators/lsrk4.jl")

# Staggered operators and discretizations
include("staggered.jl")
include("staggered/acoustic.jl")


end
