module sbp

include("config.jl")
include("grid.jl")
include("operator.jl")
include("metrics.jl")
include("sparse.jl")
include("vtk.jl")
include("test.jl")

include("staggered.jl")

include("staggered/acoustic.jl")


end
