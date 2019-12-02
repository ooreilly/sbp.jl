module sbp

include("config.jl")
include("grid.jl")
include("operator.jl")
include("metrics.jl")
include("sparse.jl")

include("staggered.jl")
include("test.jl")

include("staggered/acoustic.jl")


end # module
