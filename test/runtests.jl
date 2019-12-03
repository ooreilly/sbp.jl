using sbp
using Test

include("test_grid.jl")
include("test_operator.jl")
include("test_sparse.jl")
include("../operators/oreilly_petersson_2019/test_operators.jl")
include("test_staggered.jl")
include("test_metrics.jl")
include("test_vtk.jl")
include("test_staggered_acoustic.jl")
