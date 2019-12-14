using sbp
using Test

include("test_grid.jl")
include("test_operator.jl")
include("test_sparse.jl")
include("../operators/oreilly_petersson_2019/test_operators.jl")
include("../operators/strand_1994/test_operators.jl")
include("test_staggered.jl")
include("test_collocated.jl")
include("test_metrics.jl")
include("test_vtk.jl")
include("test_source.jl")
include("test_staggered_acoustic.jl")
include("test_lsrk4.jl")
include("test_lmv.jl")
