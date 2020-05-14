module NiLogLikeNumbers
using TropicalNumbers
using LogarithmicNumbers
using NiLang, NiLang.AD

export Tropical

include("base.jl")
include("instructs.jl")
include("blas.jl")
include("einsum.jl")

end # module
