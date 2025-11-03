module StressStrain

using DelimitedFiles, StaticArrays, DataInterpolations, Optim, LinearAlgebra


include("typedefs.jl")
include("data_prep.jl")
include("hardening_laws.jl")

end # module StressStrain
