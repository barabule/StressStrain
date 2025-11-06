module StressStrain

using DelimitedFiles, StaticArrays, DataInterpolations, Optim, LinearAlgebra
using GLMakie


include("data_prep.jl")
include("hardening_laws.jl")
include("gui.jl")


end # module StressStrain
