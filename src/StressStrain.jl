module StressStrain

using DelimitedFiles, StaticArrays, DataInterpolations, Optim, LinearAlgebra
using RegularizationTools
using GLMakie
import GLMakie.GLFW


include("data_prep.jl")
include("hardening_laws.jl")
include("data_gui.jl")
include("gui.jl")


end # module StressStrain
