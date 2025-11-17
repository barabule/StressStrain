module StressStrain

using DelimitedFiles, StaticArrays, DataInterpolations, Optim, LinearAlgebra
using RegularizationTools
using GLMakie
import GLMakie.GLFW


include("data_prep.jl")
include("hardening_laws.jl")
include("bezier.jl")
include("data_gui.jl")
include("gui.jl")

export main
export engineering_to_true, true_to_engineering
export get_modulus
export get_hardening_portion
export toein_compensate
export read_stress_strain_data

end # module StressStrain
