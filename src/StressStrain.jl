module StressStrain

using DelimitedFiles, StaticArrays, DataInterpolations, Optim, LinearAlgebra
using CSV, DataFrames
using RegularizationTools
using GLMakie
using Observables
import GLMakie.GLFW


include("data_prep.jl")
include("hardening_laws.jl")
include("bezier.jl")
include("data_gui.jl")
include("gui.jl")

export launch_gui
export engineering_to_true, true_to_engineering
export get_modulus
export get_hardening_portion
export toein_compensate
export read_stress_strain_data

end # module StressStrain
