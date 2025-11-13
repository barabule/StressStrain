using Test
using DelimitedFiles, StaticArrays, DataInterpolations, Optim, LinearAlgebra
using RegularizationTools

using StressStrain
SS = StressStrain

@testset "Test Set 1" begin
    include("hardening_tests.jl")
end