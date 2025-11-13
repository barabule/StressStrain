using Test
using DelimitedFiles, StaticArrays, DataInterpolations, Optim, LinearAlgebra
using RegularizationTools

@test "Test Set 1" begin
    include("hardening_tests.jl")
end