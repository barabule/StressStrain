

abstract type AbstractStressStrainSignal<:AbstractVector end


struct EngineeringStressStrain{T<:Real} <: AbstractStressStrainSignal
    strain::Vector{T}
    stress::Vector{T}
    length::Integer
    function EngineeringStressStrain(strain::AbstractVector, stress::AbstractVector)
        T1 = eltype(strain)
        T2 = eltype(stress)
        T = promote_type(T1, T2)
        N1= length(strain)
        N2 = length(stress)
        @assert N1 == N2 "Strain and Stress must have the same length!"
        return new{T}(T.(strain), T.(stress), N1)
    end
end

struct TrueStressStrain{T<:Real} <:AbstractStressStrainSignal
    strain::Vector{T}
    stress:Vector{T}
    length::Integer
    function TruetressStrain(strain::AbstractVector, stress::AbstractVector)
        T1 = eltype(strain)
        T2 = eltype(stress)
        T = promote_type(T1, T2)
        N1= length(strain)
        N2 = length(stress)
        @assert N1 == N2 "Strain and Stress must have the same length!"
        return new{T}(T.(strain), T.(stress), N1)
    end
end


function engineering_to_true(SS::EngineeringStressStrain)
    logstrain = log.(1 .+ SS.strain)
    truestress = SS.stress .* (1 .+ SS.strain)
    return TrueStressStrain(logstrain, truestress)
end


function true_to_engineering(SS:TrueStressStrain)
    strain = exp.(SS.strain) .- 1
    stress = SS.stress ./ (1 .+ strain)
    return EngineeringStressStrain(strain, stress)
end

import Base:getindex, length, size

function length(SS::AbstractStressStrainSignal)
    SS.length
end

function size(SS::AbstractStressStrainSignal)
    (SS.length)
end


function getindex(SS::AbstractStressStrainSignal, i)
    (SS.strain[i], SS.stress[i])
end



function nose(SS::AbstractStressStrainSignal)


end


function modulus(SS::AbstractStressStrainSignal; max_strain = 1e-3)
    idx = findall(s -> s<= max_strain, SS.strain)
    if !isempty(idx)
        return SS.stress[idx] ./ SS.strain[idx]
    end
    nothing
end