

abstract type AbstractStressStrainSignal end


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
    stress::Vector{T}
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


function true_to_engineering(SS::TrueStressStrain)
    strain = exp.(SS.strain) .- 1
    stress = SS.stress ./ (1 .+ strain)
    return EngineeringStressStrain(strain, stress)
end

####interface for AbstractArray

function length(SS::AbstractStressStrainSignal)
    SS.length
end

function size(SS::AbstractStressStrainSignal)
    (SS.length)
end


function getindex(SS::AbstractStressStrainSignal, i)
    @assert 1<=i<=SS.length "Error index must be between 1 and $(SS.length)!"
    (SS.strain[i], SS.stress[i])
end

firstindex(SS::AbstractStressStrainSignal) = (SS.strain[1], SS.stress[1])
lastindex(SS::AbstractStressStrainSignal) = (SS.strain[SS.length], SS.stress[SS.length])




function modulus(SS::AbstractStressStrainSignal; max_strain = 1e-3, sigdigits = 2)
    idx = findall(s -> s<= max_strain, SS.strain)
    if !isempty(idx)
        E = round(SS.stress[idx] \ SS.strain[idx]; sigdigits)
        return E
    end
    nothing
end