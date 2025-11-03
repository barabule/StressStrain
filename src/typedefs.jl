

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




function get_modulus(SS::AbstractStressStrainSignal; max_strain = 1e-3, sigdigits = 2)
    
    N = SS.length
    N >= 2 || error("Stress Strain data needs to have at least 2 points!")
    idx = findall(s -> s<= max_strain, SS.strain)
    if !isempty(idx)
        E = round(SS.stress[idx] \ SS.strain[idx]; sigdigits)
        return E
    end
    E = (SS.stress[2] - SS.stress[1]) / ( SS.strain[2] - SS.strain[1]) #fallback
end

function get_hardening_portion(SS::TrueStressStrain, modulus = nothing;offset = 2e-3)
    if isnothing(modulus)
        E = get_modulus(SS)
    else
        E = modulus
    end

    for i in 1:SS.length
        ϵ, σ = SS.strain[i], SS.stress[i]
        ϵel = σ / E
        ϵpl = ϵ - ϵel
        if ϵpl > offset
            return TrueStressStrain(SS.strain[i:end], SS.stress[i:end])
        end
        
    end
    return nothing
end