

# abstract type AbstractStressStrainSignal end


# struct EngineeringStressStrain{T<:Real, I1<:Integer} <: AbstractStressStrainSignal
#     strain::Vector{T}
#     stress::Vector{T}
#     length::Integer
#     function EngineeringStressStrain(strain::AbstractVector, stress::AbstractVector, L)
#         T1 = eltype(strain)
#         T2 = eltype(stress)
#         T = promote_type(T1, T2)
#         I1 = typeof(L)
#         N1 = length(strain)
#         N2 = length(stress)
#         @assert N1 == N2 == L "Strain and Stress must have the same length!"
#         return new{T, I1}(T.(strain), T.(stress), L)
#     end
# end

# struct TrueStressStrain{T<:Real, I1<:Integer} <:AbstractStressStrainSignal
#     strain::Vector{T}
#     stress::Vector{T}
#     length::Integer
#     function TruetressStrain(strain::AbstractVector, stress::AbstractVector, L)
#         T1 = eltype(strain)
#         T2 = eltype(stress)
#         T = promote_type(T1, T2)
#         I1 = typeof(L)
#         N1 = length(strain)
#         N2 = length(stress)
#         @assert N1 == N2 ==L "Strain and Stress must have the same length!"
#         return new{T, I1}(T.(strain), T.(stress), L)
#     end
# end

# function TrueStressStrain(v1::Vector{T}, v2::Vector{T}) where T
#     N = length(v1)
#     @assert N == length(v2) "Lengths of v1 and v2 must be the same"
#     return TrueStressStrain(v1, v2, N)
# end

# function EngineeringStressStrain(v1::Vector{T}, v2::Vector{T}) where T
#     N = length(v1)
#     @assert N == length(v2) "Lengths of v1 and v2 must be the same"
#     return EngineeringStressStrain(v1, v2, N)
# end


function engineering_to_true(SS)
    logstrain = log.(1 .+ SS.strain)
    truestress = SS.stress .* (1 .+ SS.strain)
    
    return (;strain = logstrain, 
             stress = truestress)
end


function true_to_engineering(SS)
    strain = exp.(SS.strain) .- 1
    stress = SS.stress ./ (1 .+ strain)
    
    return (;strain, stress)
end




function get_modulus(SS; max_strain = 1e-3, sigdigits = 2)
    
    N = length(SS.strain)
    N >= 2 || error("Stress Strain data needs to have at least 2 points!")
    idx = findall(s -> s<= max_strain, SS.strain)
    if !isempty(idx) && length(idx) >= 2
        E = SS.strain[idx] \ SS.stress[idx]
    else
        E = (SS.stress[2] - SS.stress[1]) / ( SS.strain[2] - SS.strain[1]) #fallback
    end

    return round(E; sigdigits)
end

function get_hardening_portion(SS, modulus = nothing;offset = 2e-3)
    if isnothing(modulus)
        E = get_modulus(SS)
    else
        E = modulus
    end
    @assert haskey(SS, :strain) && haskey(SS, :stress) "SS must have fields strain and stress !"
    for i in 1:length(SS.strain)
        ϵ, σ = SS.strain[i], SS.stress[i]
        ϵel = σ / E
        ϵpl = ϵ - ϵel
        if ϵpl > offset
            return (;strain = SS.strain[i:end] .- ϵel,
                     stress = SS.stress[i:end])
        end
        
    end
    return nothing
end