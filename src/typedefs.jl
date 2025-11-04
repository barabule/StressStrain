


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
    @assert haskey(SS, :strain) && haskey(SS, :stress) "SS must have fields strain and stress !"
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