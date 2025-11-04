


function toein_compensate(ss;
                cut = 0.0,
                )

    @assert haskey(ss, :strain) && haskey(ss, :stress) "ss must have fields strain and stress !"

    icut = findfirst(s -> s>cut, ss.strain)

    icut == lastindex(ss.strain) && return ss #fail
    #this works only for smooth data...
    next_slope = (ss.stress[icut+1] - ss.stress[icut]) / 
                 (ss.strain[icut+1] - ss.strain[icut])

    el_strain = ss.stress[icut] / next_slope
    offset_strain = ss.strain[icut] - el_strain
    
    strainout = ss.strain[icut:end] .- offset_strain
    stressout = ss.stress[icut:end]
    #add 0,0 as first point, only if there's no 0,0 point...
    if el_strain >0 
        z1 = zero(eltype(ss.strain))
        z2 = zero(eltype(ss.stress))
        prepend!(strainout, z1)
        prepend!(stressout, z2)
    end
    return (;strain = strainout,
             stress = stressout)
end

function moving_average(A::AbstractArray, m::Int = 3)
    @assert m > 0 && isodd(m) "m must be greater than 0 and odd"
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = m÷2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    return out
end