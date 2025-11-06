



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
    #we reached the end and no point was found
end

function toein_compensate(ss;
                cut = 0.0,
                )

    @assert haskey(ss, :strain) && haskey(ss, :stress) "ss must have fields strain and stress !"

    icut = findfirst(s -> s>cut, ss.strain)

    icut == lastindex(ss.strain) || isnothing(icut) && return ss #fail
    #this works only for smooth data...
    next_slope = (ss.stress[icut+1] - ss.stress[icut]) / 
                 (ss.strain[icut+1] - ss.strain[icut])
    next_slope = clamp(next_slope, zero(next_slope), Inf) #no negative slope allowed

    el_strain = ss.stress[icut] / next_slope
    offset_strain = ss.strain[icut] - el_strain
    
    strainout = ss.strain[icut:end] .- offset_strain
    stressout = ss.stress[icut:end]
    #add 0,0 as first point, only if there's no 0,0 point...
    if offset_strain > 0 
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



function decimate(xp, yp; rescale = true, tolerance=1e-3)
    if rescale #rescale x and y to be within [0,1]
        x , y = copy(xp), copy(yp)
        
        xmin, xmax = extrema(x)
        ymin, ymax = extrema(y)

        x .-= xmin
        x .*= 1/(xmax - xmin)
        y .-= ymin
        y .*= 1/(ymax - ymin)
    else
        x, y = xp, yp
    end

    resultlist = [firstindex(x)]
    DouglasPeucker!(resultlist, (x, y), firstindex(x), lastindex(x); tolerance)
    idx =  sort(resultlist)

    xout = xp[idx]
    yout = yp[idx]

    return (xout, yout)
end


function DouglasPeucker!(resultlist, points, idx1, idx2; tolerance = 1e-3)
    x ,y = points
    dmax = zero(eltype(x))
    index = idx1

    for i in idx1+1 : idx2
        d = distance((x[i], y[i]), (x[idx1], y[idx1]), (x[idx2], y[idx2]))
        if d > dmax
            dmax = d
            index = i
        end
    end

    if dmax>tolerance
        DouglasPeucker!(resultlist, points, idx1, index; tolerance)
        DouglasPeucker!(resultlist, points, index, idx2; tolerance)
    else
        push!(resultlist, idx2)
    end
    
end

function distance(P1, L1, L2)
    a = P1 - L1
    d = L2 - L1
    nd = norm(d)

    return norm(a .- dot(d, a) .* d ./ norm(d))
end


function resample_curve(x, y, N::Integer=100;
                            resampler = LinearInterpolation, 
                            tolerance = 1e-3, #tolerance for decimate
                            d = 3, #degree for BSplineApprox or RegularizationSmooth
                            h = 4, #number of control pts for BSplineApprox
                            )

    @assert length(x) == length(y) "x and y must  have the same length"

    
    xi = collect(LinRange(extrema(x)..., N))
    if any(resampler .== (LinearInterpolation, CubicSpline, AkimaInterpolation, QuadraticSpline))
        interpolator = resampler(y, x)
        return (xi, interpolator.(xi))
    elseif resampler == decimate #first ecimate, then interpolate
        xii, yii = decimate(x, y; tolerance)
        return resample(xii, yii, N;resampler = LinearInterpolation)
    elseif resampler == BSplineApprox
        interpolator = BSplineApprox(y, x, d, h, :ArcLen, :Average)
        return (xi, interpolator.(xi))
    elseif resampler == RegularizationSmooth
        interpolator = RegularizationSmooth(y, x, d; alg = :gcv_svd)
        return (xi, interpolator.(xi))
    
    end

    return nothing ###

end
