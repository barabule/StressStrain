



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


function get_modulus(SS; 
                max_strain = 1e-3, #consider strains up to max_strain for fitting
                sigdigits = 2, #makes no sense to return E moduli with too much precision
                ) 



    @assert haskey(SS, :strain) && haskey(SS, :stress) "SS must have fields strain and stress !"
    
    N = length(SS.strain)
    @assert N >= 2 "Stress Strain data needs to have at least 2 points!"

    idx = findall(s -> 0 <= s <= max_strain, SS.strain)
    if !isempty(idx) && length(idx) >= 2
        E = SS.strain[idx] \ SS.stress[idx]
    else
        E = (SS.stress[2] - SS.stress[1]) / ( SS.strain[2] - SS.strain[1]) #fallback
    end

    Emin, Emax = bracket_modulus(SS)
    E = clamp(E, Emin, Emax)
    
    return round(E; sigdigits)
end


function bracket_modulus(data;
                    sigdigits = 2)
    @assert haskey(data, :strain) && haskey(data, :stress) "Data must have fields strain and stress!"

    Emin , Emax = Inf, -Inf

    for i in eachindex(data.strain)
        ϵ, σ = data.strain[i], data.stress[i]
        # ϵ <= 0 || σ <= 0 && continue #ignore invalid data
        slp = σ / ϵ #secant modulus
        slp <=0 && continue 
        Emin = slp < Emin ? slp : Emin
        Emax = slp > Emax ? slp : Emax
    end
    @assert Emin > 0 && Emax > 0
    return (;Emin = round(Emin; sigdigits), 
             Emax = round(Emax; sigdigits))
end


function get_hardening_portion(SS, modulus = nothing;offset = 2e-3)
    if isnothing(modulus)
        E = get_modulus(SS)
    else
        E = modulus
    end
    @assert haskey(SS, :strain) && haskey(SS, :stress) "SS must have fields strain and stress !"

    plastic_strain = Vector{eltype(SS.strain)}()
    effective_stress = Vector{eltype(SS.stress)}()

    
    for i in eachindex(SS.strain)
        ϵ, σ = SS.strain[i], SS.stress[i]
        ϵel = σ / E
        ϵpl = ϵ - ϵel
        if ϵpl > offset
            push!(plastic_strain, ϵpl)
            push!(effective_stress, σ)
        end    
    end
    #this is sometimes needed for noisy data
    sorted_idx = sortperm(plastic_strain)

    return (;strain =  plastic_strain[sorted_idx],
            stress = effective_stress[sorted_idx])
end


function toein_compensate(ss;
                cut = 0.0, #cut at this strain
                elastic_strain_offset = 1e-3, #how much of the new curve to use to get a slope
                min_slope = -Inf, #smallest allowed slope 
                )

    @assert haskey(ss, :strain) && haskey(ss, :stress) "ss must have fields strain and stress !"
    
    min_slope = max(min_slope, last(ss.stress)/last(ss.strain)) #E modulus cannot be < than last secant modulus

    cut ≈ 0 && return ss #no need to do anything


    icut = findfirst(s -> s >= cut, ss.strain)

    icut == lastindex(ss.strain) || isnothing(icut) && return ss #fail
    
    #get the elastic portion of the cut curve
    elastic_strain_offset = clamp(elastic_strain_offset, 0, Inf)
    islope = findlast( e -> e <= cut + elastic_strain_offset, ss.strain)
       
    if elastic_strain_offset ≈ 0 || (islope <= icut) #just get the next slope
        next_slope = (ss.stress[icut+1] - ss.stress[icut]) / 
                    (ss.strain[icut+1] - ss.strain[icut])
    else
        next_slope = get_modulus((;strain = view(ss.strain, icut:islope), 
                                    stress = view(ss.stress, icut:islope));
                                    max_strain = ss.stress[islope])
    end

    next_slope = clamp(next_slope, min_slope, Inf) #no negative or zero slope allowed    
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
                            navg = 3, #number of points for moving average filter
                            offset = 2e-3,
                            E = 1.0,
                            alg = NelderMead(),
                            )

    @assert length(x) == length(y) "x and y must  have the same length"

    
    xi = collect(LinRange(extrema(x)..., N))
    if any(resampler .== (LinearInterpolation, CubicSpline, AkimaInterpolation, 
                QuadraticSpline, QuadraticInterpolation))
        interpolator = resampler(y, x)
        return (xi, interpolator.(xi))
    elseif resampler == decimate #first ecimate, then interpolate
        xii, yii = decimate(x, y; tolerance)
        return resample(xii, yii, N;resampler = LinearInterpolation)
    elseif resampler == BSplineApprox
        interpolator = BSplineApprox(y, x, d, h, :ArcLen, :Average)
        return (xi, interpolator.(xi))
    elseif resampler == RegularizationSmooth
        xi, yi = resample_curve(x, y, N; resampler = LinearInterpolation)
        interpolator = RegularizationSmooth(yi, xi, xi, d; alg = :gcv_svd)
        return (xi, interpolator.(xi))
    elseif resampler == moving_average #first smooth then resample
        newy = moving_average(y, navg)
        return resample_curve(x, newy, N; resampler = LinearInterpolation)
    elseif resampler == Hollomon
        p0 = [100.0, 0.2]
        lb = zeros(2)
        ub = [Inf, Inf]
        interpolator =  Curvefit(y, x, resampler, p0, NelderMead(), true, lb, ub; extrapolate = true)
        return (xi, interpolator(xi))
    elseif resampler == RamberOsgoodAlternativeReparametrized
        # σy, ϵy, b, r
        sy = 0.5 * maximum(y)
        ey = 0.1 * maximum(x) + 0.9 * minimum(x)
        b = 0.2
        r = 1.0

        p0 = [sy, ey, b, r]
        lb = [1e-8, 1e-8, 0.0, 0.0]
        ub = [Inf, Inf, Inf, Inf]
        
        interpolator = Curvefit(y, x, resampler, p0, alg, true, lb, ub;extrapolate = true)
        return (xi, interpolator(xi))
    end

    return nothing ###

end


function read_stress_strain_data(fn::AbstractString;
                        delim='\t',
                        skipstart = 0,
                        strain_multiplier = 1.0,
                        stress_multiplier = 1.0,
                        strain_col = 1,
                        stress_col = 2,
                        T = Float64,
                        clean_sort = true
                        )

    rawdata = readdlm(fn, delim; skipstart)
    folder = dirname(fn)

    @assert size(rawdata, 2) >= max(strain_col, stress_col) "Error"
    strain = T.(rawdata[:, strain_col]) .* strain_multiplier
    stress = T.(rawdata[:, stress_col]) .* stress_multiplier

    if clean_sort
        sortedidx = sortperm(strain)
        strain = strain[sortedidx]
        stress = stress[sortedidx]
    end

    return (;strain, stress, folder)
end


function cutoff(data, val)
    @assert val > 0 "Cutoff value must be positive!" 
   @assert haskey(data, :strain) && haskey(data, :stress) "ss must have fields strain and stress !"
    for i in reverse(eachindex(data.strain))
        if data.strain[i] < val
            return (;strain = data.strain[1:i],
                    stress = data.stress[1:i])
        end
    end
    return data

end


function find_abnormal_points(data; tooclose = 1e-6)

    @assert haskey(data, :strain) && haskey(data, :stress) "data must have strain and stress fields !"


    abnormal_indices = Int[]
    unique_strains = Dict{Float64, Int}()
    lastidx = -1
    for i in eachindex(data.strain)
        eps, sig = data.strain[i], data.stress[i]

        if i != firstindex(data.strain) 

            if data.strain[lastidx] > eps #non monotonic strain
                push!(abnormal_indices, i)
                lastidx = i
                continue
            end

        end
        lastidx = i

        eps_key = round(Int, eps / tooclose)
        if haskey(unique_strains, eps_key)#find repeated strain pts
            push!(abnormal_indices,  i)
            continue
        else
            push!(unique_strains, eps_key => i)
        end

        if eps < 0 || sig < 0#negative strain / stress
            push!(abnormal_indices, i)
            continue
        end
        
        
    end

    return abnormal_indices
end