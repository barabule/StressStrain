

### interface f(t, p) where p is a tuple of parameters
###none of these should work with negative t

function Swift(t, p)
    @assert length(p) >= 3 "p must have at least 3 elements!"
    @assert all(t .>= 0) "t must be positive" 
    K, ϵ0, n = p[1:3]
    return  @. K * abs(ϵ0 + t)^n
end


function Voce(t, p)
    @assert length(p) >= 3 "p must have at least 3 elements!"
    @assert all(t .>= 0) "t must be positive" 
    σ0, Rsat, ζ = p[1:3]

    return @. σ0 + Rsat * (1 - exp(-ζ * t))
end

function HockettSherby(t, p)
    @assert length(p) >= 4 "p must have at least 4 elements!"
    @assert all(t .>= 0) "t must be positive" 
    A, B, C, H = p[1:4]
    return @. A - B * exp(-C * abs(t)^H)
end

function StoughtonYoon(t, p)
    @assert length(p) >= 5 "p must have at least 5 elements!"
    @assert all(t .>= 0) "t must be positive" 
    A, B, C, m, D = p[1:5]
    return @. A - B * exp(-C * abs(t)^m) + D * t
end

function Hollomon(t ,p)
    @assert length(p)>=2 "p must have at least 2 elements!"
    @assert all(t .>= 0) "t must be positive" 
    K, n = p[1:2]
    return @. K * abs(t) ^ n
end


function Bilinear(t, p)
    @assert length(p)>=2 "p must have at least 2 elements!"
    @assert all(t .>= 0) "t must be positive" 
    sy, Etan = p[1:2]
    return @. sy + Etan * t
end


function SwiftVoce(t, p)
    @assert length(p) >= 8 "p must have at least 8 elements"
    @assert all(t .>= 0) "t must be positive" 
    w1, w2, K, ϵ0, n, σ0, Rsat, ζ = p[1:8]
    return @. w1 * (K * abs(ϵ0 + t)^n) + 
              w2 * (σ0 + Rsat * (1 - exp(-ζ * t)))
end


function RamberOsgoodAlternativeReparametrized(t, p)
    @assert length(p)>=4 "p must have at least 4 elements!"
    σy, ϵy, b, r = p[1:4]

    ϵbar = @. t / ϵy
    σbar = @. b * ϵbar + (1 - b) * ϵbar / (1 + abs(ϵbar) ^ r) ^ (1 / r)
    return @. σbar * σy   
end


function RambergOsgood(t, p; offset = 2e-3, E = 1.0, maxiter = 50)
    @assert length(p) >= 2 "p must have at least 2 elements!"
    n, sigy = p[1:2]

    alpha = offset * E/sigy
    return @. find_RO_stress(t, E, sigy, n, alpha; maxiter)
end

function find_RO_stress(t, E, sigy, n, alpha; 
        maxiter = 20, 
        abstol = 1e-8, #abs error tolerance
        dtol = 1e-12, #derivative tolerance
        ctol = 1e-4, #convergence tolerance
        sig0 = E * t, #first guess
        ) 
    
    @assert t>=0 "Strain must be positive!"
    @assert E > 0 "E modulus must be positive!"
    @assert sigy >0 "Yield stress must be positive!"
    

    sig = sig0
    for it in 1:maxiter
        fval, dval = RO_with_derivative(sig, E, sigy, n, alpha)
        if abs(fval) < abstol
            return sig
        end
        if abs(dval) < dtol
            @warn "Derivate almost 0; cannot converge!"
            return sig
        end

        sig_new = sig  - fval / dval

        relative_change = abs(sig_new - sig) / sig
        if relative_change < ctol #convergence
            return sig_new
        end
        sig = sig_new
    end
    @warn "Maximum iterations ($max_iter) reached!"
    return sig 
end


function RO_with_derivative(sigma, E, sigy, n, alpha)
    elastic_strain = sigma / E
    plastic_strain = alpha * (sigma / E) * (sigma / sigy)^n
    total_strain = elastic_strain + plastic_strain

    d_elastic = 1/E
    d_plastic = (alpha * n / (E * sigy ^ (n -1))) * sigma ^ (n-1)
    d_total = d_elastic + d_plastic

    return (total_strain, d_total)
end



function make_interpolant(func, data; alg = NelderMead())
    @assert haskey(data, :strain) && haskey(data, :stress) "data must have strain and stress fields !"
    if func == Swift || func == Voce
        p0 = [100.0, 1e-3, 0.2]
        lb = zeros(3)
        ub = [Inf, Inf, Inf]
        return Curvefit(data.stress, data.strain, func, p0, alg, true, lb, ub; extrapolate = true)
    elseif func == HockettSherby
        p0 = [50.0, 100.0, 1.0, 0.3]
        lb = zeros(4)
        ub = fill(Inf, 4)
        return Curvefit(data.stress, data.strain, func, p0, alg, true, lb, ub; extrapolate = true)
    elseif func == StoughtonYoon
        p0 = [50.0, 100.0, 1.0, 1.0, 1e-5]
        lb = zeros(5)
        ub = fill(Inf, 5)
        return Curvefit(data.stress, data.strain, func, p0, alg, true, lb, ub; extrapolate = true)
    elseif func == LinearInterpolation || func == CubicSpline || func == PCHIPInterpolation
        return func(data.stress, data.strain; extrapolation = ExtrapolationType.Linear)
    elseif func == BSplineApprox
        return func(data.stress, data.strain, 3, 4, :ArcLen, :Average; extrapolation = ExtrapolationType.Linear)
    elseif func == Bilinear
        p0 = [30.0, 1000.0]
        lb = zeros(2)
        ub=  [Inf, Inf]
        return Curvefit(data.stress, data.strain, func, p0, alg, true, lb, ub; extrapolate = true)
    elseif func == SwiftVoce
        w1, w2 = 0.5, 0.5
        psw = [100, 1e-3, 0.2]
        pvc = [50, 100, 0.3]
        p0 = [w1, w2, psw..., pvc...]
        lb = [-Inf, -Inf, zeros(6)...]
        ub= fill(Inf, 8)
        return Curvefit(data.stress, data.strain, func, p0, alg, true, lb, ub; extrapolate = true)
    end
    return nothing
end


function interpolant_label(interpolant, func; sigdigits = 3)
    p = zeros(10)
    try
        p = round.(interpolant.pmin; sigdigits)
    catch
        # @info "Didn't work", func
        sy = round(interpolant(0); sigdigits)
        return "σy = $sy, No parameters"
    end
    sy = round(func(0, p); sigdigits)
    if func == Swift
        return "σy = $sy, Swift: K = $(p[1]), ϵ0 = $(p[2]), n = $(p[3])"
    elseif func == Voce
        return "σy = $sy, Voce: σ0 = $(p[1]), Rsat = $(p[2]), ζ = $(p[3])"
    elseif func == HockettSherby
        return "σy = $sy, Hocket-Sherby: A = $(p[1]), B = $(p[2]), C = $(p[3]), H = $(p[4])"
    elseif func == StoughtonYoon
        return "σy = $sy, Stoughton Yoon: A = $(p[1]), B = $(p[2]), C = $(p[3]), m = $(p[4]), D = $(p[5])"
    elseif func == Bilinear
        return "σy = $(p[1]), Bilinear: Etan = $(p[2])"
    elseif func == SwiftVoce
        return "σy = $sy, Swift-Voce: w1 = $(p[1]), w2 = $(p[2]), K = $(p[3]), ϵ0 = $(p[4]), n = $(p[5]), σ0 = $(p[6]), Rsat = $(p[7]), ζ = $(p[8])"
    elseif func == RamberOsgoodAlternativeReparametrized
        σy, ϵy, b, r = p[1:4]
        return "σy = $σy, Alternative Ramberg: ϵy = $ϵy, b = $b, r = $r"
    else
        return "This was triggered"
    end
    
end