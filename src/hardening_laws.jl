

### interface f(t, p) where p is a tuple of parameters


function Swift(t, p)
    @assert length(p) >= 3 "p must have at least 3 elements!"
    K, ϵ0, n = p[1:3]
    return  @. K * abs(ϵ0 + t)^n

end


function Voce(t, p)
    @assert length(p) >= 3 "p must have at least 3 elements!"
    σ0, Rsat, ζ = p[1:3]

    return @. σ0 + Rsat * (1 - exp(-ζ * t))
end

function HockettSherby(t, p)
    @assert length(p) >= 4 "p must have at least 4 elements!"
    A, B, C, H = p[1:4]
    return @. A - B * exp(-C * abs(t)^H)
end

function StoughtonYoon(t, p)
    @assert length(p) >= 5 "p must have at least 5 elements!"
    A, B, C, m, D = p[1:5]
    return @. A - B * exp(-C * abs(t)^m) + D * t
end

function Hollomon(t ,p)
    @assert length(p)>=2 "p must have at least 2 elements!"
    K, n = p[1:2]
    return @. K * abs(t) ^ n
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
    end
    return nothing
end