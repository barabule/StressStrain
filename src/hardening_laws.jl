

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
    return @. K * t ^ n
end