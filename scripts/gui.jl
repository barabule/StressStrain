import Pkg; Pkg.activate("scripts")
using StressStrain
SS = StressStrain
using GLMakie
using DataInterpolations
using DelimitedFiles
using Optim

function make_gui(; N=1000)

    fig = Figure()
    strain  = collect(LinRange(0, 1, 50))
    SSE = (;strain = strain, 
                    stress  = 100. * strain .^ 0.2)

    axss = Axis(fig[1,1], title = "Stress Strain",
                    xlabel = "Strain [-]",
                    ylabel = "Stress [MPa]")


    

    

    fitfuncs = [LinearInterpolation, 
                CubicSpline,
                BSplineApprox, 
                SS.Swift, 
                SS.Voce, 
                SS.HockettSherby, 
                SS.StoughtonYoon]

    interpolant = fitfuncs[1](SSE.stress, SSE.strain)
    
    plot_stress!(axss, SSE, interpolant; N=1000)


    fit_menu = Menu(fig, options = zip(["Linear", 
                                        "CSplines", 
                                        "BSplineApprox",
                                        "Swift", 
                                        "Voce", 
                                        "HockettSherby",
                                        "StoughtonYoon"], fitfuncs),
                        default = "Linear")

    

    fig[1,2] = vgrid!(
            Label(fig, "Fitting Function", width = nothing),
            fit_menu;
            tellheight = false, width = 200,
    )


    on(fit_menu.selection) do s
        #fit the function on the data

        # interpolant = s(SSE.stress ,SSE.strain)
        interpolant = get_interpolant(s, SSE)
        plot_stress!(axss, SSE, interpolant; N)
    end

    on(events(fig).dropped_files) do files
        isempty(files) && return nothing
        
        f1 = first(files)
        println(f1)
        SSE = try_open(f1)
        
        plot_stress!(axss, SSE, interpolant;N)
        # if isnothing(data)
        #     println("Data is nothing!")
        # else
        #     plot_stress!(axss, SSE, interpolant)
        # end
    end



    return fig

end


function try_open(fn::AbstractString,
                    skipstart = 0)


    rawdata = readdlm(fn)
    return (;strain = rawdata[:,1],
            stress = rawdata[:, 2])
    # delims = [',','\t',';']

    # for delim in delims
    #     rawdata = readdlm(fn, delim;skipstart)
    #     size(rawdata, 2) >= 2 || continue
    #     isa(rawdata[1, 1], Number) && isa(rawdata[1,2], Number) && continue
    #     return (;strain = rawdata[:,1], 
    #             stress = rawdata[:, 2])
    # end
    # return nothing
end

function plot_stress!(ax, data, interpolant; N = 100)

    empty!(ax)
    scatter!(ax, data.strain, data.stress, color = :black, marker = :cross)
    t = LinRange(extrema(data.strain)..., N)
    lines!(ax, t , interpolant.(t), color = :black)
    return nothing
end


function get_interpolant(func, data)
    if func == SS.Swift || func == SS.Voce
        p0 = [100.0, 1e-3, 0.2]
        lb = zeros(3)
        ub = [Inf, Inf, Inf]
        return Curvefit(data.stress, data.strain, func, p0, NelderMead(), true, lb, ub)
    elseif func == SS.HockettSherby
        p0 = [50.0, 100.0, 1.0, 0.3]
        lb = zeros(4)
        ub = fill(Inf, 4)
        return Curvefit(data.stress, data.strain, func, p0, NelderMead(), true, lb, ub)
    elseif func == SS.StoughtonYoon
        p0 = [50.0, 100.0, 1.0, 1.0, 1e-5]
        lb = zeros(5)
        ub = fill(Inf, 5)
        return Curvefit(data.stress, data.strain, func, p0, NelderMead(), true, lb, ub)
    elseif func == LinearInterpolation || func == CubicSpline
        return func(data.stress, data.strain)
    elseif func == BSplineApprox
        return func(data.stress, data.strain, 3, 4, :ArcLen, :Average)
    end
    return nothing
end