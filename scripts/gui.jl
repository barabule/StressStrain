import Pkg; Pkg.activate("scripts")
using .StressStrain
using GLMakie
using DataInterpolations
using DelimitedFiles


function make_gui()

    fig = Figure()
    strain  = collect(LinRange(0, 1, 10))
    SSE = (;strain = strain, 
                    stress  = 100. * strain .^ 0.2)

    axss = Axis(fig[1,1], title = "Stress Strain",
                    xlabel = "Strain [-]",
                    ylabel = "Stress [MPa]")


    

    

    fitfuncs = [LinearInterpolation, CubicSpline]

    interpolant = fitfuncs[1](SSE.stress, SSE.strain)
    
    plot_stress!(axss, SSE, interpolant; N=1000)


    fit_menu = Menu(fig, options = zip(["Linear", "CSplines"], fitfuncs),
                        default = "Linear")

    

    fig[1,2] = vgrid!(
            Label(fig, "Fitting Function", width = nothing),
            fit_menu;
            tellheight = false, width = 200,
    )


    on(fit_menu.selection) do s
        #fit the function on the data
        interpolant = s(SSE.stress ,SSE.strain)
        plot_stress!(axss, SSE, interpolant)
    end

    on(events(fig).dropped_files) do files
        isempty(files) && return nothing
        
        f1 = first(files)
        println(f1)
        data = try_open(f1)
        if isnothing(data)
            println("Data is nothing!")
        else
            plot_stress!(axss, SSE, interpolant)
        end
    end



    return fig

end


function try_open(fn::AbstractString,
                    skipstart = 0)

    delims = [',','\t',';']

    for delim in delims
        rawdata = readdlm(fn, delim;skipstart)
        size(rawdata, 2) >= 2 || continue
        isa(rawdata[1, 1], Number) && isa(rawdata[1,2], Number) && continue
        return (;strain = rawdata[:,1], 
                stress = rawdata[:, 2])
    end
    return nothing
end

function plot_stress!(ax, data, interpolant; N = 100)

    empty!(ax)
    scatter!(ax, data.strain, data.stress, color = :black, marker = :cross)
    t = LinRange(extrema(data.strain)..., N)
    lines!(ax, t , interpolant.(t), color = :black)
    return nothing
end