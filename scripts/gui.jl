import Pkg; Pkg.activate("scripts")
using StressStrain
SS = StressStrain
using GLMakie
using DataInterpolations
using DelimitedFiles
using Optim



function make_gui(; 
                N=1000,
                sidebar_width = 300)

    fig = Figure()
    
    SSE, fitfuncs, fitfunclabels = initialize()

    axss = Axis(fig[1,1], title = "Stress Strain",
                    xlabel = "Engineering Strain [-]",
                    ylabel = "Engineering Stress [MPa]")


    plot_stress!(axss, SSE; N)
    axislegend(axss, position = :rb, merge = true)
    fit_menu = Menu(fig, options = zip(fitfunclabels, fitfuncs),
                    default = "Linear")
    

    

    # sliderobservables = [s.value for s in sg.sliders]
    gl = GridLayout(fig[1,2], width = sidebar_width, tellheight = false)

    gl_bot = GridLayout(fig[2, :], height = 100, tellwidth = false)

    subgl1 = GridLayout(gl[1,1])
    subgl2 = GridLayout(gl[2,1])

    cb_true = Checkbox(subgl1[1,2], checked = false, 
                            )
    Label(subgl1[1,1], "True", halign = :left)


    Label(subgl1[2,:], "Extrapolation Strain", width = nothing)       
    
    tb_extrapolation_strain = Textbox(subgl1[3,:], 
            validator = Float64,
            )



    Label(subgl1[5,:], "Material Name")
    tb_name = Textbox(subgl1[6,:], 
             width = sidebar_width)


    sldvals = initialize_slider(SSE)
    sld_modulus = Slider(subgl1[8, :], 
                    range = LinRange(sldvals.vmin, sldvals.vmax, 1001),
                    startvalue= sldvals.value,
                    update_while_dragging =true,
                    width = sidebar_width,
                    )

    lab_modulus = Label(subgl1[7, :], "E = $(round(sld_modulus.value[]; sigdigits= 3))MPa")

    subgl2[1,1] = vgrid!(
            Label(fig, "Fitting Function", width = nothing),
            fit_menu,
            # (Label("True", alignmode = :right), Checkbox(checked = false)),
            ;
            tellheight = false, width = 200,
    )    




    ####################EVENTS########################################################################        
    on(cb_true.checked) do val
        SSE["is true"] = val
        update_SSE!(SSE)
        plot_stress!(axss, SSE; N)
        axss.xlabel = val ? "True Strain [-]" : "Engineering Strain [-]"
        axss.ylabel = val ? "True Stress [MPa]" : "Engineering Stress [MPa]"
    end

    
    on(sld_modulus.value) do val
        update_SSE!(SSE, nothing; modulus = round(val; sigdigits =3))
        plot_stress!(axss, SSE; N)
        lab_modulus.text = "E = $(round(sld_modulus.value[]; sigdigits= 3))MPa"
    end
    




    on(fit_menu.selection) do s
        #fit the function on the data

        # interpolant = s(SSE.stress ,SSE.strain)
        SSE["interpolant"] = s
        update_SSE!(SSE)
        plot_stress!(axss, SSE; N)
    end




    on(events(fig).dropped_files) do files
        isempty(files) && return nothing
        
        f1 = first(files)
        println(f1)
        data = try_open(f1)
        
        update_SSE!(SSE, data)
        
        sldvals = initialize_slider(SSE)
        sld_modulus.range = LinRange(sldvals.vmin, sldvals.vmax, 1001)
        _ = set_close_to!(sld_modulus, sldvals.value)

        plot_stress!(axss, SSE; N)
        
    end

    on(tb_extrapolation_strain.stored_string) do s
        SSE["export max strain"] = clamp(parse(Float64, s), last(SSE["hardening"].strain), Inf)
        plot_stress!(axss, SSE; N)
    end

    on(tb_name.stored_string) do s
        SSE["name"] = s
        axss.title = s
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

function plot_stress!(ax, SSE; N = 100, tmax = SSE["export max strain"])

    empty!(ax)
    #true stress
    tss = SSE["true stress"]
    scatter!(ax, tss.strain, tss.stress, 
                    label = "True Stress (Exp)",
                    color = :black, marker = :cross, alpha = 0.3)
    lines!(ax, tss.strain, tss.stress, 
                    label = "True Stress (Exp)",
                    color = :black, linewidth = 0.5)
    t = LinRange(0, tmax, N)
    hss = SSE["hardening"]
    scatter!(ax, hss.strain, hss.stress, 
                label = "Hardening Portion (Exp)",
                color = :red, marker = :cross, alpha = 0.3)

    lines!(ax, t , SSE["hardening fit"](t), color = :black, label = "Hardening Law (Fit)", linewidth = 2.0)
    
    #plot the linear portion
    E = SSE["modulus"]
    smax = maximum(SSE["hardening"].stress)
    tmax = smax / E
    tlin = LinRange(0, tmax, 10)
    slin = E .* tlin
    lines!(ax, tlin, slin, linestyle = :dash, color = :red)


    return nothing
end


function get_interpolant(func, data)
    if func == SS.Swift || func == SS.Voce
        p0 = [100.0, 1e-3, 0.2]
        lb = zeros(3)
        ub = [Inf, Inf, Inf]
        return Curvefit(data.stress, data.strain, func, p0, NelderMead(), true, lb, ub; extrapolate = true)
    elseif func == SS.HockettSherby
        p0 = [50.0, 100.0, 1.0, 0.3]
        lb = zeros(4)
        ub = fill(Inf, 4)
        return Curvefit(data.stress, data.strain, func, p0, NelderMead(), true, lb, ub; extrapolate = true)
    elseif func == SS.StoughtonYoon
        p0 = [50.0, 100.0, 1.0, 1.0, 1e-5]
        lb = zeros(5)
        ub = fill(Inf, 5)
        return Curvefit(data.stress, data.strain, func, p0, NelderMead(), true, lb, ub; extrapolate = true)
    elseif func == LinearInterpolation || func == CubicSpline
        return func(data.stress, data.strain; extrapolation = ExtrapolationType.Linear)
    elseif func == BSplineApprox
        return func(data.stress, data.strain, 3, 4, :ArcLen, :Average; extrapolation = ExtrapolationType.Linear)
    end
    return nothing
end



# function get_initial_ss()
#     strain = LinRange(0, 0.1, 10)
#     stress = 100.0 .* strain .^ 0.2
#     return (;strain, stress)
# end


function get_initial_SSE(strain, stress)

    SSE = Dict{String, Any}(
        "name" => "Material", 
        "is true" => false,
        "rawdata" => (;strain, stress),
        "hardening offset" => 2e-3,
        "export density"=> 100,
        
        )
    
    push!(SSE, "modulus" => SS.get_modulus(SSE["rawdata"]))
    push!(SSE, "true stress" => SS.engineering_to_true(SSE["rawdata"]))
    
    push!(SSE, "hardening" => SS.get_hardening_portion(SSE["true stress"], 
                                                       SSE["modulus"]; 
                                                       offset = SSE["hardening offset"])
        )

    push!(SSE, "export max strain" => maximum(SSE["rawdata"].strain))


    return SSE
end


function initialize()
    strain = LinRange(0, 0.1, 20)
    stress = 100.0 .* strain .^ 0.2
    SSE = get_initial_SSE(strain ,stress)


    fitfuncs = [LinearInterpolation, 
                CubicSpline,
                BSplineApprox, 
                SS.Swift, 
                SS.Voce, 
                SS.HockettSherby, 
                SS.StoughtonYoon]

    fitfunclabels = ["Linear", 
                    "CSplines", 
                    "BSplineApprox",
                    "Swift", 
                    "Voce", 
                    "HockettSherby",
                    "StoughtonYoon"]

    
    push!(SSE, "interpolant" => fitfuncs[1])
    push!(SSE, "hardening fit" => get_interpolant(SSE["interpolant"], SSE["hardening"]))
    
    return (;SSE, fitfuncs, fitfunclabels)
end


function update_SSE!(SSE, data = nothing; modulus = nothing)
    
    if !isnothing(data)
        SSE["rawdata"] = (;strain = data.strain,
                           stress = data.stress) 
    end

    
    SSE["true stress"] = SSE["is true"] ? SSE["rawdata"] : SS.engineering_to_true(SSE["rawdata"])

    SSE["modulus"] = isnothing(modulus) ? SS.get_modulus(SSE["true stress"]) : modulus
    SSE["hardening"] = SS.get_hardening_portion(SSE["true stress"], 
                                                SSE["modulus"]; 
                                                offset = SSE["hardening offset"]
                                                )

    if isnothing(SSE["hardening"])
        @show SSE["hardening offset"]
        @show extrema(SSE["true stress"].strain)
        @show SSE["modulus"]
        error("Hardening Curve could not be extracted2!")
    end
    SSE["hardening fit"] = get_interpolant(SSE["interpolant"], SSE["hardening"])

    return nothing
end

function change_true_status!(val::Bool, SSE)

    SSE["is_true"] = val
    update_SSE!(SSE)
    return !val
end

function export_data(SSE)
    

end

function min_slope(SSE)
    ss = SSE["true stress"]
    return last(ss.stress) / last(ss.strain)
end

function initialize_slider(SSE)
    E = SSE["modulus"]
    Emin = round(min_slope(SSE); sigdigits =3)
    Emax = round(5 * E; sigdigits = 3)
    return (;value = E, vmin = Emin, vmax = Emax)
end