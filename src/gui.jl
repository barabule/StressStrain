

function main(; 
                N=1000,
                sidebar_width = 300,
                alg = NelderMead(),
                resample_density = 20,
                )

    fig = Figure()
    
    SSE, fitfuncs, fitfunclabels = initialize(; resample_density)

    axss = Axis(fig[1,1], title = "Stress Strain",
                    xlabel = "True Strain [-]",
                    ylabel = "True Stress [MPa]")


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


    sldvals = get_slider_range_values(SSE)
    sld_modulus = Slider(subgl1[8, :], 
                    range = LinRange(sldvals.vmin, sldvals.vmax, 1001),
                    startvalue= sldvals.value,
                    update_while_dragging =true,
                    width = sidebar_width,
                    )

    lab_modulus = Label(subgl1[7, :], "E = $(round(sld_modulus.value[]; sigdigits= 3))MPa")


    sld_toein = Slider(subgl1[10, :],
                range = LinRange(0, maximum(SSE["rawdata"].strain), 100),
                startvalue = 0.0,
                update_while_dragging = true,
                width = sidebar_width,
                )
   
    label_toein = Label(subgl1[9, :], "??")
    label_toein.text[] = "Toein = " * string(round(sld_toein.value[]; sigdigits = 3))

    subgl2[1,1] = vgrid!(
            Label(fig, "Fitting Function", width = nothing),
            fit_menu,
            # (Label("True", alignmode = :right), Checkbox(checked = false)),
            ;
            tellheight = false, width = 200,
    )    


    label_status = Label(gl_bot[1,1], "Status", tellwidth = false)

    ####################EVENTS########################################################################        
    on(cb_true.checked) do val
        SSE["is true"] = val
        update_SSE!(SSE; alg)
        plot_stress!(axss, SSE; N)
        update_status_label!(label_status, SSE)
       
    end

    
    on(sld_modulus.value) do val
        update_SSE!(SSE; modulus = round(val; sigdigits =3), alg)
        plot_stress!(axss, SSE; N)
        lab_modulus.text = "E = $(round(sld_modulus.value[]; sigdigits= 3))MPa"
        update_status_label!(label_status, SSE)
    end
    

    on(sld_toein.value) do val
        SSE["toein"] = val
        update_SSE!(SSE)
        plot_stress!(axss, SSE; N)
        update_status_label!(label_status, SSE)
    end


    on(fit_menu.selection) do s
        
        SSE["interpolant"] = s
        update_SSE!(SSE; alg)
        update_status_label!(label_status, SSE)
        plot_stress!(axss, SSE; N)
    end

    on(sld_toein.value) do val
        label_toein.text[] = "Toein = " * string(round(val; sigdigits = 3))
    end



    on(events(fig).dropped_files) do files
        isempty(files) && return nothing
        
        f1 = first(files)
        println(f1)
        data = try_open(f1)
        if !isnothing(data)
            SSE["rawdata"] = data
            update_SSE!(SSE; alg)
            
            sldvals = get_slider_range_values(SSE)
            sld_modulus.range = LinRange(sldvals.vmin, sldvals.vmax, 1001)
            _ = set_close_to!(sld_modulus, sldvals.value)
            
            plot_stress!(axss, SSE; N)
            update_status_label!(label_status, SSE)
        end
    end

    on(tb_extrapolation_strain.stored_string) do s
        SSE["export max strain"] = clamp(parse(Float64, s), last(SSE["hardening"].strain), Inf)
        plot_stress!(axss, SSE; N)
        update_status_label!(label_status, SSE)
    end

    on(tb_name.stored_string) do s
        SSE["name"] = s
        axss.title = s
    end

    return fig

end


function try_open(fn::AbstractString,
                    skipstart = 0)


    # rawdata = readdlm(fn)
    # return (;strain = rawdata[:,1],
    #         stress = rawdata[:, 2])
    delims = [',','\t',';']

    for delim in delims
        rawdata = readdlm(fn, delim;skipstart)
        size(rawdata, 2) >= 2 || continue
        eltype(rawdata) == Any && continue
        return (;strain = rawdata[:,1], 
                stress = rawdata[:, 2])
    end
    return nothing
end

function plot_stress!(ax, SSE; N = 100, tmax = SSE["export max strain"])

    empty!(ax)
    #true stress
    tss = SSE["true stress"]
    if !SSE["is true"]
        RD = SSE["rawdata"]
        scatterlines!(ax, RD.strain, RD.stress, color = (:grey10, 0.5), markersize = 10, label = "Raw Data")

    end
    scatterlines!(ax, tss.strain, tss.stress,
                    label= "True Stress (Exp)",
                        color = :black, marker = '◉', alpha = 0.3, 
                        markercolor = :black, 
                        markersize = 20)
    t = LinRange(0, tmax, N)
    hss = SSE["hardening"]
    scatter!(ax, hss.strain, hss.stress, 
                label = "Hardening Portion (Exp)",
                color = :red, marker = '◉', alpha = 0.3, markersize = 20)

    lines!(ax, t , SSE["hardening fit"](t), color = :black, label = "Hardening Law (Fit)", linewidth = 2.0)
    
    #plot the linear portion
    E = SSE["modulus"]
    smax = maximum(SSE["hardening"].stress)
    tmax = smax / E
    tlin = LinRange(0, tmax, 10)
    slin = E .* tlin
    lines!(ax, tlin, slin, linestyle = :dash, color = :red)

    vlines!(ax, [SSE["toein"]], color= :black, linestyle = :dash)

    return nothing
end


function get_initial_SSE(strain, stress; resample_density = 20)

    SSE = Dict{String, Any}(
        "name" => "Material", 
        "is true" => false, #is rawdata true stress ? 
        "rawdata" => (;strain, stress), #do not change
        "hardening offset" => 2e-3, #offset to use to extract hardening curve
        "export density"=> 100, #number of points to export
        "resample density" => resample_density, #number of points to resample rawdata
        "resample mode"=> "decimate", #how to resample ; decimate, uniform linear, uniform cspline, approx
        "toein" => 0.0,        
        )
    
    push!(SSE, "modulus" => get_modulus(SSE["rawdata"]))
    push!(SSE, "true stress" => engineering_to_true(SSE["rawdata"]))
    
    push!(SSE, "hardening" => get_hardening_portion(SSE["true stress"], 
                                                       SSE["modulus"]; 
                                                       offset = SSE["hardening offset"])
        )

    push!(SSE, "export max strain" => maximum(SSE["rawdata"].strain))


    return SSE
end


function initialize(;
                    alg = NelderMead(), 
                    resample_density = 20,
                    )

    strain = LinRange(0, 0.1, 20)
    stress = 100.0 .* strain .^ 0.2
    SSE = get_initial_SSE(strain ,stress; resample_density)


    fitfuncs = [LinearInterpolation, 
                CubicSpline,
                BSplineApprox,
                PCHIPInterpolation,
                Bilinear,
                Swift, 
                Voce, 
                HockettSherby, 
                StoughtonYoon]

    fitfunclabels = ["Linear", 
                    "CSplines", 
                    "BSplineApprox",
                    "PCHIPInterpolation",
                    "Bilinear",
                    "Swift", 
                    "Voce", 
                    "HockettSherby",
                    "StoughtonYoon"]

    
    push!(SSE, "interpolant" => fitfuncs[1])
    push!(SSE, "hardening fit" => make_interpolant(SSE["interpolant"], SSE["hardening"]; alg))
    
    return (;SSE, fitfuncs, fitfunclabels)
end


function update_SSE!(SSE; modulus = nothing, alg = NelderMead())
    
    

    
    if SSE["toein"]>0.0
        RD = toein_compensate(SSE["rawdata"]; cut = SSE["toein"])
    else
        RD = SSE["rawdata"]
    end
    SSE["true stress"] = SSE["is true"] ? RD : engineering_to_true(RD)

    SSE["modulus"] = isnothing(modulus) ? get_modulus(SSE["true stress"]) : modulus
    SSE["hardening"] = get_hardening_portion(SSE["true stress"], 
                                                SSE["modulus"]; 
                                                offset = SSE["hardening offset"]
                                                )

    if isnothing(SSE["hardening"])
        @show SSE["hardening offset"]
        @show extrema(SSE["true stress"].strain)
        @show SSE["modulus"]
        error("Hardening Curve could not be extracted2!")
    end
    SSE["hardening fit"] = make_interpolant(SSE["interpolant"], SSE["hardening"]; alg)

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
    smin = Inf
    for i in eachindex(ss.strain)
        i==firstindex(ss.strain) && continue
        slp = (ss.stress[i] - ss.stress[i-1]) / (ss.strain[i] - ss.strain[i-1])
        if smin > slp
            smin = slp
        end
    end
    return smin
end

function get_slider_range_values(SSE)
    E = SSE["modulus"]
    Emin = round(min_slope(SSE); sigdigits =3)
    Emax = round(5 * E; sigdigits = 3)
    return (;value = E, vmin = Emin, vmax = Emax)
end


function update_status_label!(label, SSE)

    label.text[] = interpolant_label(SSE["hardening fit"], SSE["interpolant"])

end
