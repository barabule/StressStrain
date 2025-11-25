

function main(data = nothing; 
                clean_data = true, #if called directly with some data
                N=1000,
                sidebar_width = 300,
                bottom_panel_height = 100,
                subscale = 0.9,
                alg = NelderMead(),
                resample_density = 20,
                )

    #the default 'data'
    if isnothing(data)
        ___strain___ = collect(LinRange(0, 0.1, 100))
        data = (;strain = ___strain___,
                stress = @. 100 * ___strain___ ^ 0.16 + rand() * 2
                )
    end

    if clean_data
            clean!(data)
    end


    #########################GLOBALS####################################################################################
    
    sidebar_sub_width = subscale * sidebar_width
    bottom_panel_sub_height = subscale * bottom_panel_height

    SSE, fitfuncs, fitfunclabels, resamplefuncs, resamplefunclabels = initialize(data; resample_density)
    
    
    #holds the state
    CURVEDATA = Dict{Symbol, Any}(
            :name => "Material",
            # 
            :sidebar_width => sidebar_width,
            :sidebar_sub_width => sidebar_sub_width,
            # resampling 
            :resample_menu_options => zip(resamplefunclabels, resamplefuncs),
            :menu_hardening_fit_options => zip(fitfunclabels, fitfuncs),
            :resampler => first(resamplefuncs), #resampling function
            # elastic properties
            :e_modulus => 1.0,
            :max_elastic_range => 5e-3,
            :is_fixed_modulus => false, #is the modulus fixed during recomputations ?
            #data source
            :is_true_stress => false,
            :rawdata => data, #this is the actual raw data from the import- must not be changed
            :base_data => deepcopy(data), #this is the basis for all the computations - can be changed
            # tunable
            :hardening_offset => 2e-3,
            :toe_in => 0.0, #cut the beginning of base_data and shift
            :cutoff => last(data.strain), # cut at the end
            #recomputable
            :true_stress =>(;strain = [], stress = []),
            :hardening =>(;strain = [], stress = []),
            # export
            :export_density => 100, # how many pts to export
            :export_true_stress => true,
            :export_hardening => true,
            :export_plot => true,
            :export_format_delim => ',',
            :export_format_ext => ".csv",
            :export_folder => nothing, #where to export
            
    )


    fig = Figure(title = "Elasto Plastic Fitter")
    
    axss = initialize_axis(fig)


    update_stress_plot!(axss, SSE; N)
    
    # xlims!(axss, low = 0)
    ylims!(axss, low = 0) #don't care about negative crap
    
    sidebar = GridLayout(fig[1,2], width = sidebar_width, tellheight = false) #holds most controls

    gl_bot = GridLayout(fig[2, 1], height = bottom_panel_height, tellwidth = false)
    gl_bot_sub = GridLayout(gl_bot[1,1], alignmode = Outside(20), height = bottom_panel_sub_height)


    ################# SIDEBAR BLOCKS #############################################################################

    for (i, entry) in enumerate(
                                    (   
                                    (draw_overview_controls!, "Overview"),
                                    (draw_true_stress_controls!, "True Stress"),
                                    (draw_emodulus_controls!, "E modulus"),
                                    (draw_hardening_controls!, "Hardening"),
                                    (draw_export_controls!, "Export"),
                                    )
                                )

        draw_controls!, btn_label = entry
        gl_block = GridLayout(sidebar[i, 1])

        controls = GridLayout()
        draw_controls!(fig, controls, CURVEDATA)
        make_button_block!(gl_block, controls;
                            btn_label,
                            btn_width = CURVEDATA[:sidebar_sub_width]/3)
                            
        Box(gl_block[:, :], linestyle = :solid,
            color = :grey90,
            width = CURVEDATA[:sidebar_width],
            tellwidth = true,
            tellheight = false,
            alignmode = Outside(-5),
            cornerradius = 10,
            z = -100)

    end

    ################## Bottom Panel ####################################################################################

    sld_int = IntervalSlider(gl_bot_sub[1,2], range = LinRange(0, last(SSE["true stress"].strain), 1000),
                    startvalues = (0.0, last(SSE["true stress"].strain)))
    
    
    
    labeltext_sld_int = lift(sld_int.interval) do int
        t1 = "Toe in\n" * string(round(int[1], digits=3))
        t2 = "Cut off\n" * string(round(int[2], digits=3))
        (t1, t2)
    end
    
    Label(gl_bot_sub[1,1],  @lift $labeltext_sld_int[1]) 
    Label(gl_bot_sub[1,3],  @lift $labeltext_sld_int[2])

    label_status = Label(gl_bot_sub[2,:], "Status", tellwidth = false)

    Box(gl_bot[1,1], 
            linestyle = :solid,
            height = bottom_panel_height,
            color = :grey90,
            tellwidth = false,
            z = -100)

    ################# EVENTS ###########################################################################################
    
    on(events(fig).dropped_files) do files
        isempty(files) && return nothing
        
        f1 = first(files)
        println(f1)
        
        data_gui(screen, f1)
        update_stress_plot!(axss, SSE)
    end


    on(sld_int.interval) do intv
        lo, hi = intv
        SSE["toein"] = lo
        SSE["cut off"] = hi
        update_SSE!(SSE; recompute_modulus = true)
        
        update_modulus_slider!(sld_modulus, SSE)
        update_stress_plot!(axss, SSE;N)
        update_status_label!(label_status, SSE)
    end


    ##################    WINDOW   ################################################################
    screen = GLMakie.Screen()
    GLFW.SetWindowTitle(screen.glscreen, "Stress Strain Fitter")
    return display(screen, fig)

end




function update_stress_plot!(ax, SSE; 
                        N = 100, 
                        tmax = SSE["export max strain"],
                        )

    empty!(ax)
    fig = ax.parent
    #delete the legend
    for (i, block) in enumerate(fig.content)
        if block isa Makie.Legend
            # 3. If found, delete it using the delete! function
            Makie.delete!(block)
        end
    end


    ax.title = SSE["name"]
    #original data
    if haskey(SSE, "original data")
        OD = SSE["original data"]
        scatter!(ax, OD.strain, OD.stress, color = (:grey90, 0.7), markersize = 10, label = "original")
    end
    #rawdata
    tss = SSE["true stress"]
    if !SSE["is true"]
        ### Engineering Stress Plot
        RD = SSE["rawdata"]
        scatterlines!(ax, RD.strain, RD.stress, 
                    color = (:grey10, 0.5), 
                    marker = 'o', markersize = 10, label = "Raw Data (Engineering)")

    end
    ### True Stress Plot
    scatterlines!(ax, tss.strain, tss.stress,
                    label= "True Stress (Exp)",
                        color = (:grey50, 0.5), marker = :diamond, 
                        markercolor = (:black, 0.5), 
                        markersize = 20)

    ### Hardening Curve Plot
    t = LinRange(0.0, tmax, N)
    hss = SSE["hardening"]
    scatter!(ax, hss.strain, hss.stress, 
                label = "Hardening Portion (Exp)",
                color = :red, marker = :cross, alpha = 0.9, markersize = 20)

    lines!(ax, t , SSE["hardening fit"](t), color = :black, label = "Hardening Law (Fit)", linewidth = 4.0)
    
    #E modulus plot
    E = SSE["modulus"]
    smax = maximum(SSE["hardening"].stress)
    tmax = smax / E
    tmax = min(last(SSE["true stress"].strain), tmax)
    tlin = LinRange(0, tmax, 10)
    slin = E .* tlin
    lines!(ax, tlin, slin, linestyle = :dash, color = :red)
    #hardening offset
    lines!(ax, tlin .+ SSE["hardening offset"],slin, color = (:red, 0.3), linestyle = :solid )
    #elastic range
    vlines!(ax, [SSE["elastic range"]], color = (:grey50, 0.5), linestyle = :solid)
    ### Cutoff limits
    vlines!(ax, [SSE["toein"], SSE["cut off"]], color= :black, linestyle = :dash)

    axislegend(ax, position = :rb, merge = true)

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
        "resampler"=> LinearInterpolation,
        "toein" => 0.0,    
        "cut off" => last(strain), #cut off value for true stress curve    
        "export folder" => nothing,
        "fixed modulus"=> false, #if this is set, modulus cannot change
        "elastic range" => 1e-3, #how much of the starting portion to use to compute modulus
        "max elastic range" => last(strain), #max strain value for elastic behavior 
        "BSpline approximation knots" => 4, #how many control points for BSpline approx
        "BSpline degree" => 3, 
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


function initialize(data = nothing;
                    alg = NelderMead(), 
                    resample_density = 20,
                    )

    if isnothing(data) #make up something
        strain = LinRange(0, 0.1, 20)
        stress = 100.0 .* strain .^ 0.2
    else
        strain, stress = data
    end

    SSE = get_initial_SSE(strain ,stress; resample_density)

    if haskey(data, :folder)
        SSE["export folder"] = data.folder
    end

    fitfuncs = [LinearInterpolation, 
                AkimaInterpolation,
                BSplineApprox,
                PCHIPInterpolation,
                Bilinear,
                Swift, 
                Voce, 
                HockettSherby, 
                StoughtonYoon,
                SwiftVoce,
                ]

    fitfunclabels = ["Linear", 
                    "Akima", 
                    "BSplineApprox",
                    "PCHIPInterpolation",
                    "Bilinear",
                    "Swift", 
                    "Voce", 
                    "HockettSherby",
                    "StoughtonYoon",
                    "Swift+Voce",
                    ]

    resamplefuncs = [
            LinearInterpolation,
            Hollomon,
            CubicSpline,
            BSplineApprox,
            QuadraticSpline,
            AkimaInterpolation,
            moving_average,
            RegularizationSmooth,
            RamberOsgoodAlternativeReparametrized,
    ]

    resamplefunclabels = [
        "Linear",
        "Hollomon",
        "CSplines",
        "BSpline",
        "Quadratic",
        "Akima",
        "Moving Average",
        "Reg Smooth",
        "RambergOsgoodAlt",
    ]

    push!(SSE, "interpolant" => fitfuncs[1])

    if !issorted(SSE["hardening"].strain)
        @warn "strain is not sorted", SSE["hardening"].strain[1:100]
        @info "E modulus", SSE["modulus"]
    end

    push!(SSE, "hardening fit" => make_interpolant(SSE["interpolant"], SSE["hardening"]; alg))
    
    return (;SSE, fitfuncs, fitfunclabels, resamplefuncs, resamplefunclabels)
end


function update_SSE!(SSE; 
                        modulus = nothing,  #impose modulus
                        recompute_modulus = true, #also triggers recalc
                        alg = NelderMead(),
                        resample = false, #modifies rawdata.. 
    )
    
    if resample #modifies the original data
        if !haskey(SSE, "original data") #first time
            push!(SSE, "original data" => SSE["rawdata"])
        end
        RD = SSE["rawdata"]
        if SSE["resampler"] != RambergOsgood || SSE["resampler"] 
            resampled_data = resample_curve(RD.strain, RD.stress, SSE["resample density"];
                                    resampler = SSE["resampler"],
                                    d =3,
                                    h = SSE["BSpline approximation knots"],
                                    )

            SSE["rawdata"] = (;strain = resampled_data[1], 
                            stress = resampled_data[2])
        end
    end
    
    if SSE["toein"]>0.0
        RD = toein_compensate(SSE["rawdata"]; 
                                cut = SSE["toein"], 
                                elastic_strain_offset = SSE["elastic range"],
                            )
    else
        RD = SSE["rawdata"]
    end

    TS = SSE["is true"] ? RD : engineering_to_true(RD)

    SSE["true stress"] = cutoff(TS, SSE["cut off"])

    if !SSE["fixed modulus"] #only allowed to change if modulus is not fixed
        if !isnothing(modulus)
            SSE["modulus"] = modulus
        elseif recompute_modulus 
            SSE["modulus"] = get_modulus(SSE["true stress"]; 
                                            max_strain = SSE["elastic range"],
                                        )
        end
    end
    
    if SSE["resampler"] == RambergOsgood
        resampled_data = resample_curve(SSE["true stress"].strain, 
                                        SSE["true stress"].stress, 
                                        length(SSE["true stress"].strain); #it's actually a fit
                                        resampler = SSE["resampler"],
                                        offset = SSE["hardening offset"],
                                        E = SSE["modulus"],
                                        )
        #update the true stress
        SSE["true stress"] = (;strain = resampled_data[1],
                                stress = resampled_data[2])
                                        
    end


    SSE["hardening"] = get_hardening_portion(SSE["true stress"], 
                                                SSE["modulus"]; 
                                                offset = SSE["hardening offset"]
                                                )
    #update allowable elastic range
    SSE["max elastic range"] = last(SSE["true stress"].strain)

    if isnothing(SSE["hardening"]) #how to handle this gracefully
        @warn "Hardening portion cannot be extracted!"
        @info "Modulus ", SSE["modulus"]
        @info "Hardening offset ", SSE["hardening offset"]
        
    end
    
    SSE["hardening fit"] = make_interpolant(SSE["interpolant"], SSE["hardening"]; alg)

    return nothing
end





function change_true_status!(val::Bool, SSE)

    SSE["is_true"] = val
    update_SSE!(SSE)
    return !val
end





function get_slider_range_values(SSE; sigdigits = 2)
    E = SSE["modulus"]
    Ebracket = bracket_modulus(SSE["true stress"])
    Emin = Ebracket.Emin
    Emax = Ebracket.Emax
    if E < Emin || E > Emax
        E = Emax
        update_SSE!(SSE; modulus = E)
    end
    Emin = round(Emin;sigdigits)
    Emax = 10round(Emax;sigdigits)
    return (;value = E, vmin = Emin, vmax = Emax)
end


function update_modulus_slider!(sld, SSE)

    sldvals = get_slider_range_values(SSE)
    sld.range = LinRange(sldvals.vmin, sldvals.vmax, 1001)
    _ = set_close_to!(sld, sldvals.value)
    return nothing
end



function update_status_label!(label, SSE)
    text = interpolant_label(SSE["hardening fit"], SSE["interpolant"])
    # println(text)
    E = "E = $(SSE["modulus"]) MPa, "
    label.text = E * text
    # println("label changed")
end



function reset_SSE!(SSE)
    if haskey(SSE, "original data") 
        
        delete!(SSE, "rawdata")
        push!(SSE, "rawdata" => SSE["original data"])
        delete!(SSE,"original data")
        update_SSE!(SSE)

    end
    
    return nothing

end





function initialize_axis(fig)
    return Axis(fig[1,1], title = "Stress Strain",
                    xlabel = "True Strain [-]",
                    ylabel = "True Stress [MPa]")

end


#TODO may be nicer plot
function export_data(SSE; 
                    delim = ',',
                    export_true = true,
                    export_hardening = true,
                    export_plot = true,
                    plot_export_format = :png,
                    px_per_unit = 2, #output size scaling
                    )

    true_stress = SSE["true stress"]
    
    hardening_func = SSE["hardening fit"]
    tmax = SSE["export max strain"]
    tout = collect(LinRange(0, tmax, SSE["export density"]))
    sout = hardening_func(tout)

    output_folder =  isnothing(SSE["export folder"]) ? pwd() : SSE["export folder"]

    fname = SSE["name"]
    if export_hardening
        fname1 = joinpath(output_folder, fname * "_hardening.csv")
        writedlm(fname1,
                    [tout sout], 
                    delim,
                    )
    end

    if export_true
        fname2 = joinpath(output_folder, fname * "_true stress.csv")             
        writedlm(fname2,
                    [true_stress.strain true_stress.stress],
                    delim,
                    )
    end

    if export_plot
        # fext = plot_export_format==:png ? ".png" : ".svg" #doesn't work in GLMakie...
        fext = ".png"
        fname3 = joinpath(output_folder, fname * "_plot" * fext)
        
        fig = Figure(size = (1200,800))
        ax = initialize_axis(fig)
        update_stress_plot!(ax, SSE; N = SSE["export density"])

        label_status = Label(fig[2,1],tellwidth = false)
        update_status_label!(label_status, SSE)
        
        screen = GLMakie.Screen()
        display(screen, fig)
        save(fname3, fig; px_per_unit)
        close(screen)
    end


    @info  "Done"
end

function draw_overview_controls!(fig::Figure, Lay::GridLayout, D::Dict{Symbol, Any})
    @assert hasallkeys(D, [:is_true_stress, :name, :sidebar_sub_width])

    cb_true_stress = Checkbox(fig, checked = false)
    
    tb_name = Textbox(fig, 
             width = D[:sidebar_sub_width],
             placeholder = "Material Name",
             boxcolor = :white)

    Lay[1, 1] = vgrid!(
                            hgrid!(cb_true_stress, Label(fig, "True"), halign = :left),
                            # Label(fig, "Material Name", halign = :left),
                            tb_name,
                        )

    ############### BEHAVIOR ###########################################################################################

    on(cb_true_stress.checked) do val
        D[:is_true_stress] = val
        #recalc
        #update plot
        #update status
    end

    on(tb_name.stored_string) do s
        D[:name] = s
        #update plot title
    end

    ####################################################################################################################

                        
    return nothing

end


function draw_true_stress_controls!(fig::Figure, Lay::GridLayout, D::Dict{Symbol, Any})

    @assert hasallkeys(D, [:resample_menu_options, :sidebar_sub_width])

    resample_menu = Menu(fig, options = D[:resample_menu_options],
                                default = "Linear", width = D[:sidebar_sub_width] * 0.5)

    tb_resample = Textbox(fig, placeholder = "Enter number",
                    validator = Int, tellwidth = false,
                    boxcolor = :white)

    btn_manual = Button(fig, label = "Manual")
    bezier_result_ref = Ref(Dict{String, Any}("status" => -1)) # holds the resampling result for bezier

    
    btn_reset = Button(fig, label = "Reset!")

    btn_bspline_plus = Button(fig, label="+")
    btn_bspline_minus = Button(fig, label = "-")
    lab_bspline_control_pts = Label(fig, "4")

    
    Lay[1,1] = vgrid!(
                    Label(fig, "Resample Function", width = nothing),
                    hgrid!(resample_menu, btn_bspline_plus, lab_bspline_control_pts, btn_bspline_minus),
                    hgrid!(Label(fig, "Resample"), tb_resample),
                    hgrid!(btn_manual, btn_reset),
                    ;
                    tellheight = false, 
                    width = D[:sidebar_sub_width],

    )

    ###### BEHAVIOR #######
    on(resample_menu.selection) do s
        # #TODO better handling of modulus update when fitting on true stress
        # SSE["resampler"] = s
        # # @info "SSE", SSE["resampler"]
        # update_SSE!(SSE;resample = true, recompute_modulus = false)
        # update_stress_plot!(axss, SSE; N)
        # update_status_label!(label_status, SSE)
    end
    on(btn_bspline_plus.clicks) do _
        # nc = SSE["BSpline approximation knots"]
        # nc += 1
        # SSE["BSpline approximation knots"] = nc
        # lab_bspline_control_pts.text = string(nc)
        # update_SSE!(SSE)
        # update_stress_plot!(axss, SSE)
        # update_status_label!(label_status, SSE)
    end

    on(btn_bspline_minus.clicks) do _ 
        # nc = SSE["BSpline approximation knots"]
        # nc = max(4, nc - 1)
        # SSE["BSpline approximation knots"] = nc
        # lab_bspline_control_pts.text = string(nc)
        # update_SSE!(SSE)
        # update_stress_plot!(axss, SSE)
        # update_status_label!(label_status, SSE)
    end

    on(tb_resample.stored_string) do s
        # SSE["resample density"] = clamp(parse(Int, s), 2, 10_000)
        # update_SSE!(SSE;alg, resample = true, recompute_modulus = false)
        # update_stress_plot!(axss, SSE; N)
        # update_status_label!(label_status, SSE)
    end

    on(btn_manual.clicks) do _
        # #open a new window, do the fitting,
        # #then close the main window and restart with new data
        # # handle_bezier_fit(SSE, screen)
        # @async begin
        #     @info "Opening Bezier Fitting Window"
        #     bezier_result_ref[]["status"] = 0 #reset the status
        #     bezier_result_ref[]["data"] = SSE["rawdata"]
        #     CubicPiecewiseBezier.bezier_fit_fig(bezier_result_ref)
            
        #     @info "Closed Bezier Fitting Window"
            
            
        #     results = bezier_result_ref[] #get the last state
        #     # @info "Status", results["status"]
        #     bezier_fit = results["bezier fit"]
        #     #no longer necessary, this returns directly what we need
        #     # strain = [pt[1] for pt in bezier_fit]
        #     # stress = [pt[2] for pt in bezier_fit]
        #     if !haskey(SSE, "original data")#backup
        #         push!(SSE, "original data" => SSE["rawdata"])
        #     end
        #     SSE["rawdata"] = bezier_fit
        #     update_SSE!(SSE)
        #     update_stress_plot!(axss, SSE; N)
        #     update_status_label!(label_status, SSE)

        # end
         
    end

    on(btn_reset.clicks) do _
        # reset_SSE!(SSE)
        # update_stress_plot!(axss, SSE; N)
    end

    return nothing
end



function draw_emodulus_controls!(fig::Figure, Lay::GridLayout, D::Dict{Symbol, Any})

    @assert hasallkeys(D, [:max_elastic_range, :sidebar_sub_width])

    # sldvals = get_slider_range_values(SSE)
    sldvals = (;vmin = 10.0, vmax = 1e6, value = 50_000.0)


    sld_modulus = Slider(fig, 
                    range = LinRange(sldvals.vmin, sldvals.vmax, 1001),
                    startvalue= sldvals.value,
                    update_while_dragging =true,
                    width = D[:sidebar_sub_width],
                    )

    max_elastic_range = D[:max_elastic_range]

    sld_elastic_range = Slider(fig,
                    range = LinRange(0.0, max_elastic_range, 1000),
                    startvalue = 1e-3,
                    update_while_dragging = true,
                    width = 0.4 * D[:sidebar_sub_width],
                    )

    lab_elastic_range_val = Label(fig, "0.001")

    lab_modulus = Label(fig, "E = ")

    tb_modulus = Textbox(fig, 
                            placeholder = "$(round(sld_modulus.value[]; sigdigits= 3))MPa",
                            validator = Float64,
                            width = 0.4 * D[:sidebar_sub_width],
                            )


    sld_offset = Slider(fig,
                range= LinRange(1e-3, 0.1, 100),
                startvalue = 2e-3,
                update_while_dragging = false,
                
                )

    cb_fixed = Checkbox(fig, checked = false)
    
    lab_offset = Label(fig, "Offset = $(round(sld_offset.value[]; sigdigits = 3))")
    

    Lay[1, 1] = vgrid!(
        hgrid!(lab_modulus, tb_modulus, Label(fig, "MPa")),
        sld_modulus,
        hgrid!(Label(fig, "Fix Modulus"), cb_fixed),
        hgrid!(Label(fig, "Elastic range "), sld_elastic_range, lab_elastic_range_val),
        hgrid!(lab_offset, sld_offset),
        ;
        width = D[:sidebar_sub_width],
    )

#########BEHAVIOR ####################################################################
    on(sld_modulus.value) do val
        # update_SSE!(SSE; modulus = round(val; sigdigits =3), alg)
        # update_stress_plot!(axss, SSE; N)
        # tb_modulus.displayed_string = "$(round(sld_modulus.value[]; sigdigits= 3))"
        # update_status_label!(label_status, SSE)
    end

    on(sld_elastic_range.value) do val
        # SSE["elastic range"] = val
        # lab_elastic_range_val.text = string(round(val;sigdigits = 4))
        # update_SSE!(SSE)
        # update_stress_plot!(axss, SSE; N)
        # update_status_label!(label_status, SSE)
    end

    on(tb_modulus.stored_string) do s
        #update modulus only in the checkbox is unchecked
        # if !cb_fixed.checked[]
        #     modulus = round(parse(Float64, s); sigdigits = 3)
        #     update_SSE!(SSE; modulus)
        #     update_stress_plot!(axss, SSE;N)
        #     update_status_label!(label_status, SSE)
        # end
    end

    
    on(cb_fixed.checked) do val
        # SSE["fixed modulus"] = val
    end

    on(sld_offset.value) do val
        # SSE["hardening offset"] = val
        # update_SSE!(SSE;recompute_modulus =false)
        # update_stress_plot!(axss, SSE; N)
        # lab_offset.text = "Offset = $(round(sld_offset.value[]; sigdigits = 3))"
    end

end


function draw_hardening_controls!(fig::Figure, Lay::GridLayout, D::Dict{Symbol, Any})

    @assert hasallkeys(D, [:sidebar_sub_width, :menu_hardening_fit_options])

    fit_menu = Menu(fig, options = D[:menu_hardening_fit_options],
                    default = "Linear",
                    width = 0.6 * D[:sidebar_sub_width])

    tb_hardening_pts = Textbox(fig, placeholder = "Enter number",
                        validator = Int,
                        width = 0.5 * D[:sidebar_sub_width],
                        )

    tb_extrapolation_strain = Textbox(fig, 
                        validator = Float64,
                        width = 50,
                        halign=:left,
                        boxcolor = :white,
                        )

    

    Lay[1, 1] = vgrid!(
            # Label(fig, "Hardening Curve", fontsize = 20, font =:italic),
            hgrid!(Label(fig, "Method"), fit_menu),
            hgrid!(Label(fig, "Num Pts"), tb_hardening_pts),
            hgrid!(Label(fig, "Extrapolate strain to"), tb_extrapolation_strain),
            # (Label("True", alignmode = :right), Checkbox(checked = false)),
            ;
            halign = :left,
            tellheight = false, 
            width = D[:sidebar_sub_width],
    )    
    ########### BEHAVIOR #############
    on(fit_menu.selection) do s
        
        # SSE["interpolant"] = s
        # update_SSE!(SSE; alg, recompute_modulus = false) #don't recalc modulus
        # update_status_label!(label_status, SSE)
        # update_stress_plot!(axss, SSE; N)
    end

    on(tb_extrapolation_strain.stored_string) do s
        
        # SSE["export max strain"] = clamp(parse(Float64, s), last(SSE["hardening"].strain), Inf)
        # update_stress_plot!(axss, SSE; N)
        # update_status_label!(label_status, SSE)
        # update_status_label!(label_status, SSE)
    end

    on(tb_hardening_pts.stored_string) do s
        # num_hardening_pts = Int(clamp(parse(Int, s), 2, Inf))
        # SSE["export density"] = num_hardening_pts
        # update_SSE!(SSE;alg)
        # update_stress_plot!(axss, SSE; N = num_hardening_pts)
    end


end


function draw_export_controls!(fig::Figure, Lay::GridLayout, D::Dict{Symbol, Any})

    @assert hasallkeys(D, [:sidebar_sub_width])

    btn_export = Button(fig, label = "Export")

    cb_exp_true = Checkbox(fig, checked = true)
    cb_exp_hardening = Checkbox(fig, checked = true)
    cb_exp_plot = Checkbox(fig, checked = true)

    Lay[1,1] = vgrid!(
                    hgrid!(cb_exp_true, Label(fig, "True Stress"), halign = :left),
                    hgrid!(cb_exp_hardening, Label(fig, "Hardening"),halign = :left),
                    hgrid!(cb_exp_plot, Label(fig, "Plot"), halign = :left),
                    btn_export,
                    )

    ########## BEHAVIOR #################

    on(btn_export.clicks) do _
        # export_data(SSE; export_true = cb_exp_true.checked[],
        #                 export_hardening = cb_exp_hardening.checked[],
        #                 export_plot = cb_exp_plot.checked[],
        #                 # plot_export_format = menu_export_format.selection,
        #                 )

    end



end


"""
    make_button_block!(main_GL::GridLayout, controls::GridLayout; #holds the widget and behavior
                        btn_label= "Show / Hide", #what to show on the button
                        btn_width = 50, 
                        is_visible = false,
                        btn_kwargs...
                        )

Generates a block with toggleable visibility via a top button.
"""
function make_button_block!(main_GL::GridLayout, controls::GridLayout; #holds the widget and behavior
                        btn_label= "Show / Hide", #what to show on the button
                        btn_width = 50,
                        is_visible = false,
                        btn_halign = :left,
                        btn_fontsize = 16,
                        btn_font = :italic,
                        ) #initial visibility

    
    btn_visible = Button(main_GL[1,1], label = btn_label, 
                        width = btn_width, 
                        halign = btn_halign,
                        fontsize = btn_fontsize,
                        font = btn_font)

    visibility = Observable(is_visible)

    sub_GL = main_GL[2, :] = GridLayout()

    hidden_GL = GridLayout(bbox = (-100, -100, 0, 0)) #offscreen
    empty_GL = GridLayout()
    sub_GL[1,1] = controls
    hidden_GL[1,1] = empty_GL
    
    on(btn_visible.clicks) do  _
        #toggle visibility
        visibility[] = !visibility[]
    end

    on(visibility) do visible
        if visible
            sub_GL[1,1] = controls
            hidden_GL[1,1] = empty_GL
            rowsize!(main_GL, 2, Auto(false))
        else 
            sub_GL[1,1] = empty_GL
            hidden_GL[1,1] = controls
            rowsize!(main_GL, 2, Fixed(0))
        end
    end
    visibility[] = is_visible #trigger
    return nothing
end