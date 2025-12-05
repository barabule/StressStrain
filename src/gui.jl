

function launch_gui(data = nothing; 
                clean_data = true, #if called directly with some data
                N=1000,
                sidebar_width = 300,
                bottom_panel_height = 100,
                subscale = 0.9,
                alg = NelderMead(),
                precompile_run = false, #only used when precompile run to close window
                screen = nothing::Union{Nothing, GLMakie.Screen} , #
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

    export_folder, fitfuncs, fitfunclabels, resamplefuncs, resamplefunclabels = initialize(data)
    
    fig = Figure(title = "Elasto Plastic Fitter")
    
    axss = initialize_axis(fig)


    #holds the state
    CURVEDATA = Dict{Symbol, Any}(
            #general
            :name => "Material",
            # gui
            :sidebar_width => sidebar_width,
            :sidebar_sub_width => sidebar_sub_width,
            :status_label => nothing, #here goes the status label handle
            :modulus_slider => nothing, #handle to the modulus_slider
            # resampling 
            :resample_menu_options => zip(resamplefunclabels, resamplefuncs),
            :menu_hardening_fit_options => zip(fitfunclabels, fitfuncs),
            :resampler => first(resamplefuncs), #resampling function
            :resample_density => 100, #how many points to resample the true stress curve
            :hardening_interpolant => first(fitfuncs), #function to fit hardening curve
            :alg => alg, #optimization algorithm for fitting hardening portion
            :max_plastic_strain => 0.0, #plastic strain to extrapolate to
            # elastic properties
            :e_modulus => 1.0,
            :max_elastic_range => 5e-3,
            :is_fixed_modulus => false, #is the modulus fixed during recomputations ?
            #data source
            :is_true_stress => false,
            :rawdata => data, #this is the actual raw data from the import- must not be changed
            # tunable
            :hardening_offset => 2e-3,
            :toein => 0.0, #cut the beginning of base_data and shift
            :cutoff => last(data.strain), # cut at the end
            #recomputable
            :base_data => deepcopy(data), #this is the basis for all the computations - can be changed by resampling
            :resample_base_data => false, #resample base_data with resampler function
            # these control the computation
            :recompute_modulus => true, 
            :recompute_true_stress => true, 
            :recompute_hardening_portion => true,
            :recompute_hardening_fit => true,
            #to be computed
            :true_stress =>(;strain = [], stress = []), #modified from base_data
            :hardening_portion =>(;strain = [], stress = []),#modified from true_stress
            :hardening_fitter => nothing, #fitted hardening function, points are not stored
            # export settings
            :export_density => 100, # how many pts to export
            :export_true_stress => true,
            :export_hardening => true,
            :export_plot => true,
            :export_format_delim => ',',
            :export_format_ext => ".csv",
            :export_folder => export_folder, #where to export
            :export_density => 100, #how many points to export
            :export_px_per_unit => 2, #plot output resolution scaling
            # plot
            :figure => fig,
            :axis => axss,  #can be changed to draw somewhere else...  
            :plot_rawdata => true, #plot the rawdata ?
            :plot_modulus => true, #plot E modulus line ?
            :plot_elastic_range =>true, #plot the elastic range ?
            :plot_base_data => true, #plot the base data ?
            :plot_true_stress => true, #plot the true stress ?
            :plot_hardening_portion => true, #plot the hardening portion ?
            :plot_hardening_fit => true, #plot the hardening fit ?
            :plot_hardening_offset => true,  #plot the hardening offset line ?
            :plot_cutoff_limits => true,  #plot the toein and cutoff vlines ?
            :plot_density => N, #how many points to plot 
            :plot_format => :png,  
            # approx bspline params
            :bspline_degree => 3,
            :bspline_num_pts => 4,

    )

    
    recompute_data!(CURVEDATA)
    update_stress_plot(CURVEDATA)

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
        draw_controls!(controls, CURVEDATA)
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

    sld_int = IntervalSlider(gl_bot_sub[1,2], range = LinRange(0, last(CURVEDATA[:true_stress].strain), 1000),
                    startvalues = (0.0, last(CURVEDATA[:true_stress].strain)))
    
    
    
    labeltext_sld_int = lift(sld_int.interval) do int
        t1 = "Toe in\n" * string(round(int[1], digits=3))
        t2 = "Cut off\n" * string(round(int[2], digits=3))
        (t1, t2)
    end
    
    Label(gl_bot_sub[1,1],  @lift $labeltext_sld_int[1]) 
    Label(gl_bot_sub[1,3],  @lift $labeltext_sld_int[2])

    label_status = Label(gl_bot_sub[2,:], "Status", tellwidth = false)
    CURVEDATA[:status_label] = label_status
    
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
        update_stress_plot(CURVEDATA)
    end


    on(sld_int.interval) do intv
        lo, hi = intv
        CURVEDATA[:toein] = lo
        CURVEDATA[:cutoff] = hi
        # update_SSE!(SSE; recompute_modulus = true)
        recompute_data!(CURVEDATA)
        
        update_modulus_slider!(CURVEDATA)
        update_stress_plot(CURVEDATA)
        update_status_label!(CURVEDATA)
    end


    ##################    WINDOW   ################################################################

    update_status_label!(CURVEDATA)

    fig[0, :] = Label(fig, "Stress-Strain Fitter. Drop a file to import data.")

    screen = isnothing(screen) ? GLMakie.Screen() : screen
    # @info typeof(screen)
    GLFW.SetWindowTitle(screen.glscreen, "Stress Strain Fitter")
    if precompile_run
        return display(screen, fig)
    end
    return wait(display(screen, fig)) #to not directly close the window after opening by script

end




function update_stress_plot(D::Dict{Symbol, Any})

    ax = D[:axis]
    fig = D[:figure]
    empty!(ax)
    
    #delete the legend
    for (i, block) in enumerate(fig.content)
        if block isa Makie.Legend
            # 3. If found, delete it using the delete! function
            Makie.delete!(block)
        end
    end

    ax.title = D[:name]

    if D[:plot_rawdata] #show very faintly
        RD = D[:rawdata]
        scatterlines!(ax, RD.strain, RD.stress, 
                color = (:grey60, 0.7), markercolor = (:grey60, 0.7),
                markersize = 10,
                linewidth = 1,
                label = "Original Data",
                )
    end

    if D[:plot_true_stress]
        BD = D[:base_data]
        if !D[:is_true_stress] #engineering 
            scatterlines!(ax, BD.strain, BD.stress, 
                            color = (:grey10, 0.5), 
                            marker = 'o', 
                            markersize = 10,
                            markercolor = :grey10, 
                            label = "Raw Data (Engineering)")
        end
        tss = D[:true_stress]
        scatterlines!(ax, tss.strain, tss.stress,
                        label= "True Stress (Exp)",
                        color = (:grey50, 0.5), marker = :diamond, 
                        markercolor = (:black, 0.5), 
                        markersize = 10)
    end

    if D[:plot_hardening_portion]
        HD = D[:hardening_portion]
        scatter!(ax, HD.strain, HD.stress,
                    color = :red,
                    markersize = 10,
                )
    end

    if D[:plot_hardening_fit] #on demand
        hss = D[:hardening_fitter]
        epi = LinRange(0, D[:max_plastic_strain], D[:plot_density])
        ssi = hss.(epi)
        lines!(ax, epi, ssi, 
                    color =:black, 
                    label = "Hardening Law (Fit)", 
                    linewidth = 4.0)
    end

    if D[:plot_modulus] #E modulus plot
        E = D[:e_modulus]
        
        smax = maximum(D[:hardening_portion].stress)
        tmax = smax / E
        tmax = min(last(D[:true_stress].strain), tmax)
        tlin = LinRange(0, tmax, 10)
        slin = E .* tlin
        lines!(ax, tlin, slin, linestyle = :dash, color = :red)

        if D[:plot_hardening_offset]
            lines!(ax, tlin .+ D[:hardening_offset], slin, color = (:red, 0.3), linestyle = :solid )
        end
    end
    
    if D[:plot_cutoff_limits]
        vlines!(ax, [D[:toein], D[:cutoff]], color= :black, linestyle = :dash)

    end
    
    if D[:plot_elastic_range]
        vlines!(ax, [D[:max_elastic_range]], color = (:grey50, 0.5), linestyle = :solid)
    end
    
    axislegend(ax, position = :rb, merge = true)

    return nothing
end




function initialize(data)

    @assert hasallkeys(data, [:strain, :stress])

    if haskey(data, :folder)
        export_folder = data.folder
    else
        export_folder = nothing
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

    return (;export_folder, fitfuncs, fitfunclabels, resamplefuncs, resamplefunclabels)
end




function recompute_data!(D::Dict{Symbol, Any})

    
    
    # 1. resample the base_data
    if D[:resample_base_data]
        BD = D[:base_data]
        resampler = D[:resampler]
        if resampler != RambergOsgood
            resampled_data = resample_curve(BD.strain, BD.stress, D[:resample_density];
                                            resampler,
                                            d = D[:bspline_degree],
                                            h = D[:bspline_num_pts],
                                            )
            
            D[:base_data] = (;strain = resampled_data[1],
                            stress = resampled_data[2])
        end

    end

    #2. check if toein
    if D[:toein] > 0
        BD = toein_compensate(D[:base_data];
                            cut = D[:toein],
                            elastic_strain_offset = D[:max_elastic_range])
        # @info "toein comp"
    else
        BD = D[:base_data]
    end

    #3. get true stress 
    TS = D[:is_true_stress] ? BD : engineering_to_true(BD)

    #4. cutoff
    D[:true_stress] = cutoff(TS, D[:cutoff])

    #5. compute E modulus
    if D[:recompute_modulus] && !D[:is_fixed_modulus]
        
        D[:e_modulus] = get_modulus(D[:true_stress];
                                        max_strain = D[:max_elastic_range])
        
    end
    
    #6. RambergOsgood 

    if D[:resample_base_data] && D[:resampler] == RambergOsgood
        resampled_data = resample_curve(D[:true_stress].strain,
                                        D[:true_stress].stress,
                                        length(D[:true_stress].strain); #optimization
                                        resampler = D[:resampler],
                                        offset = D[:hardening_offset],
                                        E = D[:e_modulus],
                                        )
        D[:true_stress] = (;strain = resampled_data[1],
                            stress = resampled_data[2])
    end
    
    
    #7. Hardening Portion
    if D[:recompute_hardening_portion]
        D[:hardening_portion] = get_hardening_portion(D[:true_stress],
                                                    D[:e_modulus];
                                                    offset = D[:hardening_offset],
                                                    )
        

    end
    # for plot /export 
    D[:max_plastic_strain] = max(D[:max_plastic_strain], last(D[:hardening_portion].strain))

    #8. update allowable elastic range
    D[:max_elastic_range] = last(D[:true_stress].strain)
    
    #9. hardening fit
    if D[:recompute_hardening_fit]
        D[:hardening_fitter] = make_interpolant(D[:hardening_interpolant], D[:hardening_portion]; alg = D[:alg]) 
    end

    return nothing

end



function get_slider_range_values(E, TT; sigdigits = 2)
    
    Ebracket = bracket_modulus(TT)
    Emin = Ebracket.Emin
    Emax = Ebracket.Emax

    Emin = round(min(Emin, E); sigdigits)
    Emax = round(max(Emax, E); sigdigits)

    return (;value = E, vmin = Emin, vmax = Emax)
end

function update_modulus_slider!(D::Dict{Symbol, Any})
    sld = D[:modulus_slider]
    E = D[:e_modulus]
    TT = D[:true_stress]
    sldvals = get_slider_range_values(E, TT)
    sld.range = LinRange(sldvals.vmin, sldvals.vmax, 1001)
    _ = set_close_to!(sld, sldvals.value)
    return nothing
end


function update_status_label!(D::Dict{Symbol, Any}, 
                label = nothing::Union{Nothing, String}, #can be overriden temporary
                            )

    label = isnothing(label) ? D[:status_label] : label
    text = interpolant_label(D[:hardening_fitter], D[:hardening_interpolant])
    # println(text)
    E = "E = $(D[:e_modulus]) MPa, "
    label.text = E * text
end



#TODO may be nicer plot (Cairomakie?)
function export_data(D::Dict{Symbol, Any})

    true_stress = D[:true_stress]
    delim = D[:export_format_delim]
    px_per_unit = D[:export_px_per_unit]

    hardening_func = D[:hardening_fitter]
    tmax = D[:max_plastic_strain]
    tout = collect(LinRange(0, tmax, D[:export_density]))
    sout = hardening_func(tout)

    output_folder =  isnothing(D[:export_folder]) ? pwd() : D[:export_folder]

    fname = D[:name]
    if D[:export_hardening]
        fname1 = joinpath(output_folder, fname * "_hardening.csv")
        writedlm(fname1,
                    [tout sout], 
                    delim,
                    )
    end

    if D[:export_true_stress]
        fname2 = joinpath(output_folder, fname * "_true stress.csv")             
        writedlm(fname2,
                    [true_stress.strain true_stress.stress],
                    delim,
                    )
    end

    if D[:export_plot]
        fext = ".png"
        fname3 = joinpath(output_folder, fname * "_plot" * fext)
        
        fig = Figure(size = (1200,800))
        ax = initialize_axis(fig)
        old_ax = D[:axis]
        D[:axis] = ax
        update_stress_plot(D)
        D[:axis] = old_ax

        label_status = Label(fig[2,1],tellwidth = false)
        old_st_label = D[:status_label]
        D[:status_label] = label_status
        update_status_label!(D)
        D[:status_label] = old_st_label

        screen = GLMakie.Screen()
        display(screen, fig)
        save(fname3, fig; px_per_unit)
        close(screen)
    end


    @info  "Done"
end



function draw_overview_controls!(Lay::GridLayout, D::Dict{Symbol, Any})
    @assert hasallkeys(D, [:is_true_stress, :name, :sidebar_sub_width, :figure])
    fig = D[:figure]

    cb_true_stress = Checkbox(fig, checked = false)
    
    tb_name = Textbox(fig, 
             width = D[:sidebar_sub_width],
             placeholder = "Material Name",
             boxcolor = :white)

    Lay[1, 1] = vgrid!(
                            Label(fig, "Basic properties", font=:italic, fontsize = 12, halign= :left),
                            hgrid!(cb_true_stress, Label(fig, "Data is true stress"), halign = :left),
                            # Label(fig, "Material Name", halign = :left),
                            tb_name,
                        )

    ############### BEHAVIOR ###########################################################################################

    on(cb_true_stress.checked) do val
        D[:is_true_stress] = val
        recompute_data!(D)
        update_stress_plot(D)
        update_status_label!(D)
    end

    on(tb_name.stored_string) do s
        D[:name] = s
        update_stress_plot(D)
    end

    ####################################################################################################################

                        
    return nothing

end


function draw_true_stress_controls!(Lay::GridLayout, D::Dict{Symbol, Any})

    @assert hasallkeys(D, [:resample_menu_options, :sidebar_sub_width, :figure])
    fig = D[:figure]

    resample_menu = Menu(fig, options = D[:resample_menu_options],
                                default = "Linear", width = D[:sidebar_sub_width] * 0.5)

    tb_resample = Textbox(fig, placeholder = "Enter number",
                    validator = Int, tellwidth = false,
                    width = D[:sidebar_sub_width]/4,
                    boxcolor = :white)

    btn_manual = Button(fig, label = "Manual")
    bezier_result_ref = Ref(Dict{String, Any}("status" => -1)) # holds the resampling result for bezier

    
    btn_reset = Button(fig, label = "Reset!")

    bspline_gl = GridLayout()
    
    bspline_degree = Observable{Int}(3)
    bspline_deg_gl = GridLayout()
    bspline_deg_spinner = make_spinner!(fig, bspline_deg_gl, bspline_degree, 1, 1;
                                        val_limits = (3, 10),
                                        label = "BSpline degree: ")

    bspline_num_pts = Observable{Int}(4)
    bspline_num_pts_lo = @lift($bspline_degree  + 1)
    limits_bspline_num_pts = (bspline_num_pts_lo, 100)
    bspline_num_spinner = make_spinner!(fig, bspline_deg_gl, bspline_num_pts, 1, 2;
                                        val_limits = limits_bspline_num_pts,
                                        label = "BSpline num pts: ")

    Lay[1,1] = vgrid!(
                    Label(fig, "Resample true stress", font=:italic, fontsize = 12, halign= :left),
                    Label(fig, "Resample Function", width = nothing),
                    hgrid!(resample_menu,),
                    bspline_deg_spinner, 
                    bspline_num_spinner,
                    hgrid!(Label(fig, "Resample to  "), tb_resample),
                    hgrid!(btn_manual, btn_reset),
                    ;
                    # tellheight = false,
                    halign = :left, 
                    width = D[:sidebar_sub_width],

    )

    ###### BEHAVIOR #######
    on(resample_menu.selection) do s
        # #TODO better handling of modulus update when fitting on true stress
        D[:resampler] = s
        D[:recompute_modulus] = false
        D[:resample_base_data] = true
        recompute_data!(D)
        D[:recompute_modulus] = true
        D[:resample_base_data] = false
        update_stress_plot(D)
        update_status_label!(D)
        
    end
    

    on(tb_resample.stored_string) do s
        
        D[:resample_density] = clamp(parse(Int,s), 2, 10_000)
        D[:resample_base_data] = true
        D[:recompute_modulus] = false
        recompute_data!(D)
        D[:resample_base_data] = false
        D[:recompute_modulus] = true
        
        update_stress_plot(D)
        update_status_label!(D)
    end

    on(btn_manual.clicks) do _
        # #open a new window, do the fitting,
        # #then close the launch_gui window and restart with new data
        # handle_bezier_fit(SSE, screen)
        @async begin
            @info "Opening Bezier Fitting Window"
            bezier_result_ref[]["status"] = 0 #reset the status
            bezier_result_ref[]["data"] = D[:base_data]
            CubicPiecewiseBezier.bezier_fit_fig(bezier_result_ref)
            
            @info "Closed Bezier Fitting Window"
            
            
            results = bezier_result_ref[] #get the last state
            # @info "Status", results["status"]
            bezier_fit = results["bezier fit"]
            D[:base_data] = bezier_fit
            recompute_data!(D)
            update_stress_plot(D)
            update_status_label!(D)

        end
         
    end

    on(btn_reset.clicks) do _
        D[:base_data] = deepcopy(D[:rawdata])
        recompute_data!(D)
        update_stress_plot(D)
    end

    on(bspline_degree) do deg
        deg = clamp(deg, 3, 10)
        D[:bspline_degree] = deg
        bspline_num_pts[] = clamp(bspline_num_pts[], bspline_degree[]+1, limits_bspline_num_pts[2])
        if D[:resampler] == BSplineApprox
            D[:resample_base_data] = true
            recompute_data!(D)
            D[:resample_base_data] = false
            update_stress_plot(D)
        end
        # @info "BSpline degree", D[:bspline_degree]
    end

    on(bspline_num_pts) do num
        num = clamp(num, D[:bspline_degree]+1, 100)
        D[:bspline_num_pts] = num
        if D[:resampler] == BSplineApprox
            D[:resample_base_data] = true
            recompute_data!(D)
            D[:resample_base_data] = false
            update_stress_plot(D)
        end
    end
    return nothing
end



function draw_emodulus_controls!(Lay::GridLayout, D::Dict{Symbol, Any})

    @assert hasallkeys(D, [:max_elastic_range, :sidebar_sub_width, :figure])
    fig = D[:figure]

    # sldvals = get_slider_range_values(SSE)
    sldvals = (;vmin = 10.0, vmax = 1e6, value = 50_000.0)


    sld_modulus = Slider(fig, 
                    range = LinRange(sldvals.vmin, sldvals.vmax, 1001),
                    startvalue= sldvals.value,
                    update_while_dragging =true,
                    width = D[:sidebar_sub_width],
                    )
    D[:modulus_slider] = sld_modulus

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
        Label(fig, "Modify E-modulus", font=:italic, fontsize = 12, halign= :left),
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
        if !D[:is_fixed_modulus]
            D[:e_modulus] = round(val; sigdigits =3)
            D[:recompute_modulus] = false
            recompute_data!(D)
            D[:recompute_modulus] = true
            
            update_stress_plot(D)
            tb_modulus.displayed_string = "$(round(sld_modulus.value[]; sigdigits= 3))"
            update_status_label!(D)
        end
    end

    on(sld_elastic_range.value) do val
        
        D[:max_elastic_range] = val
        lab_elastic_range_val.text = string(round(val;sigdigits = 4))
        recompute_data!(D)
        update_stress_plot(D)
        update_status_label!(D)
    end

    on(tb_modulus.stored_string) do s
        #update modulus only in the checkbox is unchecked
        if !D[:is_fixed_modulus]
            modulus = round(parse(Float64, s); sigdigits = 3)
            # update_SSE!(SSE; modulus)
            D[:e_modulus] = modulus
            D[:recompute_modulus] = false
            recompute_data!(D)
            D[:recompute_modulus] = true
            update_stress_plot(D)
            update_status_label!(D)
        end
    end

    
    on(cb_fixed.checked) do val
        D[:is_fixed_modulus] = val
    end

    on(sld_offset.value) do val
        D[:hardening_offset] = val
        D[:recompute_modulus] = false
        recompute_data!(D)
        D[:recompute_modulus] = true
        update_stress_plot(D)
        lab_offset.text = "Offset = $(round(sld_offset.value[]; sigdigits = 3))"
    end

    return sld_modulus #need this for something
end


function draw_hardening_controls!(Lay::GridLayout, D::Dict{Symbol, Any})

    @assert hasallkeys(D, [:sidebar_sub_width, :menu_hardening_fit_options, :figure])
    fig = D[:figure]

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
            Label(fig, "Change hardening curve fitting method.", font=:italic, fontsize = 12, halign= :left),
            hgrid!(Label(fig, "Method"), fit_menu),
            hgrid!(Label(fig, "Num Pts"), tb_hardening_pts),
            hgrid!(Label(fig, "Extrapolate strain to"), tb_extrapolation_strain),
            ;
            halign = :left,
            # tellheight = false, 
            width = D[:sidebar_sub_width],
    )    
    ########### BEHAVIOR #############
    on(fit_menu.selection) do s
        
        
        D[:hardening_interpolant] = s
         #only need to recalculated the interpolation fit
        D[:recompute_modulus] = false
        D[:recompute_true_stress] = false
        D[:recompute_hardening_portion] = false
        recompute_data!(D)
        D[:recompute_modulus] = true
        D[:recompute_true_stress] = true
        D[:recompute_hardening_portion] = true
        update_stress_plot(D)
        update_status_label!(D)
    end

    on(tb_extrapolation_strain.stored_string) do s
        
        # SSE["export max strain"] = clamp(parse(Float64, s), last(SSE["hardening"].strain), Inf)
        D[:max_plastic_strain] = clamp(parse(Float64, s), last(D[:hardening_portion].strain), Inf)
       
        update_stress_plot(D)
        update_status_label!(D)
        
    end

    on(tb_hardening_pts.stored_string) do s
        num_hardening_pts = Int(clamp(parse(Int, s), 2, Inf))
        # SSE["export density"] = num_hardening_pts
        D[:export_density] = num_hardening_pts
        # #block change of modulus, only hardening resampling should be changed
        # recompute_data!(D)
        update_stress_plot(D)
    end


end


function draw_export_controls!(Lay::GridLayout, D::Dict{Symbol, Any})

    @assert hasallkeys(D, [:sidebar_sub_width, :figure])
    fig = D[:figure]

    btn_export = Button(fig, label = "Export")

    cb_exp_true = Checkbox(fig, checked = true)
    cb_exp_hardening = Checkbox(fig, checked = true)
    cb_exp_plot = Checkbox(fig, checked = true)

    Lay[1,1] = vgrid!(
                    Label(fig, "Export controls", font=:italic, fontsize = 12, halign= :left),
                    hgrid!(cb_exp_true, Label(fig, "True Stress"), halign = :left),
                    hgrid!(cb_exp_hardening, Label(fig, "Hardening"),halign = :left),
                    hgrid!(cb_exp_plot, Label(fig, "Plot"), halign = :left),
                    btn_export;
                    halign = :left,
                    width =D[:sidebar_sub_width],
                    )

    ########## BEHAVIOR #################

    on(btn_export.clicks) do _
        D[:export_true_stress] = cb_exp_true.checked[]
        D[:export_hardening] = cb_exp_hardening.checked[]
        D[:export_plot]  = cb_exp_plot.checked[]
        export_data(D)

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
                        is_visible = false, #how it should look the 1st time
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
            rowsize!(main_GL, 2, Auto(true))
        else 
            sub_GL[1,1] = empty_GL
            hidden_GL[1,1] = controls
            rowsize!(main_GL, 2, Fixed(0))
        end
    end
    visibility[] = is_visible #trigger
    return nothing
end

function initialize_axis(fig)
return Axis(fig[1,1], title = "Stress Strain",
                    xlabel = "True Strain [-]",
                    ylabel = "True Stress [MPa]")
end


function make_spinner!(parent,
                    layout::GridLayout, 
                    value_observable::Observable{Int}, 
                    row::Int, 
                    col::Int;
                    btn_size = 15,
                    val_limits = (0, 40),
                    label = "Value: ",
                    )
    
    # --- 1. Create a local GridLayout to hold the three elements tightly ---
    spinner_grid = layout[row, col] = GridLayout()

    
    # Buttons: Use small, non-telling widths for compact appearance
    button_up   = Button(parent, label = "▲", tellwidth = false, tellheight = true, 
                                width = btn_size, height = btn_size)
    button_down = Button(parent, label = "▼", tellwidth = false, tellheight = true,
                                width = btn_size, height = btn_size)

    # --- 3. Internal Layout: Grouping the Buttons ---
    
    # Place the buttons into the rows of the spinner_grid (column 2)
    spinner_grid[1, 2] = button_up
    spinner_grid[2, 2] = button_down
    
    
    rowgap!(spinner_grid, 1, 0) # CRITICAL: Eliminate the spacing between the up and down buttons to make them almost touch
    
    # Ensure the button rows are sized automatically
    rowsize!(spinner_grid, 1, Auto())
    rowsize!(spinner_grid, 2, Auto())
    
   
    label_text = @lift(label * string($value_observable))
    value_label = Label(parent, label_text, 
                        tellwidth = true, 
                        tellheight = true, 
                        halign = :left,
                        padding = (5, 5, 5, 5)) # Add some padding around the text
    spinner_grid[1:2, 1] = value_label

    # Set the horizontal gap between the label and the buttons to zero
    rowgap!(spinner_grid, 0)
    
    # Configure column sizes: Col 1 (Label) takes minimum width, Col 2 (Buttons) takes minimum width
    #colsize!(spinner_grid, 1, Auto())
    colsize!(spinner_grid, 2, Auto())

    # --- 5. Wire the Logic ---

    # Increment logic
    on(button_up.clicks) do _
        v = value_observable[] + 1
        lo = isa(val_limits[1], Observable) ? val_limits[1][] : val_limits[1]
        hi = isa(val_limits[2], Observable) ? val_limits[1][] : val_limits[2]
        value_observable[] = Int(clamp(v, lo, hi))
    end

    # Decrement logic
    on(button_down.clicks) do _
        v = value_observable[] - 1
        lo = isa(val_limits[1], Observable) ? val_limits[1][] : val_limits[1]
        hi = isa(val_limits[2], Observable) ? val_limits[1][] : val_limits[2]
        value_observable[] = Int(clamp(v, lo, hi))
    end

    return spinner_grid
end