

function data_gui(parent_screen::GLMakie.Screen , fn::AbstractString, defaults = nothing::Union{Nothing, Dict};
                sidebar_width = 300,
                tooclose = 1e-6, #tolerance for too close strain points
                )


    fig = Figure()

    ax_plot = Axis(fig[1,1])
    controls = GridLayout(fig[1,2], tellheight = false, width = sidebar_width)

    delim_list = zip(
                    (",", "TAB", "SPACE", ";",),
                    (',', '\t', ' ', ';')
                    )


    menu_delim = Menu(fig, options = delim_list, default = ",")
    

    tb_skip = Textbox(fig, 
                    validator = Int, 
                    stored_string = isnothing(defaults) ? "0" : string(defaults[:skipstart]),
                    )


    tb_strain_col = Textbox(fig,
                        validator = Int, 
                        stored_string = isnothing(defaults) ? "1" : string(defaults[:strain_col]),
                        )

    tb_stress_col = Textbox(fig, 
                        validator = Int,
                        stored_string = isnothing(defaults) ? "2" : string(defaults[:stress_col]),
    )

    tb_strain_mult = Textbox(fig, 
                        validator = Int,
                        stored_string = isnothing(defaults) ? "1.0" : string(defaults[:strain_mult]),
    )

    tb_stress_mult = Textbox(fig, 
                        validator = Int,
                        stored_string = isnothing(defaults) ? "1.0" : string(defaults[:stress_mult]),
    )


    btn_import = Button(fig, label = "Import")
    btn_done   = Button(fig, label = "Done!")
    btn_clean = Button(fig, label = "Clean!")

    controls[1,1] = vgrid!(
            hgrid!(Label(fig, "Delimiter"), menu_delim),
            hgrid!(Label(fig, "Skip rows:"), tb_skip),
            hgrid!(Label(fig, "Strain column:"), tb_strain_col),
            hgrid!(Label(fig, "Stress column:"), tb_stress_col),
            hgrid!(Label(fig, "Strain multiplier:"), tb_strain_mult),
            hgrid!(Label(fig, "Stress multiplier"), tb_stress_mult),
            hgrid!(btn_import, btn_clean),
            btn_done,
    )


    ################defaults###########################################################
    if isnothing(defaults)
        delim = Observable(',')
        skipstart = Observable(0)
        strain_col = Observable(1)
        stress_col = Observable(2)

        strain_mult = Observable(1.0)
        stress_mult = Observable(1.0)
    else
        delim = defaults[:delim]
        skipstart = defaults[:skipstart]
        strain_col = defaults[:strain_col]
        stress_col = defaults[:stress_col]

        strain_mult = defaults[:strain_mult]
        stress_mult = defaults[:stress_mult]
    end
    data = Observable(nothing::Union{Nothing, NamedTuple})

    abnormal_indices = Observable{Vector{Int}}

    ###############EVENTS##############################################################

    
    on(tb_skip.stored_string) do s
            skipstart = clamp(parse(Int, s), 0, typemax(Int))
        end

    for (tb, obs) in zip(
                        (tb_strain_col, tb_stress_col),
                        (strain_col, stress_col)
                        )
        on(tb.stored_string) do s
            #strain and stresscol must be different
            vs = clamp(parse(Int, s), 1, typemax(Int))
            obs = vs
            if obs==strain_col
                stress_col = vs + 1
            elseif obs == stress_col
                strain_col = vs + 1
            end

        end
    end
    

    for (tb, obs) in zip(
                (tb_strain_mult, tb_stress_mult),
                (strain_mult, stress_mult))
        on(tb.stored_string) do s
            obs = parse(Float64, s)
        end
    end


    on(btn_import.clicks) do _
        try
            data = read_stress_strain_data(fn;
                        delim = delim[], 
                        skipstart = skipstart[],
                        strain_col = strain_col[],
                        stress_col = stress_col[],
                        strain_multiplier = strain_mult[],
                        stress_multiplier =stress_mult[],
                        )
            abnormal_indices = find_abnormal_points(data; tooclose)
        catch
            println("Could not import!")
        end

        if !isnothing(data)
            update_data_plot!(ax_plot, data, abnormal_indices)
        else
            println("Could not import!")
        end

    end
   
    on(btn_clean.clicks) do _
        if !isnothing(data)
            clean!(data, abnormal_indices)
            empty!(abnormal_indices)
            update_data_plot!(ax_plot, data, abnormal_indices)
        end

    end

    on(btn_done.clicks) do _
        if !isnothing(data) 
            
             #close the figure
            close(screen)
            #probably the cleanest
            close(parent_screen)
            defaults[:delim] = delim[]
            defaults[:skipstart] = skipstart[]
            defaults[:strain_col] = strain_col[]
            defaults[:stress_col] = stress_col[]
            defaults[:strain_mult] = strain_mult[]
            defaults[:stress_mult] = stress_mult[]

            main(data; import_defaults = defaults)
            
        end

    end

    screen = GLMakie.Screen()
    
    display(screen, fig)
end


function update_data_plot!(ax, data, abnormal = nothing)
    @assert haskey(data, :strain) && haskey(data, :stress) "Data must have strain and stress fields!"

    empty!(ax)
    fig = ax.parent
    #delete the legend
    for (i, block) in enumerate(fig.content)
        if block isa Makie.Legend
            # 3. If found, delete it using the delete! function
            Makie.delete!(block)
        end
    end

    scatterlines!(ax, data.strain, data.stress, label = "Data", color = (:black, 0.5))
    #color abnormal points in red
    if !isnothing(abnormal)
        scatter!(ax, data.strain[abnormal], data.stress[abnormal], 
                        color = :red, 
                        marker= 'X',
                        markersize = 8)
    end

    axislegend(ax, position = :rb)

end

function clean!(data, abnormal_indices; maxiter = 10)
  
    for _ in 1:maxiter
        isempty(abnormal_indices) && break
        deleteat!(data.strain, abnormal_indices)
        deleteat!(data.stress, abnormal_indices)
    
        abnormal_indices = find_abnormal_points(data)

    end   

    return nothing

end