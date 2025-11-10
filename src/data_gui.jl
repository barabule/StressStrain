

function data_gui(parent_screen::GLMakie.Screen , fn::AbstractString;
                sidebar_width = 300)


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
                    )


    tb_strain_col = Textbox(fig,
                        validator = Int)

    tb_stress_col = Textbox(fig, 
                        validator = Int)

    tb_strain_mult = Textbox(fig, 
                        validator = Int)

    tb_stress_mult = Textbox(fig, 
                        validator = Int)


    btn_import = Button(fig, label = "Import")
    btn_done   = Button(fig, label = "Done!")

    controls[1,1] = vgrid!(
            hgrid!(Label(fig, "Delimiter"), menu_delim),
            hgrid!(Label(fig, "Skip rows:"), tb_skip),
            hgrid!(Label(fig, "Strain column:"), tb_strain_col),
            hgrid!(Label(fig, "Stress column:"), tb_stress_col),
            hgrid!(Label(fig, "Strain multiplier:"), tb_strain_mult),
            hgrid!(Label(fig, "Stress multiplier"), tb_stress_mult),
            btn_import,
            btn_done,
    )


    ################defaults###########################################################

    delim = Observable(',')
    skipstart = Observable(0)
    strain_col = Observable(1)
    stress_col = Observable(2)

    strain_mult = Observable(1.0)
    stress_mult = Observable(1.0)

    data = Observable(nothing::Union{Nothing, NamedTuple})
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
        
        data = read_stress_strain_data(fn;
                    delim = delim[], 
                    skipstart = skipstart[],
                    strain_col = strain_col[],
                    stress_col = stress_col[],
                    strain_multiplier = strain_mult[],
                    stress_multiplier =stress_mult[],
                    )
        if !isnothing(data)
            update_data_plot!(ax_plot, data)
        else
            println("Could not import!")
        end

    end
   

    on(btn_done.clicks) do _
        if !isnothing(data) 
            
             #close the figure
            close(screen)
            #probably the cleanest
            close(parent_screen)
            main(data)
        end

    end

    screen = GLMakie.Screen()
    
    display(screen, fig)
end


function update_data_plot!(ax, data)
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

    scatterlines!(ax, data.strain, data.stress, label = "Data")


end

