

function data_gui(parent_screen::GLMakie.Screen , #GLFW screen to be able to close the parent
                fn::AbstractString;
                sidebar_width = 300,
                tooclose = 1e-6, #tolerance for too close strain points
                readahead = 10, #how many lines to read before import
                )


    DF = nothing
    try 
        DF = CSV.read(fn, DataFrame)
        size(DF,2)>= 2 || return nothing #not enough columns found...
        
    catch e
        return nothing 
    end

    initial_data = (;strain = DF[:,1], 
                            stress = DF[:,2])

    #global data store
    Imported = Dict{Symbol, Any}( :DataFrame => DF,
                                :extracted => initial_data,
                                :abnormal_indices => find_abnormal_points(initial_data; tooclose),
                                :num_cols => size(DF, 2),
                                :strain_col => 1,
                                :stress_col => 2,
                                :strain_mult => 1.0,
                                :stress_mult => 1.0,
                                )
    

    #######################LAYOUT#######################################################################################
    
    fig = Figure()

    ax_plot = Axis(fig[1,1])

    control_gl = GridLayout(fig[1,2], tellheight = false, width = sidebar_width)


    screen = GLMakie.Screen(;title = "Import Data")

    #######################PLOT######################################################################
    
    update_data_plot!(ax_plot, 
                    Imported[:extracted], 
                    Imported[:abnormal_indices])
    
    ####################################################################################################################
    
    slider_strain_col = Slider(fig, range = 1:Imported[:num_cols], startvalue = 1)
    slider_stress_col = Slider(fig, range = 1:Imported[:num_cols], startvalue = 2)

    tb_strain_mult = Textbox(fig, 
                        validator = Float64,
                        placeholder = "$(Imported[:strain_mult])",
    )

    tb_stress_mult = Textbox(fig, 
                        validator = Float64,
                        placeholder = "$(Imported[:stress_mult])",
                        )


    
    btn_done   = Button(fig, label = "Done!")
    btn_clean = Button(fig, label = "Clean!")

    control_gl[1,1] = vgrid!(
            hgrid!(Label(fig, "Strain column:"), slider_strain_col),
            hgrid!(Label(fig, "Stress column:"), slider_stress_col),
            hgrid!(Label(fig, "Strain multiplier:"), tb_strain_mult),
            hgrid!(Label(fig, "Stress multiplier"), tb_stress_mult),
            hgrid!(btn_clean, btn_done),
    )

    ###############EVENTS##############################################################

    for (i, sld) in enumerate((slider_strain_col, slider_stress_col))
        
        other_slider = i==1 ? slider_stress_col : slider_strain_col
        on(sld.value) do val
            idx_not = filter(i -> i != val, 1:Imported[:num_cols]) #the available col ids
            if other_slider.value[] == val
                set_close_to!(other_slider, first(idx_not))
            end
            strain_col = i==1 ? val : other_slider.value[]
            stress_col = i==1 ? other_slider.value[] : val
            update_Imported!(Imported; strain_col, stress_col)
            update_data_plot!(ax_plot, Imported[:extracted], Imported[:abnormal_indices])
        end
    end
    
    for (i, tb) in enumerate((tb_strain_mult, tb_stress_mult))

        
        on(tb.stored_string) do s
            val = parse(Float64, s)
            val == 0 && return nothing
            if i==1
                update_Imported!(Imported; strain_mult = val, tooclose)
            else
                update_Imported!(Imported; stress_mult = val, tooclose)
            end
            update_data_plot!(ax_plot, Imported[:extracted], Imported[:abnormal_indices])
        end

    end


    on(btn_clean.clicks) do _
        data = Imported[:extracted]
        if !isnothing(data)
            
            idx = Imported[:abnormal_indices]
            clean!(data, idx)
            Imported[:abnormal_indices] = idx
            Imported[:extracted] = data
            update_data_plot!(ax_plot, Imported[:extracted], Imported[:abnormal_indices])
        end

    end

    on(btn_done.clicks) do _
        data = Imported[:extracted]
        if !isnothing(data) 
            
             #close the figure
            close(screen)
            #probably the cleanest
            close(parent_screen)
            

            main(Imported[:extracted];
                clean_data = false,
                )
            
        end

    end

    
    
    display(screen, fig)
end


function update_data_plot!(ax, data, abnormal = nothing)
    isnothing(data) && return nothing #just in case
    @assert haskey(data, :strain) && haskey(data, :stress) "Data must have strain and stress fields!"

    empty!(ax)
    fig = ax.parent
    #delete the legend
    for (i, block) in enumerate(fig.content)
        if block isa Makie.Legend
            Makie.delete!(block)
        end
    end

    scatterlines!(ax, data.strain, data.stress, label = "Data", color = (:black, 0.5))
    #color abnormal points in red
    if !isnothing(abnormal)
        scatter!(ax, data.strain[abnormal], data.stress[abnormal], 
                        label = "Abnormal Points",
                        color = :red, 
                        marker= 'X',
                        markersize = 12)
    end

    axislegend(ax, position = :rb)
    return nothing
end

function clean!(data, abnormal_indices=nothing; maxiter = 10)
    if isnothing(abnormal_indices)
        abnormal_indices = find_abnormal_points(data)
    end
    for _ in 1:maxiter
        isempty(abnormal_indices) && break
        deleteat!(data.strain, abnormal_indices)
        deleteat!(data.stress, abnormal_indices)
    
        abnormal_indices = find_abnormal_points(data)

    end   

    return nothing

end

function update_Imported!(Imported; 
                        strain_col = nothing,
                        stress_col = nothing,
                        strain_mult = nothing,
                        stress_mult = nothing,
                        tooclose = 1e-6,
                    )
@assert hasallkeys(Imported, (:DataFrame, :extracted, :strain_col, :stress_col, :strain_mult, :stress_mult, :abnormal_indices))
@assert hasallkeys(Imported[:extracted], (:strain, :stress))

strain_col = isnothing(strain_col) ? Imported[:strain_col] : strain_col
stress_col = isnothing(stress_col) ? Imported[:stress_col] : stress_col
strain_mult = isnothing(strain_mult) ? Imported[:strain_mult] : strain_mult
stress_mult = isnothing(stress_mult) ? Imported[:stress_mult] : stress_mult

data = Imported[:extracted]
strain = Imported[:DataFrame][:, strain_col] .* strain_mult
stress = Imported[:DataFrame][:, stress_col] .* stress_mult
Imported[:extracted] = (;strain, stress)
Imported[:abnormal_indices] = find_abnormal_points(data;tooclose)

return nothing
end

function hasallkeys(container, keys)
    for key in keys
        !haskey(container, key) && return false 
    end
    return true
end