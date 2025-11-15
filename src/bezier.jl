using GLMakie
using Observables


#### AI generated

# --- 1. Bézier Curve Calculation ---
# A function to calculate a point on a cubic Bézier segment
# P0, P1, P2, P3 are the four control points (Points)
function cubic_bezier_point(t::Real, P0, P1, P2, P3)
    # The standard cubic Bézier formula: B(t) = (1-t)³P₀ + 3(1-t)²tP₁ + 3(1-t)t²P₂ + t³P₃
    a = 1 - t
    b = a * a
    c = b * a
    w0 = c
    w1 = 3 * b * t
    w2 = 3 * a * t^2
    w3 = t^3
    return P0 * w0 + P1 * w1 + P2 * w2 + P3 * w3
end

# Function to generate the piecewise cubic Bézier curve points
function piecewise_cubic_bezier(control_points::Vector{PT}; 
                                N_segments=50) where PT
    curve_points = PT[]
    n_points = length(control_points)

    if n_points < 4
        return curve_points
    end

    num_segments = div(n_points - 1, 3) 

    for k in 0:(num_segments - 1)
        
        idx = 3 * k + 1
        P0, P1, P2, P3 = control_points[idx], 
                        control_points[idx + 1], 
                        control_points[idx + 2], 
                        control_points[idx + 3]

        # Generate points for this segment
        for i in 0:N_segments
            t = i / N_segments
            push!(curve_points, cubic_bezier_point(t, P0, P1, P2, P3))
        end
    end

    return curve_points
end


function interactive_bezier_curve()
    # Initial control points (must be a multiple of 3 + 1 for cubic segments, e.g., 4, 7, 10...)
    initial_cpoints = Point2f[
        Point2f(0, 0), Point2f(1, 3), Point2f(3, -1), Point2f(5, 0), # First segment
        Point2f(7, 1), Point2f(8, 2), Point2f(10, 0)                  # Second segment
    ]
    
    # Observable to store and reactively update the control points
    cpoints = Observable(initial_cpoints) 
    
    # Reactive transformation: recalculate the curve whenever cpoints changes
    bezier_curve = lift(cpoints) do pts
        piecewise_cubic_bezier(pts)
    end

    # Setup the figure, axis
    fig = Figure()
    ax = Axis(fig[1, 1], title = "Interactive Piecewise Cubic Bézier Curve")

    # Plot the final Bézier curve
    lines!(ax, bezier_curve, color = :blue, linewidth = 4, label = "Bézier Curve")

    # Plot the control polygon (connecting control points)
    lines!(ax, cpoints, color = (:grey, 0.5), linestyle = :dash, label = "Control Polygon")

    # Plot the interactive control points as a scatter plot
    scatter_plot = scatter!(ax, cpoints, markersize = 15, color = :red, strokecolor = :black, strokewidth = 1, marker = :circle, label = "Control Points")

    # --- 3. Interactivity: Dragging Control Points ---
    # Find the index of the closest control point when the mouse is pressed
    dragged_index = Observable{Union{Nothing, Int}}(nothing)
    
    # Threshold for point selection (in pixels)
    const PICK_THRESHOLD = 20

    # Interaction for pressing the mouse button
    on(events(fig).mousebutton, priority = 10) do event
        # Only react to left mouse button press
        if event.button == Mouse.left && event.action == Mouse.press
            # Pick the closest control point on the scatter plot
            # `pick` returns the plot object and the index of the picked element
            plot, index = pick(ax.scene, events(ax).mouseposition[])
            
            # Check if a control point was picked
            if plot === scatter_plot && index !== nothing
                dragged_index[] = index
                # Consume the event so the default interaction (e.g., pan) doesn't run
                return Consume(true) 
            end
        end
        return Consume(false)
    end

    # Interaction for mouse movement (dragging)
    on(events(fig).mouseposition, priority = 10) do mp
        if dragged_index[] !== nothing && ispressed(fig, Mouse.left)
            # Convert mouse position (in pixels) to data coordinates
            # This is key for plotting on an Axis
            new_data_pos = Makie.mouseposition(ax.scene)
            
            # Update the specific control point's position
            current_points = cpoints[]
            current_points[dragged_index[]] = new_data_pos
            cpoints[] = current_points # Notify the Observable of the change
            
            # Consume the event to prevent other interactions from running
            return Consume(true)
        end
        return Consume(false)
    end

    # Interaction for releasing the mouse button
    on(events(fig).mousebutton, priority = 10) do event
        if event.button == Mouse.left && event.action == Mouse.release
            # Stop dragging
            dragged_index[] = nothing
            return Consume(true)
        end
        return Consume(false)
    end

    # --- 4. Functionality: Add, Remove, End Curve ---
    # Key binding for adding a control point (e.g., 'a')
    # This will place the new point at the current mouse position
    on(events(fig).keyboardbutton, priority = 10) do event
        if event.action == Keyboard.press
            current_points = cpoints[]
            
            if event.key == Keyboard.a # Add point
                # Check if we have a valid mouse position in data coordinates
                data_pos = try Makie.mouseposition(ax.scene) catch; return Consume(false) end
                
                # Check for space to add a new *cubic segment*. 
                # A new cubic segment requires 3 new points (P1, P2, P3) if starting from P0 of prev seg.
                # A simpler "add control point" might just append, or insert after the last segment end.
                # For this skeleton, we'll append a point where P3 of the last segment is.
                # A proper implementation for adding a whole cubic segment is more complex:
                # you'd likely want to add P1, P2, and P3, or calculate P1/P2 for a smooth C1/C2 continuation.
                
                # Simple append (not strictly piecewise cubic friendly unless adding P1, P2, P3)
                # Let's add a placeholder to guide the user to add *three* points at once for a new segment.
                # The user would have to hit 'a' three times to complete a segment.
                # For this skeleton, we will implement adding a point at the end, 
                # which will only form a new valid cubic segment if the user adds 3 more points.
                
                # Let's just add a point at the mouse position for flexibility. 
                # The user is responsible for keeping the count (4, 7, 10...)
                push!(current_points, data_pos)
                cpoints[] = current_points
                println("Added control point at $data_pos. Current points: $(length(current_points))")
                return Consume(true)
            
            elseif event.key == Keyboard.r # Remove last point
                if length(current_points) > 0
                    pop!(current_points)
                    cpoints[] = current_points
                    println("Removed last control point. Current points: $(length(current_points))")
                    return Consume(true)
                end
            
            elseif event.key == Keyboard.e # End / Clear curve
                 # A basic "end" is just to clear the points for a new curve
                cpoints[] = Point2f[] 
                println("Curve ended/cleared.")
                return Consume(true)
            end
        end
        return Consume(false)
    end

    # Set up the Axis limits to encompass the initial points
    autolimits!(ax)
    
    # Add a Legend (optional, but helpful)
    fig[1, 2] = Legend(fig, ax)

    # Display the figure
    display(fig)
end


function add_bezier_segment!(vertices, mousepos)
    #identify where to put new point
    
    
    tclosest = find_closest_point(vertices, mousepos) #closest point on control polygon to mouseposition
    i1, i2 = get_enclosing_segment(crv, tclosest) #indices of control vertices before and after
    V1, V2 = vertices[i1], vertices[i2]
    
    C1 = 0.25 * V1 + 0.75 * V2
    C2 = 0.50 * V1 + 0.50 * V2
    C3 = 0.75 * V1 + 0.25 * V2 

    vertices = vcat(vertices[1:i1], C1, C2, C3, vertices[i2:end])
    return nothing
end


function remove_bezier_segment!(crv, mousepos)
    vertices = crv.vertices
    length(vertices) <= 4 && return nothing #
    pt_on_curve, i1, i2 = find_closest_point(vertices, mousepos)
    

    iclosest = norm(pt_on_curve - vertices[i1])>= norm(pt_on_curve - vertices[i2]) ? i1 : i2
    if iclosest==firstindex(vertices)
        deleteat!(vertices, (1,2))
    elseif iclosest==lastindex(vertices)
        deleteat!(vertices, (iclosest-1, iclosest))
    else
        deleteat!(vertices, [iclosest-1, iclosest, iclosest+1])
    end
    return nothing
end


function find_closest_point(vertices, point) 
    
    CP = @view vertices[1:3:end]
    tclosest = -Inf
    bestfound = nothing
    seg = (-1, -1)
    for i in 2:length(CP)
        C1, C2 = CP[i-1], CP[i]
        C3, t = closest_point_to_line_segment(C1, C2, point)
        if 0<=t<=1
            return (C3, i-1, i)
        else
            if abs(t)<tclosest
                tclosest, bestfound = abs(t), C3
                seg = (i-1, i)
            end
        end

    end
    return (bestfound, seg...)
end

function closest_point_to_line_segment(L1, L2, P) #closest pt 
    PL1 = P - L1
    L21 = L2 - L1

    nL21 = norm(L21)
    nPL1 = norm(PL1)

    Q = L1 + dot(PL1, L21) * L21 / (nL21 * nPL1)
    t = (Q[1] - L2[1]) / (L1[1] - L2[1]) 
    if isinf(t)
        t = (Q[2] - L2[2]) / (L1[2] - L2[2])
    end
    return (Q, t)
end

function move_control_vertices(crv, idx, new_pos)
    vertices = crv.vertices
    old_pos = vertices[idx]
    vertices[idx] = new_pos #move the vertex to new_pos

    #case 1  - the moved vertex is a main control vertex
    # move all attached vertices by the same amount to preserve tangency
    if is_main_vertex(vertices, idx)
        if idx== firstindex(vertices)
            attached = (idx+1)
        elseif idx == lastindex(vertices)
            attached = (idx-1)
        else
            attached = (idx-1, idx+1)
        end
        dmove = new_pos - old_pos
        for i in attached
            vertices[i] += dmove
        end
        return nothing
    end

    #case 2 the moved vertex is a secondary control vertex
    #rotate around the nearest main vertex the sec. vertex across the main 
    if idx-1 == firstindex(vertices) || idx + 1 == lastindex(vertices) #nothing to do
        return nothing
    end

    if is_main_vertex(vertices, idx-1)
        center = vertices[idx-1]
        V = vertices[idx-2]
        idnext = idx-2
    else
        center = vertices[idx+1]
        V = vertices[idx+2]
        idnext = idx+2
    end
    L = norm(V-center)
    vertices[idnext] = normalize(center - vertices[idx]) * L 
    return nothing
end


function is_main_vertex(vertices, idx)
    @assert firstindex(vertices) <= idx <= lastindex(vertices)
    i = mod(idx, 3) #1 2 0 1 2 0 1 -> 1, 4, 7
    return i==1
end