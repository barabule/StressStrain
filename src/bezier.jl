using GLMakie
using Observables

# --- 1. Bézier Curve Calculation ---
# A function to calculate a point on a cubic Bézier segment
# P0, P1, P2, P3 are the four control points (Points)
function cubic_bezier_point(t::Real, P0, P1, P2, P3)
    # The standard cubic Bézier formula: B(t) = (1-t)³P₀ + 3(1-t)²tP₁ + 3(1-t)t²P₂ + t³P₃
    w0 = (1 - t)^3
    w1 = 3 * (1 - t)^2 * t
    w2 = 3 * (1 - t) * t^2
    w3 = t^3
    return P0 * w0 + P1 * w1 + P2 * w2 + P3 * w3
end

# Function to generate the piecewise cubic Bézier curve points
function piecewise_cubic_bezier(control_points::Vector{Point2f}, N_segments=50)
    curve_points = Point2f[]
    n_points = length(control_points)

    # A piecewise cubic Bézier requires 4 points (P₀, P₁, P₂, P₃) for each segment.
    # The total number of points must be 4 + 3k where k is the number of extra segments
    # (i.e., n_points = 4, 7, 10, ...).
    # The start/end of a segment are the end points (P0, P3). The two middle points (P1, P2) are the
    # control handles.

    # Check for enough points to form at least one cubic segment (4 points)
    if n_points < 4
        return curve_points
    end

    # The number of segments is (n_points - 1) / 3. The control point indices for segment 'k'
    # will be 3k, 3k+1, 3k+2, 3k+3 (assuming 0-indexing for segments).
    num_segments = div(n_points - 1, 3) # integer division

    for k in 0:(num_segments - 1)
        # 1-based indexing for control_points array
        idx = 3 * k + 1
        P0, P1, P2, P3 = control_points[idx], control_points[idx + 1], control_points[idx + 2], control_points[idx + 3]

        # Generate points for this segment
        for i in 0:N_segments
            t = i / N_segments
            push!(curve_points, cubic_bezier_point(t, P0, P1, P2, P3))
        end
    end

    return curve_points
end

# --- 2. Setup Plot and Observables ---
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

# To run the function:
# interactive_bezier_curve()