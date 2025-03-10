using Plots
using Random

# Constants for the CST parametrization
const N1 = 0.5  # Class parameter for airfoil leading edge
const N2 = 1.0  # Class parameter for airfoil trailing edge

# Structure to hold blade section data
struct BladeSection
    radial_position::Float64  # r/R - normalized radial position
    chord::Float64            # chord length in meters
    twist::Float64            # twist angle in degrees
    airfoil_params::Vector{Float64}  # CST parameters [upper_params; lower_params]
end

# Function to create airfoil using CST parametrization
function create_airfoil_cst(upper_params::Vector{Float64}, lower_params::Vector{Float64}, n_points::Int=100)
    # Normalized x-coordinates from 0 to 1
    x = LinRange(0, 1, n_points)
    
    # Calculate class function C(x) = x^N1 * (1-x)^N2
    C_x = x .^ N1 .* (1 .- x) .^ N2
    
    # Shape function S(x) for upper and lower surfaces
    function shape_function(x, params)
        n = length(params)
        S = zeros(length(x))
        for i in 1:n
            # Bernstein polynomial basis
            K = binomial(n - 1, i - 1)
            S .+= params[i] .* K .* x .^ (i - 1) .* (1 .- x) .^ (n - i)
        end
        return S
    end
    
    upper_surface = C_x .* shape_function(x, upper_params)
    lower_surface = -C_x .* shape_function(x, lower_params)
    
    # Return x, y coordinates for the airfoil
    return x, upper_surface, lower_surface
end

function generate_random_blade_sections(n_sections::Int=5)
    # Radial positions for each section (0 at hub, 1 at tip)
    radial_positions = LinRange(0.2, 1.0, n_sections)
    
    # Store random sections
    Random.seed!(42)  # For reproducibility
    
    sections = BladeSection[]
    
    for (i, r) in enumerate(radial_positions)
        # Generate random CST parameters
        # Upper surface parameters typically positive
        upper_params = 0.2 .+ 0.3 .* rand(4)
        # Lower surface parameters typically negative
        lower_params = -0.1 .- 0.2 .* rand(4)
        
        # Calculate chord and twist (simplified distribution)
        chord = 0.2 * (1.0 - 0.6 * r)  # Example chord distribution
        twist = 25.0 * (1.0 - r)       # Example twist distribution
        
        # Create and store the section
        params = vcat(upper_params, lower_params)
        section = BladeSection(r, chord, twist, params)
        push!(sections, section)
    end
    
    return sections
end

function visualize_blade(sections::Vector{BladeSection})
    # Create a plot for all airfoil sections
    p1 = plot(
        title = "Blade Airfoil Sections",
        xlabel = "x/c",
        ylabel = "y/c",
        aspect_ratio = :equal,
        legend = :topleft,
        size = (800, 600),
        framestyle = :box
    )
    
    # Colors for different sections
    colors = [:blue, :red, :green, :purple, :orange]
    
    # Plot each airfoil section
    for (i, section) in enumerate(sections)
        upper_params = section.airfoil_params[1:4]
        lower_params = section.airfoil_params[5:8]
        
        x, y_upper, y_lower = create_airfoil_cst(upper_params, lower_params)
        
        color = colors[mod1(i, length(colors))]
        
        # Plot with appropriate scaling and offset based on radial position
        plot!(p1, x, y_upper, label="Section $i (r/R = $(round(section.radial_position, digits=2)))", 
              linewidth=2, color=color)
        plot!(p1, x, y_lower, label="", linewidth=2, color=color)
        
        # Connect trailing edge
        plot!(p1, [x[end], x[end]], [y_upper[end], y_lower[end]], linewidth=1, color=color, label="")
    end
    
    # Create 3D visualization of the full blade
    p2 = plot(
        title = "3D Blade Visualization",
        xlabel = "Spanwise (r/R)",
        ylabel = "Chordwise (x/c)",
        zlabel = "y/c",
        camera = (30, 30),
        size = (800, 600),
        legend = false
    )
    
    # Plot 3D representation with triangles to create surfaces
    for (i, section) in enumerate(sections)
        upper_params = section.airfoil_params[1:4]
        lower_params = section.airfoil_params[5:8]
        
        x, y_upper, y_lower = create_airfoil_cst(upper_params, lower_params)
        
        # Scale by chord
        x_scaled = x .* section.chord
        y_upper_scaled = y_upper .* section.chord
        y_lower_scaled = y_lower .* section.chord
        
        # Apply twist
        twist_rad = deg2rad(section.twist)
        cos_t = cos(twist_rad)
        sin_t = sin(twist_rad)
        
        x_upper_twisted = x_scaled .* cos_t - y_upper_scaled .* sin_t
        y_upper_twisted = x_scaled .* sin_t + y_upper_scaled .* cos_t
        
        x_lower_twisted = x_scaled .* cos_t - y_lower_scaled .* sin_t
        y_lower_twisted = x_scaled .* sin_t + y_lower_scaled .* cos_t
        
        # Create z coordinates based on radial position
        z_upper = fill(section.radial_position, length(x))
        z_lower = z_upper
        
        color = colors[mod1(i, length(colors))]
        
        # Plot in 3D
        plot!(p2, z_upper, x_upper_twisted, y_upper_twisted, linecolor=color, linewidth=2)
        plot!(p2, z_lower, x_lower_twisted, y_lower_twisted, linecolor=color, linewidth=2)
        
        # Add a connection at trailing edge
        plot!(p2, [z_upper[end], z_lower[end]], 
                 [x_upper_twisted[end], x_lower_twisted[end]], 
                 [y_upper_twisted[end], y_lower_twisted[end]], 
                 linecolor=color, linewidth=1)
        
        # Connect to next section if not the last one
        if i < length(sections) && i < length(sections)
            next_section = sections[i+1]
            next_upper_params = next_section.airfoil_params[1:4]
            next_lower_params = next_section.airfoil_params[5:8]
            
            _, next_y_upper, next_y_lower = create_airfoil_cst(next_upper_params, next_lower_params)
            
            # Scale by chord
            next_x_scaled = x .* next_section.chord
            next_y_upper_scaled = next_y_upper .* next_section.chord
            next_y_lower_scaled = next_y_lower .* next_section.chord
            
            # Apply twist
            next_twist_rad = deg2rad(next_section.twist)
            next_cos_t = cos(next_twist_rad)
            next_sin_t = sin(next_twist_rad)
            
            next_x_upper_twisted = next_x_scaled .* next_cos_t - next_y_upper_scaled .* next_sin_t
            next_y_upper_twisted = next_x_scaled .* next_sin_t + next_y_upper_scaled .* next_cos_t
            
            # Sample some points to connect the sections (e.g., leading edge, trailing edge, mid-point)
            sample_points = [1, div(length(x), 2), length(x)]
            
            for j in sample_points
                plot!(p2, [z_upper[j], fill(next_section.radial_position, 1)[1]], 
                         [x_upper_twisted[j], next_x_upper_twisted[j]], 
                         [y_upper_twisted[j], next_y_upper_twisted[j]], 
                         linecolor=:gray, linewidth=1, linestyle=:dash)
            end
        end
    end
    
    # Create a 3D view showing the full turbine with 3 blades
    p3 = plot(
        title = "3-Blade Turbine (Top View)",
        aspect_ratio = :equal,
        legend = false,
        size = (600, 600),
        framestyle = :box
    )
    
    # Plot top view of complete turbine
    radius = 1.0  # Normalized radius
    
    # Draw hub
    hub_radius = 0.2 * radius
    θ = LinRange(0, 2π, 100)
    plot!(p3, hub_radius .* cos.(θ), hub_radius .* sin.(θ), linecolor=:black, fillcolor=:gray, fill=true)
    
    # Plot the three blades at 120 degree intervals
    for blade_angle in [0, 120, 240]
        blade_angle_rad = deg2rad(blade_angle)
        blade_cos = cos(blade_angle_rad)
        blade_sin = sin(blade_angle_rad)
        
        for section in sections
            # Get chord length and radial position
            r = section.radial_position * radius
            chord = section.chord
            
            # Calculate points for a simplified blade section (just a line for top view)
            x_center = r * blade_cos
            y_center = r * blade_sin
            
            # Rotate the blade section according to the blade angle
            rotation = blade_angle_rad + π/2  # Perpendicular to radial direction
            
            # Draw a simple line representing the chord at this section
            half_chord = chord / 2
            x1 = x_center - half_chord * cos(rotation)
            y1 = y_center - half_chord * sin(rotation)
            x2 = x_center + half_chord * cos(rotation)
            y2 = y_center + half_chord * sin(rotation)
            
            plot!(p3, [x1, x2], [y1, y2], linewidth=2, color=:blue)
        end
    end
    
    # Add turbine radius circle
    plot!(p3, radius .* cos.(θ), radius .* sin.(θ), linecolor=:black, linestyle=:dash)
    
    # Set plot limits
    plot!(p3, xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
    
    # Display all plots
    return plot(p1, p2, p3, layout = (2, 2), size = (1200, 1000))
end

# Generate random sections and visualize
random_sections = generate_random_blade_sections(5)
plot_result = visualize_blade(random_sections)

# Display the plot
display(plot_result)

# Print the parameters used
println("Random blade section parameters:")
for (i, section) in enumerate(random_sections)
    println("Section $i (r/R = $(section.radial_position)):")
    println("  Chord = $(section.chord)m")
    println("  Twist = $(section.twist)°")
    println("  Upper params: ", section.airfoil_params[1:4])
    println("  Lower params: ", section.airfoil_params[5:8])
    println()
end