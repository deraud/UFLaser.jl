module Visualization

using Plots
using StaticArrays
using Statistics

# Import necessary modules
using ..FreeSurfaceLBM
using ..Constants
using ..Structs

# Temperature thresholds for phase transitions (gold)
const MELTING_TEMP = 1337.0
const VAPORIZATION_TEMP = 3129.0

"""
    create_basic_visualizations(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, frame_indices::Vector{Int}, save_dir::String="visualizations")

Create basic visualizations for specific frames and save them to the specified directory.
"""
function create_basic_visualizations(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, 
                                    frame_indices::Vector{Int}, save_dir::String="visualizations")
    # Create directory if it doesn't exist
    mkpath(save_dir)
    
    for frame in frame_indices
        # Temperature plot
        t_plot = visualize_temperature(ttm_lattice, frame)
        savefig(t_plot, joinpath(save_dir, "temperature_frame_$(frame).png"))
        
        # Phase plot
        p_plot = visualize_phase(vof_domain, frame)
        savefig(p_plot, joinpath(save_dir, "phase_frame_$(frame).png"))
        
        # Combined plot
        c_plot = visualize_combined(ttm_lattice, vof_domain, frame)
        savefig(c_plot, joinpath(save_dir, "combined_frame_$(frame).png"))
        
        # Velocity field
        v_plot = visualize_velocity_field(vof_domain, frame)
        savefig(v_plot, joinpath(save_dir, "velocity_frame_$(frame).png"))
        
        # Dashboard
        d_plot = create_dashboard(ttm_lattice, vof_domain, frame)
        savefig(d_plot, joinpath(save_dir, "dashboard_frame_$(frame).png"))
        
        println("Created visualizations for frame $frame")
    end
end

"""
    create_animations(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, 
                     frame_indices::Vector{Int}, save_dir::String="visualizations")

Create animations for the simulation and save them to the specified directory.
"""
function create_animations(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, 
                          frame_indices::Vector{Int}, save_dir::String="visualizations")
    # Create directory if it doesn't exist
    mkpath(save_dir)
    
    # Temperature animation
    t_anim = @animate for frame in frame_indices
        visualize_temperature(ttm_lattice, frame)
    end
    gif(t_anim, joinpath(save_dir, "temperature_animation.gif"), fps=10)
    
    # Phase animation
    p_anim = @animate for frame in frame_indices
        visualize_phase(vof_domain, frame)
    end
    gif(p_anim, joinpath(save_dir, "phase_animation.gif"), fps=10)
    
    # Combined animation
    c_anim = @animate for frame in frame_indices
        visualize_combined(ttm_lattice, vof_domain, frame)
    end
    gif(c_anim, joinpath(save_dir, "combined_animation.gif"), fps=10)
    
    # Velocity field animation
    v_anim = @animate for frame in frame_indices
        visualize_velocity_field(vof_domain, frame)
    end
    gif(v_anim, joinpath(save_dir, "velocity_animation.gif"), fps=10)
    
    # Dashboard animation
    d_anim = @animate for frame in frame_indices
        create_dashboard(ttm_lattice, vof_domain, frame)
    end
    gif(d_anim, joinpath(save_dir, "dashboard_animation.gif"), fps=10)
    
    println("Created all animations")
end

"""
    visualize_temperature(ttm_lattice::Lattice, frame::Int)

Create a heatmap of the lattice temperature.
"""
function visualize_temperature(ttm_lattice::Lattice, frame::Int)
    # Get electron and lattice temperatures
    # For real-time visualization, use current temperatures
    # For saved frames, use the temperature history if available
    Te = ttm_lattice.Tₑ
    Tl = ttm_lattice.Tₗ
    
    # If we have temperature history, use it instead
    if hasfield(typeof(ttm_lattice), :Te_history) && length(ttm_lattice.Te_history) >= frame
        Te = ttm_lattice.Te_history[frame]
    end
    
    if hasfield(typeof(ttm_lattice), :Tl_history) && length(ttm_lattice.Tl_history) >= frame
        Tl = ttm_lattice.Tl_history[frame]
    end
    
    # Calculate time in femtoseconds
    time_fs = frame * 10 * dt * 1e15  # Assuming 10 timesteps between frames
    
    # Create heatmap of lattice temperature
    p = heatmap(1:Nx, 1:Ny, transpose(Tl), 
               color=:inferno, aspect_ratio=:equal, 
               title="Lattice Temperature at t = $(round(time_fs, digits=1)) fs",
               xlabel="X position", ylabel="Y position",
               colorbar_title="Temperature (K)")
    
    # Add vaporization and melting temperature contours
    if maximum(Tl) > MELTING_TEMP
        contour!(p, 1:Nx, 1:Ny, transpose(Tl), 
                levels=[MELTING_TEMP], linewidth=2, 
                color=:cyan, alpha=0.7, linestyle=:dash,
                label="Melting Point ($(MELTING_TEMP) K)")
    end
    
    if maximum(Tl) > VAPORIZATION_TEMP
        contour!(p, 1:Nx, 1:Ny, transpose(Tl), 
                levels=[VAPORIZATION_TEMP], linewidth=2, 
                color=:white, alpha=0.7, linestyle=:solid,
                label="Vaporization Point ($(VAPORIZATION_TEMP) K)")
    end
    
    return p
end

"""
    visualize_phase(vof_domain::FreeSurfaceLBM.LBMDomain, frame::Int)

Create a visualization of material phases from the VOF simulation.
"""
function visualize_phase(vof_domain::FreeSurfaceLBM.LBMDomain, frame::Int)
    # Calculate time in femtoseconds
    time_fs = frame * 10 * dt * 1e15  # Assuming 10 timesteps between frames
    
    # Get cell types from history or current state
    cell_types = vof_domain.cell_type
    if length(vof_domain.cell_type_history) >= frame
        cell_types = vof_domain.cell_type_history[frame]
    end
    
    # Create phase map
    phases = zeros(Int, vof_domain.nx, vof_domain.ny)
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            if (cell_types[ix, iy] & FreeSurfaceLBM.CellType.FLUID) != 0
                phases[ix, iy] = 1  # Fluid
            elseif (cell_types[ix, iy] & FreeSurfaceLBM.CellType.INTERFACE) != 0
                phases[ix, iy] = 2  # Interface
            elseif (cell_types[ix, iy] & FreeSurfaceLBM.CellType.GAS) != 0
                phases[ix, iy] = 3  # Gas
            else
                phases[ix, iy] = 4  # Solid/Obstacle
            end
        end
    end
    
    # Create phase plot
    p = heatmap(1:vof_domain.nx, 1:vof_domain.ny, transpose(phases),
               color=[:blue, :cyan, :white, :gray], 
               aspect_ratio=:equal,
               title="Material Phase at t = $(round(time_fs, digits=1)) fs",
               xlabel="X position", ylabel="Y position")
    
    # Add custom legend
    plot!(p, legend=true, 
          labels=["Fluid" "Interface" "Gas" "Solid"],
          line_z=[1 2 3 4],
          linecolor=[:blue :cyan :white :gray],
          linewidth=10,
          linetype=:hline)
    
    return p
end

"""
    visualize_combined(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, frame::Int)

Create a combined visualization showing temperature with phase boundaries.
"""
function visualize_combined(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, frame::Int)
    # Calculate time in femtoseconds
    time_fs = frame * 10 * dt * 1e15  # Assuming 10 timesteps between frames
    
    # Get temperature data
    Tl = ttm_lattice.Tₗ
    if hasfield(typeof(ttm_lattice), :Tl_history) && length(ttm_lattice.Tl_history) >= frame
        Tl = ttm_lattice.Tl_history[frame]
    end
    
    # Get phase data
    cell_types = vof_domain.cell_type
    if length(vof_domain.cell_type_history) >= frame
        cell_types = vof_domain.cell_type_history[frame]
    end
    
    # Create phase boundaries (for contour overlay)
    phase_boundaries = zeros(vof_domain.nx, vof_domain.ny)
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            if (cell_types[ix, iy] & FreeSurfaceLBM.CellType.INTERFACE) != 0
                phase_boundaries[ix, iy] = 1
            end
        end
    end
    
    # Create temperature plot
    p = heatmap(1:Nx, 1:Ny, transpose(Tl), 
               color=:inferno, aspect_ratio=:equal, 
               title="Temperature with Phase Boundaries at t = $(round(time_fs, digits=1)) fs",
               xlabel="X position", ylabel="Y position",
               colorbar_title="Temperature (K)")
    
    # Overlay phase boundaries
    if sum(phase_boundaries) > 0
        # Find coordinates of interface cells for scatter plot
        interface_x = Int[]
        interface_y = Int[]
        for ix in 1:vof_domain.nx
            for iy in 1:vof_domain.ny
                if (cell_types[ix, iy] & FreeSurfaceLBM.CellType.INTERFACE) != 0
                    push!(interface_x, ix)
                    push!(interface_y, iy)
                end
            end
        end
        
        scatter!(p, interface_x, interface_y, 
                markersize=2, markershape=:circle, 
                color=:cyan, alpha=0.7, label="Interface")
    end
    
    return p
end

"""
    visualize_velocity_field(vof_domain::FreeSurfaceLBM.LBMDomain, frame::Int, scale::Float64=5.0)

Create a visualization of the velocity field with vectors.
"""
function visualize_velocity_field(vof_domain::FreeSurfaceLBM.LBMDomain, frame::Int, scale::Float64=5.0)
    # Calculate time in femtoseconds
    time_fs = frame * 10 * dt * 1e15  # Assuming 10 timesteps between frames
    
    # Get velocity data
    velocities = vof_domain.u
    if length(vof_domain.velocity_history) >= frame
        velocities = vof_domain.velocity_history[frame]
    end
    
    # Get cell types
    cell_types = vof_domain.cell_type
    if length(vof_domain.cell_type_history) >= frame
        cell_types = vof_domain.cell_type_history[frame]
    end
    
    # Extract velocity components
    vx = velocities[1, :, :]
    vy = velocities[2, :, :]
    
    # Calculate velocity magnitude
    vmag = sqrt.(vx.^2 .+ vy.^2)
    
    # Create fluid mask (fluid or interface cells)
    fluid_mask = zeros(vof_domain.nx, vof_domain.ny)
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            if (cell_types[ix, iy] & (FreeSurfaceLBM.CellType.FLUID | FreeSurfaceLBM.CellType.INTERFACE)) != 0
                fluid_mask[ix, iy] = 1.0
            end
        end
    end
    
    # Apply mask to velocity magnitude
    masked_vmag = vmag .* fluid_mask
    
    # Sample points for quiver plot (too many arrows makes it unreadable)
    skip = 5  # Display every 5th point
    x = 1:skip:vof_domain.nx
    y = 1:skip:vof_domain.ny
    
    # Create velocity magnitude plot
    p = heatmap(1:vof_domain.nx, 1:vof_domain.ny, transpose(masked_vmag),
               color=:thermal, aspect_ratio=:equal,
               title="Velocity Field at t = $(round(time_fs, digits=1)) fs",
               xlabel="X position", ylabel="Y position",
               colorbar_title="Velocity magnitude")
    
    # Add velocity vectors
    quiver!(p, repeat(x, length(y)), [y[j] for i in 1:length(x), j in 1:length(y)][:],
            quiver=(scale .* vx[1:skip:end, 1:skip:end][:], scale .* vy[1:skip:end, 1:skip:end][:]),
            color=:white, alpha=0.7)
    
    return p
end

"""
    create_dashboard(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, frame::Int)

Create a comprehensive dashboard with multiple plots for a given frame.
"""
function create_dashboard(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, frame::Int)
    # Calculate time in femtoseconds
    time_fs = frame * 10 * dt * 1e15  # Assuming 10 timesteps between frames
    
    # Get lattice and electron temperatures
    Tl = ttm_lattice.Tₗ
    Te = ttm_lattice.Tₑ
    if hasfield(typeof(ttm_lattice), :Tl_history) && length(ttm_lattice.Tl_history) >= frame
        Tl = ttm_lattice.Tl_history[frame]
    end
    if hasfield(typeof(ttm_lattice), :Te_history) && length(ttm_lattice.Te_history) >= frame
        Te = ttm_lattice.Te_history[frame]
    end
    
    # Get cell types
    cell_types = vof_domain.cell_type
    if length(vof_domain.cell_type_history) >= frame
        cell_types = vof_domain.cell_type_history[frame]
    end
    
    # Get velocities
    velocities = vof_domain.u
    if length(vof_domain.velocity_history) >= frame
        velocities = vof_domain.velocity_history[frame]
    end
    
    # Create phase map
    phases = zeros(Int, vof_domain.nx, vof_domain.ny)
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            if (cell_types[ix, iy] & FreeSurfaceLBM.CellType.FLUID) != 0
                phases[ix, iy] = 1  # Fluid
            elseif (cell_types[ix, iy] & FreeSurfaceLBM.CellType.INTERFACE) != 0
                phases[ix, iy] = 2  # Interface
            elseif (cell_types[ix, iy] & FreeSurfaceLBM.CellType.GAS) != 0
                phases[ix, iy] = 3  # Gas
            else
                phases[ix, iy] = 4  # Solid/Obstacle
            end
        end
    end
    
    # Extract velocity components
    vx = velocities[1, :, :]
    vy = velocities[2, :, :]
    vmag = sqrt.(vx.^2 .+ vy.^2)
    
    # Lattice temperature plot
    p1 = heatmap(1:Nx, 1:Ny, transpose(Tl), 
                color=:inferno, aspect_ratio=:equal, 
                title="Lattice Temperature (K)",
                colorbar_title="Temperature (K)")
    
    # Electron temperature plot
    p2 = heatmap(1:Nx, 1:Ny, transpose(Te), 
                color=:viridis, aspect_ratio=:equal, 
                title="Electron Temperature (K)",
                colorbar_title="Temperature (K)")
    
    # Material phase plot
    p3 = heatmap(1:vof_domain.nx, 1:vof_domain.ny, transpose(phases),
                color=[:blue, :cyan, :white, :gray], 
                aspect_ratio=:equal,
                title="Material Phase")
    
    # Create fluid mask for velocity
    fluid_mask = zeros(vof_domain.nx, vof_domain.ny)
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            if (cell_types[ix, iy] & (FreeSurfaceLBM.CellType.FLUID | FreeSurfaceLBM.CellType.INTERFACE)) != 0
                fluid_mask[ix, iy] = 1.0
            end
        end
    end
    
    # Velocity plot
    p4 = heatmap(1:vof_domain.nx, 1:vof_domain.ny, transpose(vmag .* fluid_mask),
                color=:thermal, aspect_ratio=:equal,
                title="Velocity Magnitude",
                colorbar_title="Velocity")
    
    # Combine plots
    dashboard = plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 1000),
                    plot_title="Laser Ablation Simulation at t = $(round(time_fs, digits=1)) fs")
    
    return dashboard
end

"""
    create_ablation_analysis(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, 
                           frame_indices::Vector{Int}, save_dir::String="visualizations")

Create analytical plots showing the evolution of key metrics over time and save them.
"""
function create_ablation_analysis(ttm_lattice::Lattice, vof_domain::FreeSurfaceLBM.LBMDomain, 
                                 frame_indices::Vector{Int}, save_dir::String="visualizations")
    # Create directory if it doesn't exist
    mkpath(save_dir)
    
    # Time points in femtoseconds
    time_fs = frame_indices .* 10 .* dt .* 1e15  # Assuming 10 timesteps between frames
    
    # Maximum temperature evolution
    max_Te = Vector{Float64}(undef, length(frame_indices))
    max_Tl = Vector{Float64}(undef, length(frame_indices))
    
    for (i, frame) in enumerate(frame_indices)
        # Get temperatures
        Te = ttm_lattice.Tₑ
        Tl = ttm_lattice.Tₗ
        if hasfield(typeof(ttm_lattice), :Te_history) && length(ttm_lattice.Te_history) >= frame
            Te = ttm_lattice.Te_history[frame]
        end
        if hasfield(typeof(ttm_lattice), :Tl_history) && length(ttm_lattice.Tl_history) >= frame
            Tl = ttm_lattice.Tl_history[frame]
        end
        
        max_Te[i] = maximum(Te)
        max_Tl[i] = maximum(Tl)
    end
    
    # Plot temperature evolution
    p1 = plot(time_fs, max_Te, label="Electron", lw=2, 
             xlabel="Time (fs)", ylabel="Maximum Temperature (K)",
             title="Maximum Temperature Evolution", 
             legend=:topright)
    plot!(p1, time_fs, max_Tl, label="Lattice", lw=2)
    # Add horizontal lines for phase transitions
    hline!(p1, [MELTING_TEMP], label="Melting Point", ls=:dash, color=:cyan)
    hline!(p1, [VAPORIZATION_TEMP], label="Vaporization Point", ls=:dash, color=:red)
    
    savefig(p1, joinpath(save_dir, "temperature_evolution.png"))
    
    # Ablation depth evolution
    ablation_depths = Vector{Float64}(undef, length(frame_indices))
    
    if hasfield(typeof(ttm_lattice), :AblationEdge)
        for (i, frame) in enumerate(frame_indices)
            # Calculate average ablation depth
            if hasfield(typeof(ttm_lattice), :AblationEdge)
                ablation_depths[i] = mean(ttm_lattice.AblationEdge) * dx * 1e9  # Convert to nm
            else
                ablation_depths[i] = 0.0
            end
        end
        
        p2 = plot(time_fs, ablation_depths, 
                xlabel="Time (fs)", ylabel="Ablation Depth (nm)",
                title="Ablation Depth Evolution", 
                legend=false, lw=2)
        
        savefig(p2, joinpath(save_dir, "ablation_depth.png"))
    end
    
    # Material phase evolution
    fluid_cells = Vector{Float64}(undef, length(frame_indices))
    gas_cells = Vector{Float64}(undef, length(frame_indices))
    
    for (i, frame) in enumerate(frame_indices)
        # Get cell types
        cell_types = vof_domain.cell_type
        if length(vof_domain.cell_type_history) >= frame
            cell_types = vof_domain.cell_type_history[frame]
        end
        
        # Count cells of each type
        fluid_count = 0
        gas_count = 0
        total_count = 0
        
        for ix in 1:vof_domain.nx
            for iy in 1:vof_domain.ny
                if (cell_types[ix, iy] & FreeSurfaceLBM.CellType.OBSTACLE) == 0  # Skip obstacle cells
                    total_count += 1
                    if (cell_types[ix, iy] & (FreeSurfaceLBM.CellType.FLUID | FreeSurfaceLBM.CellType.INTERFACE)) != 0
                        fluid_count += 1
                    elseif (cell_types[ix, iy] & FreeSurfaceLBM.CellType.GAS) != 0
                        gas_count += 1
                    end
                end
            end
        end
        
        fluid_cells[i] = fluid_count / total_count * 100  # Percentage
        gas_cells[i] = gas_count / total_count * 100  # Percentage
    end
    
    p3 = plot(time_fs, fluid_cells, label="Fluid", lw=2,
             xlabel="Time (fs)", ylabel="Material Fraction (%)",
             title="Material Phase Evolution", 
             legend=:right)
    plot!(p3, time_fs, gas_cells, label="Gas", lw=2)
    plot!(p3, time_fs, 100 .- fluid_cells .- gas_cells, label="Solid", lw=2)
    
    savefig(p3, joinpath(save_dir, "phase_evolution.png"))
    
    # Combined analysis plot
    p_combined = plot(p1, p3, layout=(2,1), size=(800, 1000),
                     plot_title="Laser Ablation Analysis")
    
    if hasfield(typeof(ttm_lattice), :AblationEdge)
        p_combined = plot(p1, p2, p3, layout=(3,1), size=(800, 1200),
                         plot_title="Laser Ablation Analysis")
    end
    
    savefig(p_combined, joinpath(save_dir, "ablation_analysis.png"))
    
    println("Created ablation analysis plots")
end

end # module Visualization