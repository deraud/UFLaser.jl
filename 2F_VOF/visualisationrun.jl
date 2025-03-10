include("./visualisation.jl")

using ..Visualization

"""
    create_visualization_dashboard(ttm_lattice, vof_domain; 
                                  save_dir="visualization_results", 
                                  frame_range=1:10:500)

Create a complete set of visualizations for TTM-VOF simulation results.

This function processes the simulation data and generates individual frames,
animations, and analysis plots, organizing them into a directory structure.

Arguments:
- `ttm_lattice`: The TTM lattice data
- `vof_domain`: The VOF domain data
- `save_dir`: Directory to save results (default: "visualization_results")
- `frame_range`: Range of frame indices to visualize (default: every 10th frame up to 500)
"""
function create_visualization_dashboard(ttm_lattice, vof_domain; 
                                       save_dir="visualization_results", 
                                       frame_range=1:10:500)
    # Create output directories
    base_dir = save_dir
    viz_dir = joinpath(base_dir, "frames")
    anim_dir = joinpath(base_dir, "animations")
    analysis_dir = joinpath(base_dir, "analysis")
    
    mkpath(base_dir)
    mkpath(viz_dir)
    mkpath(anim_dir)
    mkpath(analysis_dir)
    
    println("Starting TTM-VOF Visualization Dashboard")
    println("=======================================")
    println("Creating visualizations in: $(abspath(base_dir))")
    
    # Determine available frames
    max_ttm_frames = hasfield(typeof(ttm_lattice), :Tl_history) ? 
                     length(ttm_lattice.Tl_history) : 1
    
    max_vof_frames = length(vof_domain.cell_type_history)
    
    max_frames = min(max_ttm_frames, max_vof_frames)
    if max_frames == 0
        max_frames = 1  # At least process current state
    end
    
    # Convert frame_range to a Vector{Int} to ensure compatibility
    actual_frames = Int[]
    for f in frame_range
        if f <= max_frames
            push!(actual_frames, Int(f))
        end
    end
    
    if isempty(actual_frames)
        actual_frames = [1]
    end
    
    println("Processing $(length(actual_frames)) frames out of $(max_frames) available frames")
    
    # Generate individual frame visualizations
    println("Generating individual frame visualizations...")
    Visualization.create_basic_visualizations(ttm_lattice, vof_domain, actual_frames, viz_dir)
    
    # Generate animations
    println("Creating animations...")
    Visualization.create_animations(ttm_lattice, vof_domain, actual_frames, anim_dir)
    
    # Generate analysis plots
    println("Creating analysis plots...")
    Visualization.create_ablation_analysis(ttm_lattice, vof_domain, actual_frames, analysis_dir)
    
    println("Visualization dashboard complete!")
    println("Results saved to: $(abspath(base_dir))")
end

# Call the function with your data
create_visualization_dashboard(ttm_lattice, vof_domain)