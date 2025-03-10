module FreeSurfaceLBM

using LinearAlgebra
using Plots

# Define cell types as bitflags (equivalent to CT_Enum in original)
module CellType
    const FLUID = 2^0          # 1
    const INTERFACE = 2^1      # 2
    const GAS = 2^2            # 4
    const OBSTACLE = 2^3       # 8
    const INLET = 2^4          # 16
    const OUTLET = 2^5         # 32
    
    const NO_FLUID_NEIGH = 2^6  # 64
    const NO_EMPTY_NEIGH = 2^7  # 128
    const NO_IFACE_NEIGH = 2^8  # 256
    
    const TO_FLUID = 2^9       # 512
    const TO_GAS = 2^10        # 1024
end

# Define D2Q9 lattice constants
const v = [
    [0, 0],   # 0: center
    [0, -1],  # 1: north
    [0, 1],   # 2: south
    [-1, 0],  # 3: west
    [-1, -1], # 4: northwest
    [-1, 1],  # 5: southwest
    [1, 0],   # 6: east
    [1, -1],  # 7: northeast
    [1, 1]    # 8: southeast
]

# Weights for the D2Q9 model
const t = [
    4/9,   # center
    1/9,   # north, south, east, west
    1/9,
    1/9,
    1/36,  # diagonals
    1/36,
    1/9,
    1/36,
    1/36
]

# Inverse direction mapping
function calculate_v_inv()
    v_inv = zeros(Int, 9)
    for i in 1:9
        for j in 1:9
            if v[i][1] == -v[j][1] && v[i][2] == -v[j][2]
                v_inv[i] = j
                break
            end
        end
    end
    return v_inv
end

const v_inv = calculate_v_inv()

@kwdef mutable struct LBMDomain
    nx::Int
    ny::Int
    omega::Float64 = 1.0
    gravity::Vector{Float64} = [0.0, -0.1]
    
    # Distribution functions and macroscopic quantities
    fin::Array{Float64, 3} = zeros(Float64, 9, nx, ny)
    fout::Array{Float64, 3} = zeros(Float64, 9, nx, ny)
    equi::Array{Float64, 3} = zeros(Float64, 9, nx, ny)
    fdist::Array{Float64, 3} = zeros(Float64, 9, nx, ny)
    inlet::Array{Float64, 3} = zeros(Float64, 2, nx, ny)
    u::Array{Float64, 3} = zeros(Float64, 2, nx, ny)
    rho::Matrix{Float64} = zeros(Float64, nx, ny)
    mass::Matrix{Float64} = zeros(Float64, nx, ny)
    cell_type::Matrix{Int} = fill(CellType.GAS, nx, ny)
    
    # Storage for previous timestep
    rho_prev::Matrix{Float64} = zeros(Float64, nx, ny)
    mass_prev::Matrix{Float64} = zeros(Float64, nx, ny)
    u_prev::Array{Float64, 3} = zeros(Float64, 2, nx, ny)
    cell_type_prev::Matrix{Int} = fill(CellType.GAS, nx, ny)
    
    # History for visualization
    velocity_history::Vector{Array{Float64, 3}} = Array{Float64, 3}[]
    cell_type_history::Vector{Matrix{Int}} = Matrix{Int}[]
end

# Set equilibrium distribution function for a cell
function set_equi!(domain::LBMDomain, ix::Int, iy::Int, rho_val::Float64, u_vec::Vector{Float64})
    usqr = 1.5 * (u_vec[1]^2 + u_vec[2]^2)
    for n in 1:9
        vu = 3.0 * (v[n][1] * u_vec[1] + v[n][2] * u_vec[2])
        domain.equi[n, ix, iy] = rho_val * t[n] * (1.0 + vu + 0.5 * vu^2 - usqr)
    end
end

# Update density and velocity for a cell
function update_rho_u!(domain::LBMDomain, ix::Int, iy::Int)
    domain.rho[ix, iy] = 0.0
    domain.u[1, ix, iy] = 0.0
    domain.u[2, ix, iy] = 0.0
    
    for n in 1:9
        domain.rho[ix, iy] += domain.fin[n, ix, iy]
        domain.u[1, ix, iy] += v[n][1] * domain.fin[n, ix, iy]
        domain.u[2, ix, iy] += v[n][2] * domain.fin[n, ix, iy]
    end
    
    # Divide velocity by density
    domain.u[:, ix, iy] ./= max(domain.rho[ix, iy], 1e-10)
    
    # Velocity limiting for stability
    vel = sqrt(domain.u[1, ix, iy]^2 + domain.u[2, ix, iy]^2)
    maxvel = sqrt(2/3)
    if vel > maxvel
        domain.u[:, ix, iy] .*= maxvel / vel
    end
end

# Collision step of LBM
function collide!(domain::LBMDomain, ix::Int, iy::Int)
    for n in 1:9
        domain.fout[n, ix, iy] = domain.fin[n, ix, iy] + 
                                 domain.omega * (domain.equi[n, ix, iy] - domain.fin[n, ix, iy])
        
        # Add gravitational forces
        grav_temp = v[n][1] * domain.gravity[1] + v[n][2] * domain.gravity[2]
        domain.fout[n, ix, iy] -= domain.rho[ix, iy] * t[n] * grav_temp
    end
end

# Get the fraction a cell is filled with fluid
function get_epsilon(cell_type::Int, rho::Float64, mass::Float64)
    if (cell_type & CellType.FLUID != 0) || (cell_type & CellType.OBSTACLE != 0)
        return 1.0
    elseif (cell_type & CellType.GAS) == CellType.GAS
        return 0.0
    else
        if rho > 0.0
            # Calculate and clip epsilon
            epsilon = mass / rho
            return clamp(epsilon, 0.0, 1.0)
        end
        return 0.5
    end
end

# Get the surface normal at an interface cell
function get_normal(domain::LBMDomain, ix::Int, iy::Int)
    if ix <= 1 || ix >= domain.nx || iy <= 1 || iy >= domain.ny
        return (0.0, 0.0)
    end
    
    x = 0.5 * (get_epsilon(domain.cell_type_prev[ix-1, iy], domain.rho_prev[ix-1, iy], domain.mass_prev[ix-1, iy]) -
               get_epsilon(domain.cell_type_prev[ix+1, iy], domain.rho_prev[ix+1, iy], domain.mass_prev[ix+1, iy]))
    y = 0.5 * (get_epsilon(domain.cell_type_prev[ix, iy-1], domain.rho_prev[ix, iy-1], domain.mass_prev[ix, iy-1]) -
               get_epsilon(domain.cell_type_prev[ix, iy+1], domain.rho_prev[ix, iy+1], domain.mass_prev[ix, iy+1]))
    return (x, y)
end

# Calculate mass exchange between interface cells
function interface_mass_exchange(cell_type_self::Int, cell_type_nei::Int, fout_self::Float64, fout_nei::Float64)
    if (cell_type_self & CellType.NO_FLUID_NEIGH) != 0
        if (cell_type_nei & CellType.NO_FLUID_NEIGH) != 0
            return fout_nei - fout_self
        else
            return -fout_self
        end
    elseif (cell_type_self & CellType.NO_EMPTY_NEIGH) != 0
        if (cell_type_nei & CellType.NO_EMPTY_NEIGH) != 0
            return fout_nei - fout_self
        else
            return fout_nei
        end
    else
        if (cell_type_nei & CellType.NO_FLUID_NEIGH) != 0
            return fout_nei
        elseif (cell_type_nei & CellType.NO_EMPTY_NEIGH) != 0
            return fout_self
        else
            return fout_nei - fout_self
        end
    end
end

# Initialize surrounding values for a new interface cell
function average_surround!(domain::LBMDomain, ix::Int, iy::Int)
    domain.mass[ix, iy] = 0.0
    domain.rho[ix, iy] = 0.0
    domain.u[:, ix, iy] .= 0.0
    
    count = 0
    for n in 2:9  # Skip the center direction
        ix_next = ix + v[n][1]
        iy_next = iy + v[n][2]
        
        if 1 <= ix_next <= domain.nx && 1 <= iy_next <= domain.ny
            cell = domain.cell_type_prev[ix_next, iy_next]
            if (cell & (CellType.FLUID | CellType.INTERFACE)) != 0
                count += 1
                domain.rho[ix, iy] += domain.rho_prev[ix_next, iy_next]
                domain.u[1, ix, iy] += domain.u_prev[1, ix_next, iy_next]
                domain.u[2, ix, iy] += domain.u_prev[2, ix_next, iy_next]
            end
        end
    end
    
    if count > 0
        domain.rho[ix, iy] /= count
        domain.u[:, ix, iy] ./= count
    end
    
    set_equi!(domain, ix, iy, domain.rho[ix, iy], domain.u[:, ix, iy])
end

# LBM collision step for the whole domain
function lbm_collision!(domain::LBMDomain)
    for ix in 1:domain.nx
        for iy in 1:domain.ny
            # Skip empty cells and obstacles
            if (domain.cell_type[ix, iy] & (CellType.GAS | CellType.OBSTACLE)) != 0
                continue
            end
            
            # Regular fluid or interface cells
            update_rho_u!(domain, ix, iy)
            
            # Calculate equilibrium distribution
            set_equi!(domain, ix, iy, domain.rho[ix, iy], domain.u[:, ix, iy])
            
            # Perform collision
            collide!(domain, ix, iy)
        end
    end
end

# LBM streaming step for the whole domain
function lbm_streaming!(domain::LBMDomain)
    # Copy previous state values
    domain.rho_prev .= domain.rho
    domain.mass_prev .= domain.mass
    domain.u_prev .= domain.u
    domain.cell_type_prev .= domain.cell_type
    
    for ix in 1:domain.nx
        for iy in 1:domain.ny
            # Skip obstacles and gas cells
            if (domain.cell_type[ix, iy] & (CellType.OBSTACLE | CellType.GAS)) != 0
                continue
            end
            
            # Handle fluid cells
            if (domain.cell_type[ix, iy] & (CellType.FLUID | CellType.INLET | CellType.OUTLET)) != 0
                for n in 1:9
                    ix_next = ix - v[n][1]
                    iy_next = iy - v[n][2]
                    
                    # Skip out of bounds
                    if ix_next < 1 || ix_next > domain.nx || iy_next < 1 || iy_next > domain.ny
                        continue
                    end
                    
                    if (domain.cell_type[ix_next, iy_next] & (CellType.FLUID | CellType.INTERFACE | CellType.INLET | CellType.OUTLET)) != 0
                        # Regular streaming
                        domain.fin[n, ix, iy] = domain.fout[n, ix_next, iy_next]
                        # Mass exchange
                        domain.mass[ix, iy] += domain.fout[n, ix_next, iy_next] - domain.fout[v_inv[n], ix, iy]
                    else  # Obstacle
                        # Bounce-back (no-slip boundary)
                        domain.fin[n, ix, iy] = domain.fout[v_inv[n], ix, iy]
                    end
                end
                
            # Handle interface cells
            elseif (domain.cell_type[ix, iy] & CellType.INTERFACE) != 0
                # Get the fluid fraction
                epsilon = get_epsilon(domain.cell_type[ix, iy], domain.rho_prev[ix, iy], domain.mass_prev[ix, iy])
                set_equi!(domain, ix, iy, 1.0, domain.u[:, ix, iy])
                
                for n in 2:9  # Skip center direction
                    ix_next = ix - v[n][1]
                    iy_next = iy - v[n][2]
                    
                    # Skip out of bounds
                    if ix_next < 1 || ix_next > domain.nx || iy_next < 1 || iy_next > domain.ny
                        continue
                    end
                    
                    if (domain.cell_type[ix_next, iy_next] & CellType.FLUID) != 0
                        # Streaming from fluid cell
                        domain.fin[n, ix, iy] = domain.fout[n, ix_next, iy_next]
                        # Mass exchange
                        domain.mass[ix, iy] += domain.fout[n, ix_next, iy_next] - domain.fout[v_inv[n], ix, iy]
                    elseif (domain.cell_type[ix_next, iy_next] & CellType.INTERFACE) != 0
                        # Streaming from another interface cell
                        domain.fin[n, ix, iy] = domain.fout[n, ix_next, iy_next]
                        # Mass exchange with weighted epsilon
                        epsilon_nei = get_epsilon(domain.cell_type[ix_next, iy_next], 
                                                 domain.rho_prev[ix_next, iy_next], 
                                                 domain.mass_prev[ix_next, iy_next])
                        domain.mass[ix, iy] += interface_mass_exchange(
                            domain.cell_type[ix, iy], 
                            domain.cell_type[ix_next, iy_next],
                            domain.fout[v_inv[n], ix, iy], 
                            domain.fout[n, ix_next, iy_next]
                        ) * (epsilon + epsilon_nei) * 0.5
                    elseif (domain.cell_type[ix_next, iy_next] & CellType.GAS) != 0
                        # Streaming from gas cell (use equilibrium)
                        domain.fin[n, ix, iy] = domain.equi[n, ix, iy] + domain.equi[v_inv[n], ix, iy] - domain.fout[v_inv[n], ix, iy]
                    else  # Obstacle
                        # Bounce-back (no-slip boundary)
                        domain.fin[n, ix, iy] = domain.fout[v_inv[n], ix, iy]
                    end
                end
                
                # Surface normal correction
                normal = get_normal(domain, ix, iy)
                for n in 2:9
                    if (normal[1] * v[v_inv[n]][1] + normal[2] * v[v_inv[n]][2]) > 0
                        domain.fin[n, ix, iy] = domain.equi[n, ix, iy] + domain.equi[v_inv[n], ix, iy] - domain.fout[v_inv[n], ix, iy]
                    end
                end
            end
            
            # Update density and velocity
            update_rho_u!(domain, ix, iy)
            if (domain.cell_type[ix, iy] & CellType.FLUID) != 0
                domain.rho[ix, iy] = domain.mass[ix, iy]
            end
        end
    end
end

# Update cell types after streaming
function update_types!(domain::LBMDomain)
    fill_offset = 0.003
    lonely_thresh = 0.1
    
    # First pass: mark cells to change
    for ix in 1:domain.nx
        for iy in 1:domain.ny
            if (domain.cell_type[ix, iy] & CellType.INTERFACE) != 0
                if domain.mass[ix, iy] > (1 + fill_offset) * domain.rho[ix, iy] || 
                   (domain.mass[ix, iy] >= (1 - lonely_thresh) * domain.rho[ix, iy] && 
                    (domain.cell_type[ix, iy] & CellType.NO_EMPTY_NEIGH) != 0)
                    domain.cell_type[ix, iy] = CellType.TO_FLUID
                elseif domain.mass[ix, iy] < -fill_offset * domain.rho[ix, iy] || 
                       (domain.mass[ix, iy] <= lonely_thresh * domain.rho[ix, iy] && 
                        (domain.cell_type[ix, iy] & CellType.NO_FLUID_NEIGH) != 0) || 
                       ((domain.cell_type[ix, iy] & (CellType.NO_IFACE_NEIGH | CellType.NO_FLUID_NEIGH)) != 0)
                    domain.cell_type[ix, iy] = CellType.TO_GAS
                end
            end
            
            # Remove neighborhood flags
            domain.cell_type[ix, iy] &= ~(CellType.NO_FLUID_NEIGH + CellType.NO_EMPTY_NEIGH + CellType.NO_IFACE_NEIGH)
        end
    end
    
    # Second pass: interface -> fluid
    for ix in 1:domain.nx
        for iy in 1:domain.ny
            if (domain.cell_type[ix, iy] & CellType.TO_FLUID) != 0
                for n in 2:9
                    ix_next = ix - v[n][1]
                    iy_next = iy - v[n][2]
                    
                    if ix_next < 1 || ix_next > domain.nx || iy_next < 1 || iy_next > domain.ny
                        continue
                    end
                    
                    if (domain.cell_type[ix_next, iy_next] & CellType.GAS) != 0
                        domain.cell_type[ix_next, iy_next] = CellType.INTERFACE
                        average_surround!(domain, ix_next, iy_next)
                    end
                end
            end
        end
    end
    
    # Third pass: interface -> gas
    for ix in 1:domain.nx
        for iy in 1:domain.ny
            if (domain.cell_type[ix, iy] & CellType.TO_GAS) != 0
                for n in 2:9
                    ix_next = ix - v[n][1]
                    iy_next = iy - v[n][2]
                    
                    if ix_next < 1 || ix_next > domain.nx || iy_next < 1 || iy_next > domain.ny
                        continue
                    end
                    
                    if (domain.cell_type[ix_next, iy_next] & CellType.FLUID) != 0
                        domain.cell_type[ix_next, iy_next] = CellType.INTERFACE
                    end
                end
            end
        end
    end
    
    # Fourth pass: distribute excess mass
    for ix in 1:domain.nx
        for iy in 1:domain.ny
            if (domain.cell_type[ix, iy] & CellType.OBSTACLE) != 0
                continue
            end
            
            normal = get_normal(domain, ix, iy)
            mex = 0.0
            
            if (domain.cell_type[ix, iy] & CellType.TO_FLUID) != 0
                mex = domain.mass[ix, iy] - domain.rho[ix, iy]
                domain.mass[ix, iy] = domain.rho[ix, iy]
            elseif (domain.cell_type[ix, iy] & CellType.TO_GAS) != 0
                mex = domain.mass[ix, iy]
                normal = (-normal[1], -normal[2])
                domain.mass[ix, iy] = 0.0
            else
                continue
            end
            
            eta = zeros(Float64, 9)
            isIF = zeros(Int, 9)
            eta_total = 0.0
            IF_total = 0
            
            for n in 2:9
                ix_next = ix + v[n][1]
                iy_next = iy + v[n][2]
                
                if ix_next < 1 || ix_next > domain.nx || iy_next < 1 || iy_next > domain.ny
                    continue
                end
                
                if (domain.cell_type[ix_next, iy_next] & CellType.INTERFACE) != 0
                    eta[n] = v[n][1] * normal[1] + v[n][2] * normal[2]
                    
                    if eta[n] < 0
                        eta[n] = 0
                    end
                    
                    eta_total += eta[n]
                    isIF[n] = 1
                    IF_total += 1
                end
            end
            
            if eta_total > 0
                eta_frac = 1 / eta_total
                for n in 2:9
                    domain.fdist[n, ix, iy] = mex * eta[n] * eta_frac
                end
            elseif IF_total > 0
                mex_rel = mex / IF_total
                for n in 2:9
                    domain.fdist[n, ix, iy] = isIF[n] == 1 ? mex_rel : 0.0
                end
            end
        end
    end
    
    # Fifth pass: collect distributed mass and finalize cell flags
    for ix in 1:domain.nx
        for iy in 1:domain.ny
            if (domain.cell_type[ix, iy] & CellType.INTERFACE) != 0
                for n in 2:9
                    ix_n = ix + v[n][1]
                    iy_n = iy + v[n][2]
                    
                    if 1 <= ix_n <= domain.nx && 1 <= iy_n <= domain.ny
                        domain.mass[ix, iy] += domain.fdist[n, ix_n, iy_n]
                    end
                end
            elseif (domain.cell_type[ix, iy] & CellType.TO_FLUID) != 0
                domain.cell_type[ix, iy] = CellType.FLUID
            elseif (domain.cell_type[ix, iy] & CellType.TO_GAS) != 0
                domain.cell_type[ix, iy] = CellType.GAS
            end
        end
    end
    
    # Sixth pass: set neighborhood flags
    for ix in 1:domain.nx
        for iy in 1:domain.ny
            if (domain.cell_type[ix, iy] & CellType.OBSTACLE) != 0
                continue
            end
            
            domain.cell_type[ix, iy] |= (CellType.NO_FLUID_NEIGH | CellType.NO_EMPTY_NEIGH | CellType.NO_IFACE_NEIGH)
            
            for n in 2:9
                ix_next = ix - v[n][1]
                iy_next = iy - v[n][2]
                
                if ix_next < 1 || ix_next > domain.nx || iy_next < 1 || iy_next > domain.ny
                    continue
                end
                
                if (domain.cell_type[ix_next, iy_next] & CellType.FLUID) != 0
                    domain.cell_type[ix, iy] &= ~CellType.NO_FLUID_NEIGH
                elseif (domain.cell_type[ix_next, iy_next] & CellType.GAS) != 0
                    domain.cell_type[ix, iy] &= ~CellType.NO_EMPTY_NEIGH
                elseif (domain.cell_type[ix_next, iy_next] & CellType.INTERFACE) != 0
                    domain.cell_type[ix, iy] &= ~CellType.NO_IFACE_NEIGH
                end
            end
        end
    end
end

# Main time evolution function
function evolve!(domain::LBMDomain, total_timesteps::Int)
    for time in 1:total_timesteps
        println("Timestep: $time")
        
        # Collision step
        lbm_collision!(domain)
        
        # Streaming step
        lbm_streaming!(domain)
        
        # Update cell types
        update_types!(domain)
        
        # Store history for visualization
        push!(domain.velocity_history, copy(domain.u))
        push!(domain.cell_type_history, copy(domain.cell_type))
    end
end

# Initialize a dam break simulation
function init_dam_break(nx::Int, ny::Int)
    domain = LBMDomain(nx=nx, ny=ny)
    
    # Create dam: fluid square
    minx, maxx = 75,100
    miny, maxy = 40, ny
    
    # Initialize all cells as gas
    domain.cell_type .= CellType.GAS
    
    # Create dam structure
    for ix in minx:maxx
        for iy in miny:maxy
            if ix > minx && ix < maxx && iy > miny && iy < maxy
                domain.cell_type[ix, iy] = CellType.FLUID
            else
                domain.cell_type[ix, iy] = CellType.INTERFACE
            end
        end
    end
    
    # Set boundary walls
    domain.cell_type[:, 1] .= CellType.OBSTACLE
    domain.cell_type[:, ny] .= CellType.OBSTACLE
    domain.cell_type[1, :] .= CellType.OBSTACLE
    domain.cell_type[nx, :] .= CellType.OBSTACLE
    
    # Set mass and density
    domain.mass[domain.cell_type .== CellType.FLUID] .= 1.0
    domain.mass[domain.cell_type .== CellType.INTERFACE] .= 0.5
    domain.rho .+= domain.mass
    
    # Initialize fin with equilibrium
    for ix in 1:nx
        for iy in 1:ny
            set_equi!(domain, ix, iy, domain.rho[ix, iy], domain.u[:, ix, iy])
        end
    end
    domain.fin .= domain.equi
    
    return domain
end

# Visualization helper
function visualize_dam_break(domain::LBMDomain, frame::Int)
    # Create a mask of fluid and interface cells
    fluid_mask = (domain.cell_type_history[frame] .& (CellType.FLUID | CellType.INTERFACE)) .!= 0
    
    # Plot the dam break state
    heatmap(1:domain.nx, 1:domain.ny, transpose(fluid_mask), 
            color=:blues, aspect_ratio=:equal, legend=false,
            title="Dam Break Simulation (Frame $frame)")
end

# Create animation from simulation results
# Modified create_animation function to accept StepRange or Vector
function create_animation(domain::LBMDomain, frames::Union{Vector{Int}, StepRange{Int, Int}}, filename::String="dam_break.gif")
    anim = @animate for frame in frames
        visualize_dam_break(domain, frame)
    end
    
    gif(anim, filename, fps=10)
    println("Animation saved as $filename")
end

function run_dam_break(nx::Int=170, ny::Int=170, timesteps::Int=750)
    println("Initializing dam break simulation...")
    domain = init_dam_break(nx, ny)
    
    println("Running simulation for $timesteps timesteps...")
    evolve!(domain, timesteps)
    
    # Create animation frames (every 10 timesteps)
    frames = 1:10:length(domain.cell_type_history)
    
    println("Creating animation...")
    create_animation(domain, frames)
    
    return domain
end

# Function to analyze dam break results
function analyze_dam_break(domain::LBMDomain)
    # Calculate maximum velocity at each timestep
    max_velocities = Vector{Float64}(undef, length(domain.velocity_history))
    for i in 1:length(domain.velocity_history)
        velocities = sqrt.(domain.velocity_history[i][1, :, :].^2 .+ domain.velocity_history[i][2, :, :].^2)
        max_velocities[i] = maximum(velocities)
    end
    
    # Plot maximum velocity over time
    p1 = plot(1:length(max_velocities), max_velocities, 
         xlabel="Time Step", ylabel="Maximum Velocity",
         title="Maximum Fluid Velocity", legend=false, lw=2)
    
    # Calculate fluid coverage (percentage of domain filled with fluid)
    fluid_coverage = Vector{Float64}(undef, length(domain.cell_type_history))
    for i in 1:length(domain.cell_type_history)
        fluid_cells = count((domain.cell_type_history[i] .& (CellType.FLUID | CellType.INTERFACE)) .!= 0)
        total_cells = domain.nx * domain.ny - count(domain.cell_type_history[i] .== CellType.OBSTACLE)
        fluid_coverage[i] = fluid_cells / total_cells * 100
    end
    
    # Plot fluid coverage over time
    p2 = plot(1:length(fluid_coverage), fluid_coverage,
         xlabel="Time Step", ylabel="Fluid Coverage (%)",
         title="Fluid Coverage", legend=false, lw=2)
    
    # Combine plots
    combined_plot = plot(p1, p2, layout=(2,1), size=(800, 600))
    savefig(combined_plot, "dam_break_analysis.png")
    
    return max_velocities, fluid_coverage
end

# Function to compare different gravity settings
function compare_gravity_settings(nx::Int=170, ny::Int=170, timesteps::Int=500)
    gravity_values = [
        [0.0, 0.05],
        [0.0, 0.1],
        [0.0, 0.2]
    ]
    
    labels = ["Low gravity", "Medium gravity", "High gravity"]
    results = []
    
    for (i, gravity) in enumerate(gravity_values)
        println("Running simulation with gravity = $gravity")
        domain = init_dam_break(nx, ny)
        domain.gravity = gravity
        
        evolve!(domain, timesteps)
        push!(results, domain)
        
        # Save a snapshot of the final state
        visualize_dam_break(domain, length(domain.cell_type_history))
        savefig("gravity_comparison_$(labels[i]).png")
    end
    
    # Compare final states
    p = plot(layout=(1,3), size=(1200, 400))
    for i in 1:3
        fluid_mask = (results[i].cell_type_history[end] .& (CellType.FLUID | CellType.INTERFACE)) .!= 0
        heatmap!(p[i], 1:nx, 1:ny, transpose(fluid_mask), 
                color=:blues, aspect_ratio=:equal, legend=false,
                title=labels[i])
    end
    
    savefig(p, "gravity_comparison.png")
    return results
end

# Function to run a custom obstacle simulation
function run_obstacle_simulation(nx::Int=170, ny::Int=170, timesteps::Int=750)
    domain = init_dam_break(nx, ny)
    
    # Add some obstacles in the flow path
    # Central pillar
    center_x, center_y = nx ÷ 2, ny ÷ 2
    radius = 10
    
    for ix in (center_x-radius):(center_x+radius)
        for iy in (center_y-radius):(center_y+radius)
            if (ix - center_x)^2 + (iy - center_y)^2 <= radius^2
                domain.cell_type[ix, iy] = CellType.OBSTACLE
            end
        end
    end
    
    # Add a barrier
    barrier_x = nx ÷ 2 + 30
    barrier_height = ny ÷ 3
    domain.cell_type[barrier_x, (ny÷2-barrier_height):(ny÷2+barrier_height)] .= CellType.OBSTACLE
    
    println("Running obstacle simulation for $timesteps timesteps...")
    evolve!(domain, timesteps)
    
    # Create animation
    frames = 1:10:length(domain.cell_type_history)
    create_animation(domain, frames, "obstacle_simulation.gif")
    
    return domain
end

# Function to visualize velocity field
function visualize_velocity_field(domain::LBMDomain, frame::Int, scale::Float64=5.0)
    # Extract fluid mask and velocity data for the given frame
    fluid_mask = (domain.cell_type_history[frame] .& (CellType.FLUID | CellType.INTERFACE)) .!= 0
    vx = domain.velocity_history[frame][1, :, :]
    vy = domain.velocity_history[frame][2, :, :]
    
    # Calculate velocity magnitude
    vmag = sqrt.(vx.^2 .+ vy.^2)
    
    # Sample points for quiver plot (too many arrows makes it unreadable)
    skip = 5  # Display every 5th point
    x = 1:skip:domain.nx
    y = 1:skip:domain.ny
    
    # Create the plot
    p = heatmap(1:domain.nx, 1:domain.ny, transpose(fluid_mask .* vmag),
                color=:thermal, aspect_ratio=:equal,
                title="Velocity Field (Frame $frame)")
    
    # Add velocity vectors
    quiver!(p, repeat(x, length(y)), [y[j] for i in 1:length(x), j in 1:length(y)][:],
            quiver=(scale .* vx[1:skip:end, 1:skip:end][:], scale .* vy[1:skip:end, 1:skip:end][:]),
            color=:white, alpha=0.7)
    
    return p
end

# Main function to demonstrate all capabilities
function demo()
    println("Running FreeSurfaceLBM Dam Break Demo")
    println("-------------------------------------")
    
    # Run basic dam break
    domain = run_dam_break(170, 170, 1000)
    
    # Analyze results
    max_vel, coverage = analyze_dam_break(domain)
    
    # Create velocity field visualization for selected frames
    for frame in [50, 100, 200, 300, 400]
        p = visualize_velocity_field(domain, frame)
        savefig(p, "velocity_field_frame_$frame.png")
    end
    
    println("Demo completed successfully!")
    println("Output files created:")
    println("  - dam_break.gif: Animation of the dam break")
    println("  - dam_break_analysis.png: Plots of velocity and coverage")
    println("  - velocity_field_frame_*.png: Velocity field visualizations")
    
    return domain
end

# Execute demo if this file is run directly
demo()


end # module FreeSurfaceLBM
