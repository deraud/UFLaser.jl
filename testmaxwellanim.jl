using Plots
using FFTW
using LinearAlgebra

# Define simulation parameters
struct SimParams
    w0::Float64      # Beam waist (m)
    λ::Float64       # Wavelength (m)
    t0::Float64      # Pulse center time (s)
    τ::Float64       # Pulse duration (s)
    c::Float64       # Speed of light (m/s)
    n::Float64       # Refractive index
    k::Float64       # Extinction coefficient (for absorption)
end

# Setting up parameters based on the provided equation
params = SimParams(
    2.0e-6,          # Beam waist: 2 μm
    800e-9,          # Wavelength: 800 nm
    200e-15,         # Pulse center time: 200 fs
    30e-15,          # Pulse duration: 30 fs
    3.0e8,           # Speed of light
    1.5,             # Refractive index (typical for glass)
    0.01             # Extinction coefficient
)

function gaussian_beam_electric_field(t, x, z, params)
    # Unpack parameters
    w0 = params.w0           # Beam waist
    λ = params.λ             # Wavelength
    t0 = params.t0           # Pulse center time
    τ = params.τ             # Pulse duration
    
    # Derived parameters
    k = 2π / λ               # Wave number
    zR = π * w0^2 / λ        # Rayleigh length
    
    # Beam width at position z
    w_z = w0 * sqrt(1 + (z/zR)^2)
    
    # Radius of curvature
    R_z = z == 0 ? Inf : z * (1 + (zR/z)^2)
    
    # Gouy phase
    ζ_z = atan(z/zR)
    
    # Temporal part
    temporal = exp(-(t - t0)^2 / τ^2)
    
    # Spatial part - combining all terms from the equation
    spatial = exp(-x^2 / w_z^2 - 2im*π*z/λ - im*π*x^2/(R_z*λ) + im*ζ_z)
    
    # Amplitude factor
    amplitude = w0/w_z
    
    # Complete electric field
    return amplitude * temporal * spatial
end

# Calculate intensity from electric field
function calculate_intensity(E)
    # I = (n/2)√(ε₀/μ₀)|E|² where n is refractive index
    return abs2.(E)
end

# Calculate absorbed energy based on extinction coefficient
function calculate_absorbed_energy(I, params)
    # Absorption coefficient α = 4πk/λ
    α_abs = 4π * params.k / params.λ
    
    # Absorbed energy = I * α_abs
    return I .* α_abs
end

# Create animation of the pulse propagation
function animate_pulse_propagation()
    # Space grid
    nx, nz = 100, 100
    x_range = range(-10e-6, 10e-6, length=nx)  # -10 to 10 μm
    z_range = range(-20e-6, 20e-6, length=nz)  # -20 to 20 μm
    
    # Time range for animation
    t_start = 0
    t_end = 400e-15  # 400 fs
    nt = 60          # Number of frames
    t_range = range(t_start, t_end, length=nt)
    
    # Prepare animation
    anim = @animate for (i, t) in enumerate(t_range)
        # Calculate E-field at this time
        E = [gaussian_beam_electric_field(t, x, z, params) for x in x_range, z in z_range]
        
        # Calculate intensity
        I = calculate_intensity(E)
        
        # Calculate absorbed energy
        absorbed = calculate_absorbed_energy(I, params)
        
        # Normalize for visualization
        I_norm = I ./ max(maximum(I), 1e-10)
        absorbed_norm = absorbed ./ max(maximum(absorbed), 1e-10)
        
        # Create subplots
        p1 = heatmap(z_range*1e6, x_range*1e6, real.(E), 
                title="Real(Ex) at t=$(round(t*1e15, digits=1)) fs",
                xlabel="z (μm)", ylabel="x (μm)", 
                color=:bluesreds, clims=(-1,1))
                
        p2 = heatmap(z_range*1e6, x_range*1e6, absorbed_norm, 
                title="Absorbed Energy",
                xlabel="z (μm)", ylabel="x (μm)", 
                color=:thermal)
        
        plot(p1, p2, layout=(2,1), size=(800, 600))
    end
    
    # Save animation
    gif(anim, "laser_pulse_animation.gif", fps=15)
    return anim
end

# Function to simulate periodic structure interaction
function simulate_periodic_structure()
    # Define periodic structure
    # For example, a simple grating with alternating absorption properties
    nx, nz = 200, 200
    x_range = range(-10e-6, 10e-6, length=nx)  # -10 to 10 μm
    z_range = range(-20e-6, 20e-6, length=nz)  # -20 to 20 μm
    
    # Create periodic structure: alternating high/low absorption coefficient
    period = 2e-6  # 2 μm period
    duty_cycle = 0.5
    
    # Create absorption coefficient map
    k_high = 0.05
    k_low = 0.001
    k_map = zeros(nx, nz)
    
    for i in 1:nx, j in 1:nz
        x, z = x_range[i], z_range[j]
        # Create periodic pattern along x
        if mod(x, period) < period * duty_cycle
            k_map[i,j] = k_high
        else
            k_map[i,j] = k_low
        end
    end
    
    # Time range for animation
    t_start = 0
    t_end = 400e-15  # 400 fs
    nt = 60          # Number of frames
    t_range = range(t_start, t_end, length=nt)
    
    # Prepare animation
    anim = @animate for (i, t) in enumerate(t_range)
        # Calculate E-field at this time
        E = [gaussian_beam_electric_field(t, x, z, params) for x in x_range, z in z_range]
        
        # Calculate intensity
        I = calculate_intensity(E)
        
        # Calculate position-dependent absorbed energy
        α_map = 4π * k_map ./ params.λ
        absorbed = I .* α_map
        
        # Normalize for visualization
        I_norm = I ./ max(maximum(I), 1e-10)
        absorbed_norm = absorbed ./ max(maximum(absorbed), 1e-10)
        
        # Show the structure
        structure = heatmap(z_range*1e6, x_range*1e6, k_map, 
                   title="Periodic Structure (k)",
                   xlabel="z (μm)", ylabel="x (μm)", 
                   color=:grays)
        
        # Electric field
        field = heatmap(z_range*1e6, x_range*1e6, real.(E), 
                title="Real(Ex) at t=$(round(t*1e15, digits=1)) fs",
                xlabel="z (μm)", ylabel="x (μm)", 
                color=:bluesreds, clims=(-1,1))
        
        # Absorbed energy
        absorption = heatmap(z_range*1e6, x_range*1e6, absorbed_norm, 
                    title="Absorbed Energy",
                    xlabel="z (μm)", ylabel="x (μm)", 
                    color=:thermal)
        
        plot(structure, field, absorption, layout=(3,1), size=(800, 900))
    end
    
    # Save animation
    gif(anim, "periodic_structure_interaction.gif", fps=15)
    return anim
end

# Function to visualize Gouy phase shift
function visualize_gouy_phase()
    # Space grid
    z_range = range(-40e-6, 40e-6, length=200)  # -40 to 40 μm
    x = 0.0  # On-axis
    t = params.t0  # At peak time
    
    # Calculate E-field along z-axis
    E = [gaussian_beam_electric_field(t, x, z, params) for z in z_range]
    
    # Extract phase
    phase = angle.(E)
    
    # Calculate Gouy phase shift analytically
    zR = π * params.w0^2 / params.λ  # Rayleigh length
    gouy_phase = [atan(z/zR) for z in z_range]
    
    # Plot
    plot(z_range*1e6, phase, 
         label="Total Phase", xlabel="z (μm)", ylabel="Phase (rad)",
         title="Gouy Phase Shift", lw=2)
    plot!(z_range*1e6, gouy_phase, label="Gouy Phase Term", lw=2, ls=:dash)
    plot!(size=(800, 400))
end

# Run the animations
println("Generating pulse propagation animation...")
animate_pulse_propagation()

println("Generating periodic structure interaction animation...")
simulate_periodic_structure()

println("Visualizing Gouy phase shift...")
p_gouy = visualize_gouy_phase()
savefig(p_gouy, "gouy_phase_shift.png")

println("Done! Animations saved as 'laser_pulse_animation.gif' and 'periodic_structure_interaction.gif'")