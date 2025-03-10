using Plots

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
    R_z = z * (1 + (zR/z)^2)
    R_z = z == 0 ? Inf : R_z  # Handle z=0 case
    
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

# Define simulation parameters
struct SimParams
    w0::Float64      # Beam waist (m)
    λ::Float64       # Wavelength (m)
    t0::Float64      # Pulse center time (s)
    τ::Float64       # Pulse duration (s)
end

# Setting up parameters based on the provided equation
params = SimParams(
    2.0e-6,          # Beam waist: 2 μm
    800e-9,          # Wavelength: 800 nm
    100e-15,         # Pulse center time: 100 fs
    30e-15           # Pulse duration: 30 fs
)

# Function to calculate intensity
function calculate_intensity(E)
    # I = (n/2)√(ε₀/μ₀)|E|² where n is refractive index
    # For simplicity, assuming n=1 and normalizing constants
    return abs2.(E)
end

# Create a visualization of the pulse
function visualize_pulse()
    # Space and time grid
    x_range = range(-10e-6, 10e-6, length=100)  # -10 to 10 μm
    t_range = range(0, 200e-15, length=100)     # 0 to 200 fs
    
    # Calculate field at z=0 (focus)
    z = 0.0
    
    # Calculate E-field
    E = [gaussian_beam_electric_field(t, x, z, params) for x in x_range, t in t_range]
    
    # Calculate intensity
    I = calculate_intensity(E)
    
    # Normalize for visualization
    I_norm = I ./ maximum(I)
    
    # Create heatmap
    heatmap(t_range*1e15, x_range*1e6, I_norm, 
            xlabel="Time (fs)", ylabel="x (μm)", 
            title="Intensity of Gaussian Pulse at z=0",
            color=:thermal)
end

# Function to simulate pulse propagation along z
function simulate_propagation()
    # Space and z grid
    x_range = range(-10e-6, 10e-6, length=100)   # -10 to 10 μm
    z_range = range(-40e-6, 40e-6, length=100)   # -40 to 40 μm
    
    # Fixed time (at pulse center)
    t = params.t0
    
    # Calculate E-field
    E = [gaussian_beam_electric_field(t, x, z, params) for x in x_range, z in z_range]
    
    # Calculate intensity
    I = calculate_intensity(E)
    
    # Normalize for visualization
    I_norm = I ./ maximum(I)
    
    # Create heatmap
    heatmap(z_range*1e6, x_range*1e6, I_norm, 
            xlabel="z (μm)", ylabel="x (μm)", 
            title="Beam Intensity at t=$(params.t0*1e15) fs",
            color=:thermal)
end

# To calculate absorbed energy
function calculate_absorbed_energy(E, α_abs)
    # Absorbed energy = I * α_abs, where I = |E|²
    # α_abs is the bremsstrahlung absorption coefficient
    
    I = calculate_intensity(E)
    absorbed = I .* α_abs
    
    return absorbed
end

# Example: Calculate and visualize pulse at different propagation distances
function visualize_pulse_evolution()
    t = params.t0  # At peak time
    x_range = range(-10e-6, 10e-6, length=100)
    
    z_positions = [0.0, params.w0^2*π/params.λ, 2*params.w0^2*π/params.λ]  # At focus, at Rayleigh length, at 2*Rayleigh
    
    plots = []
    for z in z_positions
        E = [gaussian_beam_electric_field(t, x, z, params) for x in x_range]
        I = calculate_intensity(E)
        I_norm = I ./ maximum(I)
        
        p = plot(x_range*1e6, I_norm, 
              xlabel="x (μm)", ylabel="Normalized Intensity", 
              title="z = $(round(z*1e6, digits=1)) μm",
              legend=false, lw=2)
        push!(plots, p)
    end
    
    plot(plots..., layout=(3,1), size=(800, 600))
end

# Run visualizations
visualize_pulse()
simulate_propagation()
visualize_pulse_evolution()

