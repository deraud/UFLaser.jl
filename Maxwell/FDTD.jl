#= 
2D FDTD code with unidirectional soft source and absorbing PML boundary conditions. 
This code can be used to simulate optical phenomena.

Users can adjust the source wavelength, as well as the object geometry
and material properties. 
        
          -------------------------------------------------------
          |                    BACK PML                         |
          -------------------------------------------------------
          |          |                              |           |
          |          |                              |           |
          | Left PML |         MAIN GRID            | Right PML |
          |          |                              |           |
          |          |                              |           |
          -------------------------------------------------------
          |                    FRONT PML                        |
          -------------------------------------------------------   

Original property of Ana Ciocoiu, University of British Columbia
Converted to Julia
=#

using Plots
using DSP  # For hilbert transform
using Statistics

# Clear environment and close plots (Julia equivalent of clear all, close all)
# In Julia, you don't need to clear variables as in MATLAB

## 1: Fundamental constants
const epsilon_0 = 8.85e-12           # permittivity of free space    
const mu_0 = 4*π*1e-7                # permeability of free space
const c = 1/sqrt(epsilon_0*mu_0)     # speed of light (m/s)
const eta_0 = sqrt(mu_0/epsilon_0)   # impedance of free space (Ohms) 

println("Please select simulator mode: 1 is TE, 0 is TM: >>")
mode = 1        # select simulator mode (1 is TE, 0 is TM)

## 2: Units
const micrometers = 1e-6             # the units used will affect the grid spacing as well as the PML boundary conditions. 

## 3: Grid parameters
# This section defines the spacing of the 2-D grid, and sets the time step
# increments to satisfy stability (Courant) conditions.
#
delta_x = 6*micrometers              # spatial step increment
delta_z = delta_x                    
delta_y = delta_z                    
delta_t = delta_z/(2*c)              # Courant condition for stability of simulator
Nz = 500                             # spatial grid size
Nx = 500
Nt = 3000                            # number of iterations in time. Nt*delta_t = total duration (seconds)

x = (1:Nx) .* delta_x                # define space and time step arrays
y = (1:Nx) .* delta_y
z = (1:Nz) .* delta_z
t = (1:Nt) .* delta_t

# Set up 2-D grid for TE or TM modes
if mode == 0
    Z, X = meshgrid(z, x)
elseif mode == 1
    Y, X = meshgrid(x, y)
end

# Initialize a set of test points and arrays
Ey_test = zeros(Nz, Nt)
Hz_test = zeros(Nz, Nt)
E_array = zeros(Nx, Nz)
E_array_new = zeros(Nx, Nz)

## 4: Source parameters. 
# This section sets up the properties of the electromagnetic source
source_z = 60                        # location of source as a function of location within the grid points (Nz)
wavelength = 300*micrometers         # source central wavelength
omega = 2 * π * c/wavelength         # source central frequency
beam_waist = 500*micrometers         # width of gaussian beam
beam_centre = Nx/2*delta_x           # location of beam center

## 5: Field initialization.
# This section sets up split-field electric and magnetic fields for a the PML boundary (absorbing layer).

# Initialize fields for both TE and TM modes
# TM with split-field for PML
Hx = zeros(Nx, Nz)
Hz = zeros(Nx, Nz)
Eyz = zeros(Nx, Nz)
Eyx = zeros(Nx, Nz)

# TE with split-field for PML
Ex = zeros(Nx, Nz)
Ey = zeros(Nx, Nz)
Hzx = zeros(Nx, Nz)
Hzy = zeros(Nx, Nz)

# Material properties initialization
epsilon_r = ones(Nx, Nz)    # relative permittivity set to default of 1 across the grid (air)
mu_r = ones(Nx, Nz)         # relative permeability set to default of 1 across the grid

# Initialize conductivity in both x and z directions (S/m) within the grid
sigma_x = zeros(Nx, Nz)
sigma_mx = zeros(Nx, Nz)
sigma_z = zeros(Nx, Nz)
sigma_mz = zeros(Nx, Nz)

## PML region
PML_width = 40               # thickness of PML region around the grid
PML_cond = 80                # condition for stable and functional PML

# This section sets PML conditions for total wave attenuation within the region
sigma_x[1:PML_width, :] .= PML_cond
sigma_x[Nx-PML_width+1:Nx, :] .= PML_cond
sigma_mx[1:PML_width, :] .= PML_cond*mu_0/epsilon_0
sigma_mx[Nx-PML_width+1:Nx, :] .= PML_cond*mu_0/epsilon_0

sigma_z[:, 1:PML_width] .= PML_cond
sigma_z[:, Nz-PML_width+1:Nz] .= PML_cond
sigma_mz[:, 1:PML_width] .= PML_cond*mu_0/epsilon_0
sigma_mz[:, Nz-PML_width+1:Nz] .= PML_cond*mu_0/epsilon_0

## 6: Object definition.
# This section sets up the geometry of the simulation.

# Currently set up to simulate a 300um diameter sphere of permittivity 2.25
front_sphere_radius = (25*delta_x)

# Define permittivity of object based on the mode
if mode == 0
    # Create a meshgrid for calculating distances
    for i in 1:Nx, j in 1:Nz
        if ((X[i,j]-Nx/2*delta_x)^2 + (Z[i,j]-Nz/2*delta_z)^2) <= front_sphere_radius^2
            epsilon_r[i,j] = 2.25
        end
    end
    
    # Display the permittivity profile
    heatmap(z.*1e3, x.*1e3, epsilon_r, aspect_ratio=:equal, xlabel="z [μm]", ylabel="x [μm]")
elseif mode == 1
    # Create a meshgrid for calculating distances
    for i in 1:Nx, j in 1:Nz
        if ((X[i,j]-Nx/2*delta_x)^2 + (Y[i,j]-Nz/2*delta_z)^2) <= front_sphere_radius^2
            epsilon_r[i,j] = 2.25
        end
    end
    
    # Display the permittivity profile
    heatmap(y.*1e3, x.*1e3, epsilon_r, aspect_ratio=:equal, xlabel="y [μm]", ylabel="x [μm]")
end

## 7: FDTD factors.
# Setting up multipliers for fields of the FDTD kernel

Hxfactor1 = 1 .- sigma_mx .* delta_t./(2 .* mu_r .* mu_0)
Hxfactor2 = 1 .+ sigma_mx .* delta_t./(2 .* mu_r .* mu_0)
Hzfactor1 = 1 .- sigma_mz .* delta_t./(2 .* mu_r .* mu_0)
Hzfactor2 = 1 .+ sigma_mz .* delta_t./(2 .* mu_r .* mu_0)

Exfactor1 = 1 .- sigma_x .* delta_t./(2 .* epsilon_r .* epsilon_0)
Exfactor2 = 1 .+ sigma_x .* delta_t./(2 .* epsilon_r .* epsilon_0)
Ezfactor1 = 1 .- sigma_z .* delta_t./(2 .* epsilon_r .* epsilon_0)
Ezfactor2 = 1 .+ sigma_z .* delta_t./(2 .* epsilon_r .* epsilon_0)

reduction = 10   # for visualization purposes

## 8: FDTD kernel for TE of TM modes
# Create empty arrays for test points
if mode == 0
    Hx_test = zeros(Nz, Nt)
    Hx_test2 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Hx_test3 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Hx_test4 = zeros(Nz, Nt)
    
    Ey_test2 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Ey_test3 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Ey_test4 = zeros(Nz, Nt)
    
    Hz_test2 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Hz_test3 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Hz_test4 = zeros(Nz, Nt)
else
    Ex_test = zeros(Nz, Nt)
    Ex_test2 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Ex_test3 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Ex_test4 = zeros(Nz, Nt)
    
    Ey_test2 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Ey_test3 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Ey_test4 = zeros(Nz, Nt)
    
    Hz_test2 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Hz_test3 = zeros(Nt, Int(0.9*Nz - 1.5*source_z) + 1)
    Hz_test4 = zeros(Nz, Nt)
end

# Run the FDTD simulation
if mode == 0
    for n in 1:Nt
        # Update E field components
        Eyz[2:Nx, 2:Nz] = Eyz[2:Nx, 2:Nz] .* Ezfactor1[2:Nx, 2:Nz] ./ Ezfactor2[2:Nx, 2:Nz] .+ 
            delta_t ./ (epsilon_0 .* epsilon_r[2:Nx, 2:Nz] .* Ezfactor2[2:Nx, 2:Nz]) .* 
            ((Hx[2:Nx, 2:Nz] - Hx[2:Nx, 1:Nz-1]) ./ delta_z)
        
        Eyx[2:Nx, 2:Nz] = Eyx[2:Nx, 2:Nz] .* Exfactor1[2:Nx, 2:Nz] ./ Exfactor2[2:Nx, 2:Nz] .+ 
            delta_t ./ (epsilon_0 .* epsilon_r[2:Nx, 2:Nz] .* Exfactor2[2:Nx, 2:Nz]) .* 
            (-(Hz[2:Nx, 2:Nz] - Hz[1:Nx-1, 2:Nz]) ./ delta_x)
        
        Ey[2:Nx, 2:Nz] = Eyz[2:Nx, 2:Nz] .+ Eyx[2:Nx, 2:Nz]
        
        # Update H field components
        Hx[1:Nx, 1:Nz-1] = Hx[1:Nx, 1:Nz-1] .* Hzfactor1[1:Nx, 1:Nz-1] ./ Hzfactor2[1:Nx, 1:Nz-1] .+ 
            delta_t ./ (mu_0 .* mu_r[1:Nx, 1:Nz-1] .* Hzfactor2[1:Nx, 1:Nz-1]) .* 
            (Ey[1:Nx, 2:Nz] - Ey[1:Nx, 1:Nz-1]) ./ delta_z
        
        Hz[1:Nx-1, 1:Nz] = Hz[1:Nx-1, 1:Nz] .* Hxfactor1[1:Nx-1, 1:Nz] ./ Hxfactor2[1:Nx-1, 1:Nz] .- 
            delta_t ./ (mu_0 .* mu_r[1:Nx-1, 1:Nz] .* Hxfactor2[1:Nx-1, 1:Nz]) .* 
            (Ey[2:Nx, 1:Nz] - Ey[1:Nx-1, 1:Nz]) ./ delta_x
        
        # Initialize source in x direction
        Hz[source_z, :] = Hz[source_z, :] .+ sin(omega*t[n]) .* exp.(- ((x .- beam_centre) ./ beam_waist).^2) ./ eta_0
        Eyx[source_z, :] = Eyx[source_z, :] .+ sin(omega*t[n]) .* exp.(- ((x .- beam_centre) ./ beam_waist).^2)
        
        # Test points for Poynting vector calculations
        Ey_test[:, n] = Ey[1:Nz, Int(1.5*source_z)]
        Hx_test[:, n] = Hx[1:Nz, Int(1.5*source_z)]
        Hz_test[:, n] = Hz[1:Nz, Int(1.5*source_z)]
        
        range_idx = Int(1.5*source_z):Int(0.9*Nz)
        Ey_test2[n, 1:length(range_idx)] = Ey[2*PML_width, range_idx]
        Hx_test2[n, 1:length(range_idx)] = Hx[2*PML_width, range_idx]
        Hz_test2[n, 1:length(range_idx)] = Hz[2*PML_width, range_idx]
        
        Ey_test3[n, 1:length(range_idx)] = Ey[Int(0.9*Nx), range_idx]
        Hx_test3[n, 1:length(range_idx)] = Hx[Int(0.9*Nx), range_idx]
        Hz_test3[n, 1:length(range_idx)] = Hz[Int(0.9*Nx), range_idx]
        
        Ey_test4[:, n] = Ey[1:Nz, Int(0.9*Nz)]
        Hx_test4[:, n] = Hx[1:Nz, Int(0.9*Nz)]
        Hz_test4[:, n] = Hz[1:Nz, Int(0.9*Nz)]
        
        # Uncomment to create real-time plot of the simulation
        # if n % reduction == 0
        #     heatmap(z, x, Ey, aspect_ratio=:equal, xlabel="z [m]", ylabel="x [m]")
        # end
        
        # Store maximum field values for visualization
        if n > 0.8*Nt
            E_array_new = copy(Ey)
            comp_arr_GT = E_array_new .> E_array
            comp_arr_LT = E_array_new .<= E_array
            E_array_new = E_array_new .* comp_arr_GT .+ E_array .* comp_arr_LT
            E_array = copy(E_array_new)
        end
    end
else # mode == 1
    for n in 1:Nt
        # Update Hz field components
        Hzx[1:Nx-1, 1:Nz-1] = Hzx[1:Nx-1, 1:Nz-1] .* Hxfactor1[1:Nx-1, 1:Nz-1] ./ Hxfactor2[1:Nx-1, 1:Nz-1] .- 
            delta_t ./ (mu_0 .* mu_r[2:Nx, 2:Nz] .* Hxfactor2[1:Nx-1, 1:Nz-1]) .* 
            (Ey[2:Nx, 1:Nz-1] - Ey[1:Nx-1, 1:Nz-1]) ./ delta_x
        
        Hzy[1:Nx-1, 1:Nz-1] = Hzy[1:Nx-1, 1:Nz-1] .* Hzfactor1[1:Nx-1, 1:Nz-1] ./ Hzfactor2[1:Nx-1, 1:Nz-1] .- 
            delta_t ./ (mu_0 .* mu_r[2:Nx, 2:Nz] .* Hzfactor2[1:Nx-1, 1:Nz-1]) .* 
            (-(Ex[1:Nx-1, 2:Nz] - Ex[1:Nx-1, 1:Nz-1]) ./ delta_y)
        
        Hz[1:Nx-1, 1:Nz-1] = Hzy[1:Nx-1, 1:Nz-1] .+ Hzx[1:Nx-1, 1:Nz-1]
        
        # Update E field components
        Ey[2:Nx, 2:Nz] = Ey[2:Nx, 2:Nz] .* Exfactor1[2:Nx, 2:Nz] ./ Exfactor2[2:Nx, 2:Nz] .- 
            delta_t ./ (epsilon_0 .* epsilon_r[2:Nx, 2:Nz] .* Exfactor2[2:Nx, 2:Nz]) .* 
            (Hz[2:Nx, 2:Nz] - Hz[1:Nx-1, 2:Nz]) ./ delta_x
        
        Ex[2:Nx, 2:Nz] = Ex[2:Nx, 2:Nz] .* Ezfactor1[2:Nx, 2:Nz] ./ Ezfactor2[2:Nx, 2:Nz] .+ 
            delta_t ./ (epsilon_0 .* epsilon_r[2:Nx, 2:Nz] .* Ezfactor2[2:Nx, 2:Nz]) .* 
            (Hz[2:Nx, 2:Nz] - Hz[2:Nx, 1:Nz-1]) ./ delta_y
        
        # Initialize source in x direction
        Hzx[source_z, :] = Hzx[source_z, :] .+ sin(omega*t[n]) .* exp.(- ((x .- beam_centre) ./ beam_waist).^2) ./ eta_0
        Ey[source_z, :] = Ey[source_z, :] .+ sin(omega*t[n]) .* exp.(- ((x .- beam_centre) ./ beam_waist).^2)
        
        # Test points for Poynting vector calculations
        Ey_test[:, n] = Ey[1:Nz, Int(1.5*source_z)]
        Ex_test[:, n] = Ex[1:Nz, Int(1.5*source_z)]
        Hz_test[:, n] = Hz[1:Nz, Int(1.5*source_z)]
        
        range_idx = Int(1.5*source_z):Int(0.9*Nz)
        Ey_test2[n, 1:length(range_idx)] = Ey[2*PML_width, range_idx]
        Ex_test2[n, 1:length(range_idx)] = Ex[2*PML_width, range_idx]
        Hz_test2[n, 1:length(range_idx)] = Hz[2*PML_width, range_idx]
        
        Ey_test3[n, 1:length(range_idx)] = Ey[Int(0.9*Nx), range_idx]
        Ex_test3[n, 1:length(range_idx)] = Ex[Int(0.9*Nx), range_idx]
        Hz_test3[n, 1:length(range_idx)] = Hz[Int(0.9*Nx), range_idx]
        
        Ey_test4[:, n] = Ey[1:Nz, Int(0.9*Nz)]
        Ex_test4[:, n] = Ex[1:Nz, Int(0.9*Nz)]
        Hz_test4[:, n] = Hz[1:Nz, Int(0.9*Nz)]
        
        # Uncomment to create real-time plot of the simulation
        # if n % reduction == 0
        #     heatmap(y, x, Hz, aspect_ratio=:equal, xlabel="y [m]", ylabel="x [m]")
        # end
        
        # Store maximum field values for visualization
        if n > 0.8*Nt
            E_array_new = copy(Hz)
            comp_arr_GT = E_array_new .> E_array
            comp_arr_LT = E_array_new .<= E_array
            E_array_new = E_array_new .* comp_arr_GT .+ E_array .* comp_arr_LT
            E_array = copy(E_array_new)
        end
    end
end

## 9: Poynting vector calculations for energy conservation test
threshold = 1e-5  # error threshold for conservation of energy calculations

# Helper function to create meshgrid similar to MATLAB's meshgrid
function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))
    return Y, X
end

# Calculate Poynting vectors for energy conservation test
if mode == 1
    # For TE mode
    H_poynt = hcat(zeros(Nz), zeros(Nz), Hz_test[:, end])
    E_poynt = hcat(Ex_test[:, end], Ey_test[:, end], zeros(Nz))
    S = 0.5 .* real.(cross.(hilbert.(E_poynt), conj.(hilbert.(H_poynt))))
    
    # Select smaller range for test points 2 and 3
    test_range = 1:201
    
    H_poynt2 = hcat(zeros(length(test_range)), zeros(length(test_range)), Hz_test2[end, test_range])
    E_poynt2 = hcat(Ex_test2[end, test_range], Ey_test2[end, test_range], zeros(length(test_range)))
    S2 = 0.5 .* real.(cross.(hilbert.(E_poynt2), conj.(hilbert.(H_poynt2))))
    
    H_poynt3 = hcat(zeros(length(test_range)), zeros(length(test_range)), Hz_test3[end, test_range])
    E_poynt3 = hcat(Ex_test3[end, test_range], Ey_test3[end, test_range], zeros(length(test_range)))
    S3 = 0.5 .* real.(cross.(hilbert.(E_poynt3), conj.(hilbert.(H_poynt3))))
    
    H_poynt4 = hcat(zeros(Nz), zeros(Nz), Hz_test4[:, end])
    E_poynt4 = hcat(Ex_test4[:, end], Ey_test4[:, end], zeros(Nz))
    S4 = 0.5 .* real.(cross.(hilbert.(E_poynt4), conj.(hilbert.(H_poynt4))))
    
else
    # For TM mode
    H_poynt = hcat(Hx_test[:, end], zeros(Nz), Hz_test[:, end])
    E_poynt = hcat(zeros(Nz), Ey_test[:, end], zeros(Nz))
    S = 0.5 .* real.(cross.(hilbert.(E_poynt), conj.(hilbert.(H_poynt))))
    
    # Select smaller range for test points 2 and 3
    test_range = 1:201
    
    H_poynt2 = hcat(Hx_test2[end, test_range], zeros(length(test_range)), Hz_test2[end, test_range])
    E_poynt2 = hcat(zeros(length(test_range)), Ey_test2[end, test_range], zeros(length(test_range)))
    S2 = 0.5 .* real.(cross.(hilbert.(E_poynt2), conj.(hilbert.(H_poynt2))))
    
    H_poynt3 = hcat(Hx_test3[end, test_range], zeros(length(test_range)), Hz_test3[end, test_range])
    E_poynt3 = hcat(zeros(length(test_range)), Ey_test3[end, test_range], zeros(length(test_range)))
    S3 = 0.5 .* real.(cross.(hilbert.(E_poynt3), conj.(hilbert.(H_poynt3))))
    
    H_poynt4 = hcat(Hx_test4[:, end], zeros(Nz), Hz_test4[:, end])
    E_poynt4 = hcat(zeros(Nz), Ey_test4[:, end], zeros(Nz))
    S4 = 0.5 .* real.(cross.(hilbert.(E_poynt4), conj.(hilbert.(H_poynt4))))
end

# Check conservation of energy
sumtotal = -sum(getindex.(S, 3)) + sum(getindex.(S4, 3)) + sum(getindex.(S3, 1)) - sum(getindex.(S2, 1))
if abs(sumtotal) < threshold
    println("stable: hooray!")
else
    println("not stable: try reducing grid spacing!")
end

## 10: Plot the maximum E or H-field array obtained as a heat map
if mode == 0
    heatmap(z./delta_x, x./delta_x, E_array_new, 
            aspect_ratio=:equal, 
            xlabel="z/dx", 
            ylabel="x/dx", 
            title="Maximum E-field Distribution (TM mode)")
else
    heatmap(y./delta_x, x./delta_x, E_array_new, 
            aspect_ratio=:equal, 
            xlabel="y/dx", 
            ylabel="x/dx", 
            title="Maximum H-field Distribution (TE mode)")
end