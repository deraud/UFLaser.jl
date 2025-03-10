using StaticArrays, Einsum, Plots

using Profile

include("././struct.jl")
include("././constants.jl")
include("././laser_source.jl")
include("././material_thermal.jl")
include("././coupling.jl")
include("././LBM_Electron.jl")
include("././LBM_Lattice.jl")
include("././VOF.jl")
include("././VOF_TTM.jl")

using .Constants
using .Structs
using .LaserEnergy
using .Thermal
using .Coupling
using .lbElectron
using .lbLattice
using .FreeSurfaceLBM
using .TTM_VOF_Coupling
# Initialize Lattice
# Initialize TTM
ttm_lattice = Lattice()
thermal = ThermalProperties()

# Initialize VOF with dimensions matching TTM
vof_domain = LBMDomain(nx=Nx, ny=Ny)

# Initialize VOF from initial TTM state
vof_domain = initializeVOFFromTTM(ttm_lattice, vof_domain)

ElectronTemperature = zeros(Nx, Ny, Nt)
SourceLaserEnergy = zeros(Nx, Ny, Nt)

# Main simulation loop
for k in 1:Nt
    t = timeToFemto[k]
    
    # Run TTM step
    calculateLaserEnergy2DGaussianAblation!(t, ttm_lattice)
    calculateG!(ttm_lattice)
    
    ElectronTemperature[:,:,k] = ttm_lattice.Tₑ[:,:]
    SourceLaserEnergy[:,:,k] = ttm_lattice.S[:,:]

    # Electron step
    CollisionElectronLB!(ttm_lattice, thermal)
    streamingElectron!(ttm_lattice)
    updateElectron!(ttm_lattice)
    boundaryElectronAblation!(ttm_lattice)
    
    # Lattice step
    CollisionLatticeLB!(ttm_lattice, thermal)
    updateLattice!(ttm_lattice)
    
    # Calculate ablation in TTM
    #if k % 10 == 0
    #    calculateAblation!(ttm_lattice, LatticeTemperature, k)
    #end
    
    # Update VOF based on TTM results
    updateVOFFromTTM(ttm_lattice, vof_domain)
    
    # Run VOF step
    lbm_collision!(vof_domain)
    lbm_streaming!(vof_domain)
    update_types!(vof_domain)
    
    # Save VOF data for visualization
    if k % 10 == 0
        push!(vof_domain.velocity_history, copy(vof_domain.u))
        push!(vof_domain.cell_type_history, copy(vof_domain.cell_type))
    end
    
    # Print progress
    if k % 10 == 0
        println("Step $k completed out of $Nt | $(k/Nt*100)%")
    end
end








