using StaticArrays, Einsum, Plots
using Revise
using Profile

include("././struct.jl")
include("././constants.jl")
include("././laser_source.jl")
include("././material_thermal.jl")
include("././coupling.jl")
include("././LBM_Electron.jl")
include("././LBM_Lattice.jl")

using .Constants
using .Structs
using .LaserEnergy
using .Thermal
using .Coupling
using .lbElectron
using .lbLattice

# Initialize Lattice
@time lattice = Lattice()
@time thermal = ThermalProperties()


# Simulate energy distribution over time
#SourceLaserEnergy = zeros(Nx, Ny, Nt)  
ElectronTemperature = zeros(Nx, Ny, trunc(Int, Nt/10.0))
LatticeTemperature = zeros(Nx, Ny, trunc(Int, Nt/10.0))
AblationData = zeros(Nx,Ny,trunc(Int, Nt/10.0))
#DataSource = zeros(Nt)
#DataCoupling = zeros(Nt)

for k in 1:Q
    lattice.fₑ[:,:,k] = w[k]*lattice.Tₑ
    lattice.fₗ[:,:,k] = w[k]*lattice.Tₗ   
end

@profile for k in 1:Nt
    t = timeToFemto[k]
    # Save Data
    #ElectronTemperature[1,1, k] = lattice.Tₑ[1,1]
    #ElectronTemperature[end,1, k] = lattice.Tₑ[end,1]
    if k % 10 == 0
        ElectronTemperature[:,:,trunc(Int, k/10.0)] = lattice.Tₑ[:,:]
        LatticeTemperature[:,:,trunc(Int, k/10.0)] = lattice.Tₗ[:,:]
    end
    #SourceLaserEnergy[:,:,k] = lattice.S[:,:]
    #DataSource[k] = sum(lattice.Source[25,1,:])
    #DataCoupling[k] = sum(lattice.Coupling[25,1,:])

    calculateLaserEnergy2DGaussianAblation!(t, lattice)
    calculateG!(lattice)

    #ELECTRON
    CollisionElectronLB!(lattice, thermal)
    streamingElectron!(lattice)
    updateElectron!(lattice)
    boundaryElectronAblation!(lattice)

    #LATTICE
    CollisionLatticeLB!(lattice, thermal)
    #streamingLattice!(lattice)
    updateLattice!(lattice)
    #boundaryLatticeAblation!(lattice)

    
    if k % 10 == 0
        println("Step $k completed out of $Nt | $(k/Nt*100)%")
        calculateAblation!(lattice, LatticeTemperature, k)
        AblationData[:,:,trunc(Int, k/10)] = lattice.Ablation
    end
end








