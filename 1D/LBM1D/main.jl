using StaticArrays, Plots
#using Profile

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
lattice = Lattice()
thermal = ThermalProperties()

#SourceLaserEnergy = zeros(Nx, Nt)  
ElectronTemperature = zeros(Nx, Nt)
LatticeTemperature = zeros(Nx, Nt)

for k in 1:Q
    lattice.fₑ[:,k] = w[k]*lattice.Tₑ
    lattice.fₗ[:,k] = w[k]*lattice.Tₗ   
end

@time for k in 1:Nt
    t = timeToFemto[k]
    # Save Data
    ElectronTemperature[:, k] = lattice.Tₑ
    LatticeTemperature[:,k] = lattice.Tₗ
    #SourceLaserEnergy[:,k] = lattice.S

    calculateLaserEnergy2DGaussianAblation!(t, lattice)
    calculateG!(lattice)

    #ELECTRON
    CollisionElectronLB!(lattice, thermal)
    StreamingElectron!(lattice)
    UpdateElectron!(lattice)
    BoundaryElectron!(lattice)

    #LATTICE
    CollisionLatticeLB!(lattice, thermal)
    StreamingLattice!(lattice)
    UpdateLattice!(lattice)
    BoundaryLattice!(lattice)

    #calculateAblation!(lattice, ElectronTemperature, k)
    if k % 10 == 0
        println("Step $k completed out of $Nt | $(k/Nt*100)%")
        #AblationData[:,:,trunc(Int, k/10)] = lattice.Ablation
    end

end
#println(maximum(LatticeTemperature[:,:]))
