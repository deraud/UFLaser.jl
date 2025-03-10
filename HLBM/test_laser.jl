using StaticArrays, Einsum, Plots
using Revise

include("./struct.jl")
include("./constants.jl")
include("./laser_source.jl")
include("./material_thermal.jl")
include("./coupling.jl")
include("./hyperbolic_collision.jl")

using .Constants
using .Structs
using .LaserEnergy
using .Thermal
using .Coupling
using .Collision

# Initialize Lattice
@time lattice = Lattice()
@time thermal = ThermalProperties()


# Simulate energy distribution over time
SourceLaserEnergy = zeros(Nx, Ny, Nt)  
ElectronTemperature = zeros(Nx, Ny, Nt)
AblationData = zeros(Nx,Ny,trunc(Int, Nt/10.0))
DataSource = zeros(Nt)
DataCoupling = zeros(Nt)

for k in 1:Q
    lattice.f[:,:,k] = w[k]*lattice.Tₑ
    lattice.fₗ[:,:,k] = w[k]*lattice.Tₗ   
end

@time for k in 1:Nt


    t = timeToFemto[k]
    SourceLaserEnergy[:,:,k] = lattice.S[:,:]
    calculateLaserEnergy2DGaussianAblation!(t, lattice)

end







