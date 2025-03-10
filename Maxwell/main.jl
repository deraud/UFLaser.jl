using Plots

include("././struct.jl")
include("././constants.jl")
include("././laser_source.jl")
include("././optics.jl")

using .Constants
using .Structs
using .LaserEnergy
using .OpticsFunctions

lattice = Lattice()
optics = Optics()

SourceLaserEnergy = zeros(Nx, Ny, Nt) 

for k in 1:200
    t = timeToFemto[k]
    calculateIndexK!(lattice, optics)
    calculateIndexN!(lattice, optics)
    calculateS!(lattice, optics, t)
    SourceLaserEnergy[:,:,k] = optics.S[:,:]
    if k % 10 == 0
        println("Step $k completed out of $Nt | $(k/Nt*100)%")
    end
end