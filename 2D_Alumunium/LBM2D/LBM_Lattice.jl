module lbLattice

using ..Structs
using ..Constants
using ..Thermal
using ..Coupling

using Einsum

export CollisionLatticeLB!
export streamingLattice!
export updateLattice!
export boundaryLattice!
export boundaryLatticeAblation!

function CollisionLatticeLB!(lb::Lattice, t::ThermalProperties)
    lb.fₗ = lb.fₗ[:,:,:] + calculateLbCoupling(lb, t)
end

function calculateLbCoupling(lb::Lattice, t::ThermalProperties)
    calculateG!(lb)
    Threads.@threads for k in 1:Q
        for i in 1:Nx
            for j in 1:Ny
                lb.Couplingl[i,j,k] = dt*lb.G[i,j]*(lb.fₑ[i,j,k]-lb.fₗ[i,j,k])/t.Cₗ[i,j]
            end
        end
    end
    return lb.Couplingl
end

function updateLattice!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            lb.Tₗ[i,j] = 0
            for k in 1:Q
                lb.Tₗ[i,j] += lb.fₗ[i, j, k]
            end
        end
    end
end

end